/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/global.h>
#include <aspect/postprocess/particle_lpo.h>
#include <aspect/particle/property/lpo.h>
#include <aspect/utilities.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <stdio.h>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    LPO<dim>::LPO ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN()),
      output_file_number (numbers::invalid_unsigned_int),
      group_files(0),
      write_in_background_thread(false)
    {}

    template <int dim>
    LPO<dim>::~LPO ()
    {
      // make sure a thread that may still be running in the background,
      // writing data, finishes
      background_thread_master.join ();
      background_thread_content_raw.join ();
      background_thread_content_draw_volume_weighting.join ();
    }

    template <int dim>
    void
    LPO<dim>::initialize ()
    {
      const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
      this->random_number_generator.seed(random_number_seed+my_rank);
      // todo: check wheter this works correctly. Since the get_random_number function takes a reference
      // to the random_number_generator function, changing the function should mean that I have to update the
      // get_random_number function as well. But I will need to test this.
    }

    template <int dim>
    std::list<std::string>
    LPO<dim>::required_other_postprocessors () const
    {
      return {"particles"};
    }


    template <int dim>
    // We need to pass the arguments by value, as this function can be called on a separate thread:
    void LPO<dim>::writer (const std::string filename, //NOLINT(performance-unnecessary-value-param)
                           const std::string temporary_output_location, //NOLINT(performance-unnecessary-value-param)
                           const std::string *file_contents,
                           const bool compress_contents)
    {
      std::string tmp_filename = filename;
      if (temporary_output_location != "")
        {
          tmp_filename = temporary_output_location + "/aspect.tmp.XXXXXX";

          // Create the temporary file; get at the actual filename
          // by using a C-style string that mkstemp will then overwrite
          std::vector<char> tmp_filename_x (tmp_filename.size()+1);
          std::strcpy(tmp_filename_x.data(), tmp_filename.c_str());
          const int tmp_file_desc = mkstemp(tmp_filename_x.data());
          tmp_filename = tmp_filename_x.data();

          // If we failed to create the temp file, just write directly to the target file.
          // We also provide a warning about this fact. There are places where
          // this fails *on every node*, so we will get a lot of warning messages
          // into the output; in these cases, just writing multiple pieces to
          // std::cerr will produce an unreadable mass of text; rather, first
          // assemble the error message completely, and then output it atomically
          if (tmp_file_desc == -1)
            {
              const std::string x = ("***** WARNING: could not create temporary file <"
                                     +
                                     tmp_filename
                                     +
                                     ">, will output directly to final location. This may negatively "
                                     "affect performance. (On processor "
                                     + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                                     + ".)\n");

              std::cerr << x << std::flush;

              tmp_filename = filename;
            }
          else
            close(tmp_file_desc);
        }

      std::ofstream out(tmp_filename.c_str(), std::ofstream::binary);

      AssertThrow (out, ExcMessage(std::string("Trying to write to file <") +
                                   filename +
                                   ">, but the file can't be opened!"))

      // now write and then move the tmp file to its final destination
      // if necessary
      if (compress_contents)
        {
          namespace bio = boost::iostreams;

          std::stringstream compressed;
          std::stringstream origin(*file_contents);

          bio::filtering_streambuf<bio::input> compress_out;
          compress_out.push(bio::zlib_compressor());
          compress_out.push(origin);
          bio::copy(compress_out, out);
        }
      else
        {
          out << *file_contents;
        }

      out.close ();

      if (tmp_filename != filename)
        {
          std::string command = std::string("mv ") + tmp_filename + " " + filename;
          int error = system(command.c_str());

          AssertThrow(error == 0,
                      ExcMessage("Could not move " + tmp_filename + " to "
                                 + filename + ". On processor "
                                 + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + "."));
        }

      // destroy the pointer to the data we needed to write
      delete file_contents;
    }

    template <int dim>
    std::pair<std::string,std::string>
    LPO<dim>::execute (TableHandler &statistics)
    {

      const Particle::Property::Manager<dim> &manager = this->get_particle_world().get_property_manager();

      bool hexagonal_plugin_exists = manager.plugin_name_exists("decompose elastic matrix");

      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - output_interval;

      // If it's not time to generate an output file or we do not write output
      // return early.
      if (this->get_time() < last_output_time + output_interval && this->get_time() != end_time)
        return std::make_pair("","");

      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else
        ++output_file_number;

      // Now prepare everything for writing the output and choose output format
      std::string particle_file_prefix_master = this->get_output_directory() +  "particle_LPO/particles-" + Utilities::int_to_string (output_file_number, 5);
      std::string particle_file_prefix_content_raw = this->get_output_directory() +  "particle_LPO/LPO-" + Utilities::int_to_string (output_file_number, 5);
      std::string particle_file_prefix_content_draw_volume_weighting = this->get_output_directory() +  "particle_LPO/weighted_LPO-" + Utilities::int_to_string (output_file_number, 5);

      const auto &particle_handler = this->get_particle_world().get_particle_handler();

      std::stringstream string_stream_master;
      std::stringstream string_stream_content_raw;
      std::stringstream string_stream_content_draw_volume_weighting;

      string_stream_master << "id x y" << (dim == 3 ? " z" : " ") << " olivine_deformation_type"
                           << (hexagonal_plugin_exists ? (std::string(" full_norm_square ")
                                                          + "triclinic_norm_square_p1 triclinic_norm_square_p2 triclinic_norm_square_p3 "
                                                          + "monoclinic_norm_square_p1 monoclinic_norm_square_p2 monoclinic_norm_square_p3 "
                                                          + "orthohombic_norm_square_p1 orthohombic_norm_square_p2 orthohombic_norm_square_p3 "
                                                          + "tetragonal_norm_square_p1 tetragonal_norm_square_p2 tetragonal_norm_square_p3 "
                                                          + "hexagonal_norm_square_p1 hexagonal_norm_square_p2 hexagonal_norm_square_p3 "
                                                          + "isotropic_norm_square") : "") << std::endl;

      // get particle data
      bool wrote_weighted_header = false;
      bool wrote_unweighted_header = false;
      for (typename dealii::Particles::ParticleHandler<dim>::particle_iterator it = particle_handler.begin(); it != particle_handler.end(); ++it)
        {

          AssertThrow(it->has_properties(),
                      ExcMessage("No particle properties found. Make sure that the LPO particle property plugin is selected."));



          unsigned int id = it->get_id();
          const ArrayView<double> &properties = it->get_properties();

          const Particle::Property::ParticlePropertyInformation &property_information = this->get_particle_world().get_property_manager().get_data_info();

          AssertThrow(property_information.fieldname_exists("cpo mineral 0 type") ,
                      ExcMessage("No LPO particle properties found. Make sure that the LPO particle property plugin is selected."));



          const unsigned int lpo_data_position = property_information.n_fields() == 0
                                                 ?
                                                 0
                                                 :
                                                 property_information.get_position_by_field_name("cpo mineral 0 type");

          Point<dim> position = it->get_location();

          std::vector<unsigned int> deformation_type;
          std::vector<double> volume_fraction_mineral;
          std::vector<std::vector<double>> volume_fractions_grains;
          std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains;

          Particle::Property::LPO<dim>::load_particle_data(lpo_data_position,
                                                           properties,
                                                           deformation_type,
                                                           volume_fraction_mineral,
                                                           volume_fractions_grains,
                                                           a_cosine_matrices_grains);

          const unsigned int lpo_hex_data_position = property_information.n_fields() == 0 || hexagonal_plugin_exists == false
                                                     ?
                                                     0
                                                     :
                                                     property_information.get_position_by_field_name("lpo elastic axis e1");

          // write master file
          string_stream_master << id << " " << position << " " << properties[lpo_data_position];
          if (hexagonal_plugin_exists == true)
            {
              string_stream_master << " " << properties[lpo_hex_data_position+12] << " " << properties[lpo_hex_data_position+13]
                                   << " " << properties[lpo_hex_data_position+14] << " " << properties[lpo_hex_data_position+15]
                                   << " " << properties[lpo_hex_data_position+16] << " " << properties[lpo_hex_data_position+17]
                                   << " " << properties[lpo_hex_data_position+18] << " " << properties[lpo_hex_data_position+19]
                                   << " " << properties[lpo_hex_data_position+20] << " " << properties[lpo_hex_data_position+21]
                                   << " " << properties[lpo_hex_data_position+22] << " " << properties[lpo_hex_data_position+23]
                                   << " " << properties[lpo_hex_data_position+24] << " " << properties[lpo_hex_data_position+25]
                                   << " " << properties[lpo_hex_data_position+26] << " " << properties[lpo_hex_data_position+27]
                                   << " " << properties[lpo_hex_data_position+28];
            }
          string_stream_master <<  std::endl;


          // write content file
          std::vector<std::vector<std::vector<double> > > euler_angles;
          if (compute_raw_euler_angles == true)
            {
              euler_angles.resize(n_minerals);
              for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
                {

                  euler_angles[mineral_i].resize(n_grains);
                  for (unsigned int i_grain = 0; i_grain < n_grains; i_grain++)
                    {
                      euler_angles[mineral_i][i_grain] = euler_angles_from_rotation_matrix(a_cosine_matrices_grains[mineral_i][i_grain]);
                    }
                }
            }

          if (write_raw_lpo.size() != 0)
            {
              // write unweighted header
              if (wrote_unweighted_header == false)
                {
                  string_stream_content_raw << "id" << " " << std::setprecision(12);
                  for (unsigned int property_i = 0; property_i < write_raw_lpo.size(); ++property_i)
                    {
                      switch (write_raw_lpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_raw << "mineral_" << write_raw_lpo[property_i].first << "_volume_fraction" << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_raw << "mineral_" << write_raw_lpo[property_i].first << "_RM_0" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_1" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_2" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_3" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_4" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_5" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_6" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_7" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_RM_8" << " ";
                            break;

                          case Output::EulerAngles:
                            string_stream_content_raw << "mineral_" << write_raw_lpo[property_i].first << "_EA_phi" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_EA_theta" << " "
                                                      << "mineral_" << write_raw_lpo[property_i].first << "_EA_z" << " ";
                            break;

                          default:
                            Assert(false, ExcMessage("Internal error: raw LPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_raw << std::endl;
                  wrote_unweighted_header = true;
                }

              // write unweighted data
              for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                {
                  string_stream_content_raw << id << " ";
                  for (unsigned int property_i = 0; property_i < write_raw_lpo.size(); ++property_i)
                    {
                      switch (write_raw_lpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_raw << volume_fractions_grains[write_raw_lpo[property_i].first][grain_i] << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_raw << a_cosine_matrices_grains[write_raw_lpo[property_i].first][grain_i] << " ";
                            break;

                          case Output::EulerAngles:
                            Assert(compute_raw_euler_angles == true,
                                   ExcMessage("Internal error: writing out raw Euler angles, without them being computed."));
                            string_stream_content_raw << euler_angles[write_raw_lpo[property_i].first][grain_i][0] << " "
                                                      <<  euler_angles[write_raw_lpo[property_i].first][grain_i][1] << " "
                                                      <<  euler_angles[write_raw_lpo[property_i].first][grain_i][2] << " ";
                            break;

                          default:
                            Assert(false, ExcMessage("Internal error: raw LPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_raw << std::endl;
                }

            }
          if (write_draw_volume_weighted_lpo.size() != 0)
            {
              std::vector<std::vector<std::vector<double> >> weighted_euler_angles;

              weighted_euler_angles.resize(n_minerals);
              for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
                {
                  weighted_euler_angles[mineral_i] = random_draw_volume_weighting(volume_fractions_grains[mineral_i], euler_angles[mineral_i]);

                  Assert(weighted_euler_angles[mineral_i].size() == euler_angles[mineral_i].size(),
                         ExcMessage("Weighted angles vector (size = " + std::to_string(weighted_euler_angles[mineral_i].size()) +
                                    ") has different size from input angles (size = " + std::to_string(euler_angles[mineral_i].size()) + ")."));

                }


              std::vector<std::vector<Tensor<2,3> > > weighted_a_cosine_matrices;

              if (compute_weighted_A_matrix == true)
                {
                  weighted_a_cosine_matrices.resize(n_minerals);
                  for (unsigned int mineral_i = 0; mineral_i < n_minerals; ++mineral_i)
                    {
                      weighted_a_cosine_matrices[mineral_i].resize(weighted_euler_angles[mineral_i].size());

                      for (unsigned int i = 0; i < weighted_euler_angles.size(); ++i)
                        {
                          weighted_a_cosine_matrices[mineral_i][i] = euler_angles_to_rotation_matrix(weighted_euler_angles[mineral_i][i][0],
                                                                                                     weighted_euler_angles[mineral_i][i][1],
                                                                                                     weighted_euler_angles[mineral_i][i][2]);
                        }
                    }
                }

              // write weighted header
              if (wrote_weighted_header == false)
                {
                  string_stream_content_draw_volume_weighting << "id" << " ";
                  for (unsigned int property_i = 0; property_i < write_draw_volume_weighted_lpo.size(); ++property_i)
                    {
                      switch (write_draw_volume_weighted_lpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_volume_fraction" << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_0" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_1" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_2" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_3" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_4" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_5" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_6" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_7" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_RM_8" << " ";
                            break;

                          case Output::EulerAngles:
                            string_stream_content_draw_volume_weighting << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_EA_phi" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_EA_theta" << " "
                                                                        << "mineral_" << write_draw_volume_weighted_lpo[property_i].first << "_EA_z" << " ";
                            break;

                          default:
                            Assert(false, ExcMessage("Internal error: raw LPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_draw_volume_weighting << std::endl;
                  wrote_weighted_header = true;
                }

              // write data
              for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
                {
                  string_stream_content_draw_volume_weighting << id << " ";
                  for (unsigned int property_i = 0; property_i < write_draw_volume_weighted_lpo.size(); ++property_i)
                    {
                      switch (write_draw_volume_weighted_lpo[property_i].second)
                        {
                          case Output::VolumeFraction:
                            string_stream_content_draw_volume_weighting << volume_fractions_grains[write_draw_volume_weighted_lpo[property_i].first][grain_i] << " ";
                            break;

                          case Output::RotationMatrix:
                            string_stream_content_draw_volume_weighting << weighted_a_cosine_matrices[write_draw_volume_weighted_lpo[property_i].first][grain_i] << " ";
                            break;

                          case Output::EulerAngles:
                            Assert(compute_raw_euler_angles == true,
                                   ExcMessage("Internal error: writing out raw Euler angles, without them being computed."));
                            string_stream_content_draw_volume_weighting << weighted_euler_angles[write_draw_volume_weighted_lpo[property_i].first][grain_i][0] << " "
                                                                        <<  weighted_euler_angles[write_draw_volume_weighted_lpo[property_i].first][grain_i][1] << " "
                                                                        <<  weighted_euler_angles[write_draw_volume_weighted_lpo[property_i].first][grain_i][2] << " ";
                            break;
                          default:
                            Assert(false, ExcMessage("Internal error: raw LPO postprocess case not found."));
                            break;
                        }
                    }
                  string_stream_content_draw_volume_weighting << std::endl;
                }
            }
        }

      std::string filename_master = particle_file_prefix_master + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD),4) + ".dat";
      std::string filename_raw = particle_file_prefix_content_raw + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD),4) + ".dat";
      std::string filename_draw_volume_weighting = particle_file_prefix_content_draw_volume_weighting + "." + Utilities::int_to_string(dealii::Utilities::MPI::this_mpi_process (MPI_COMM_WORLD),4) + ".dat";

      std::string *file_contents_master = new std::string (string_stream_master.str());
      std::string *file_contents_raw = new std::string (string_stream_content_raw.str());
      std::string *file_contents_draw_volume_weighting = new std::string (string_stream_content_draw_volume_weighting.str());

      if (write_in_background_thread)
        {
          // Wait for all previous write operations to finish, should
          // any be still active,
          background_thread_master.join ();

          // then continue with writing the master file
          background_thread_master = Threads::new_thread (&writer,
                                                          filename_master,
                                                          temporary_output_location,
                                                          file_contents_master,
                                                          false);

          if (write_raw_lpo.size() != 0)
            {
              // Wait for all previous write operations to finish, should
              // any be still active,
              background_thread_content_raw.join ();

              // then continue with writing our own data.
              background_thread_content_raw = Threads::new_thread (&writer,
                                                                   filename_raw,
                                                                   temporary_output_location,
                                                                   file_contents_raw,
                                                                   compress_lpo_data_files);
            }

          if (write_draw_volume_weighted_lpo.size() != 0)
            {
              // Wait for all previous write operations to finish, should
              // any be still active,
              background_thread_content_draw_volume_weighting.join ();

              // then continue with writing our own data.
              background_thread_content_draw_volume_weighting = Threads::new_thread (&writer,
                                                                                     filename_draw_volume_weighting,
                                                                                     temporary_output_location,
                                                                                     file_contents_draw_volume_weighting,
                                                                                     compress_lpo_data_files);
            }
        }
      else
        {
          writer(filename_master,temporary_output_location,file_contents_master, false);
          if (write_raw_lpo.size() != 0)
            writer(filename_raw,temporary_output_location,file_contents_raw, compress_lpo_data_files);
          if (write_draw_volume_weighted_lpo.size() != 0)
            writer(filename_draw_volume_weighting,temporary_output_location,file_contents_draw_volume_weighting, compress_lpo_data_files);
        }


      // up the next time we need output
      set_last_output_time (this->get_time());

      const std::string particle_lpo_output = particle_file_prefix_content_raw;

      // record the file base file name in the output file
      statistics.add_value ("Particle LPO file name",
                            particle_lpo_output);
      return std::make_pair("Writing particle lpo output:", particle_lpo_output);
    }

    template<int dim>
    std::vector<std::vector<double> >
    LPO<dim>::random_draw_volume_weighting(std::vector<double> fv,
                                           std::vector<std::vector<double>> angles) const
    {
      // Get volume weighted euler angles, using random draws to convert odf
      // to a discrete number of orientations, weighted by volume
      // 1a. Get index that would sort volume fractions AND
      //ix = np.argsort(fv[q,:]);
      // 1b. Get the sorted volume and angle arrays
      std::vector<double> fv_to_sort = fv;
      std::vector<double> fv_sorted = fv;
      std::vector<std::vector<double>> angles_sorted = angles;

      unsigned int n_grain = fv_to_sort.size();


      /**
       * ...
       */
      for (int i = n_grain-1; i >= 0; --i)
        {
          unsigned int index_max_fv = std::distance(fv_to_sort.begin(),max_element(fv_to_sort.begin(), fv_to_sort.end()));

          fv_sorted[i] = fv_to_sort[index_max_fv];
          angles_sorted[i] = angles[index_max_fv];
          Assert(angles[index_max_fv].size() == 3, ExcMessage("angles vector (size = " + std::to_string(angles[index_max_fv].size()) +
                                                              ") should have size 3."));
          Assert(angles_sorted[i].size() == 3, ExcMessage("angles_sorted vector (size = " + std::to_string(angles_sorted[i].size()) +
                                                          ") should have size 3."));
          fv_to_sort[index_max_fv] = -1;
        }


      // 2. Get cumulative weight for volume fraction
      std::vector<double> cum_weight(n_grains);
      std::partial_sum(fv_sorted.begin(),fv_sorted.end(),cum_weight.begin());
      // 3. Generate random indices
      boost::random::uniform_real_distribution<> dist(0, 1);
      std::vector<double> idxgrain(n_grains);
      for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
        idxgrain[grain_i] = dist(this->random_number_generator); //random.rand(ngrains,1);

      // 4. Find the maximum cum_weight that is less than the random value.
      // the euler angle index is +1. For example, if the idxGrain(g) < cumWeight(1),
      // the index should be 1 not zero)
      std::vector<std::vector<double>> angles_out(n_grains,std::vector<double>(3));
      for (unsigned int grain_i = 0; grain_i < n_grains; ++grain_i)
        {
          unsigned int counter = 0;
          for (unsigned int grain_j = 0; grain_j < n_grains; ++grain_j)
            {
              if (cum_weight[grain_j] < idxgrain[grain_i])
                {
                  counter++;
                }
              else
                {
                  break;
                }
            }
          angles_out[grain_i] = angles_sorted[counter];
        }
      return angles_out;
    }


    template <int dim>
    double
    LPO<dim>::wrap_angle(const double angle) const
    {
      return angle - 360.0*std::floor(angle/360.0);
    }

    template <int dim>
    std::vector<double>
    LPO<dim>::euler_angles_from_rotation_matrix(const Tensor<2,3> &rotation_matrix) const
    {
      // ZXZ Euler angles
      std::vector<double> euler_angles(3);
      for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 3; j++)
          Assert(abs(rotation_matrix[i][j]) <= 1.0,
                 ExcMessage("rotation_matrix[" + std::to_string(i) + "][" + std::to_string(j) +
                            "] is larger than one: " + std::to_string(rotation_matrix[i][j]) + " (" + std::to_string(rotation_matrix[i][j]-1.0) + "). rotation_matrix = \n"
                            + std::to_string(rotation_matrix[0][0]) + " " + std::to_string(rotation_matrix[0][1]) + " " + std::to_string(rotation_matrix[0][2]) + "\n"
                            + std::to_string(rotation_matrix[1][0]) + " " + std::to_string(rotation_matrix[1][1]) + " " + std::to_string(rotation_matrix[1][2]) + "\n"
                            + std::to_string(rotation_matrix[2][0]) + " " + std::to_string(rotation_matrix[2][1]) + " " + std::to_string(rotation_matrix[2][2])));


      AssertThrow(rotation_matrix[2][2] <= 1.0, ExcMessage("rot_matrix[2][2] > 1.0"));
      //std::cout<<"The rotation matrix is: "<<rotation_matrix << std::endl;
      /*for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  std::cout << rotation_matrix[i][j] << ", ";
                }
                std::cout << std::endl;
              }*/

      const double theta = std::acos(rotation_matrix[2][2]);//Mod Agi 2022.06.22
      double phi1 = 0.0;
      double phi2 = 0.0;

      if (theta != 0.0 && theta != dealii::numbers::PI)
        {
          //
          phi1  = std::atan2(rotation_matrix[2][0]/sin(theta),rotation_matrix[2][1]/-sin(theta));//Mod Agi 2022.06.22
          phi2  = std::atan2(rotation_matrix[0][2]/sin(theta),rotation_matrix[1][2]/sin(theta));//Mod Agi 2022.06.22
        }
      else
        {
          // note that in the case theta is 0 or phi a dimension is lost
          // see: https://en.wikipedia.org/wiki/Gimbal_lock. We set phi1
          // to 0 and compute the corresponding phi2. The resulting direction
          // (cosine matrix) should be the same.
          if (theta == 0.0)
            {
              phi2 = - phi1 - std::atan2(-rotation_matrix[0][1],rotation_matrix[0][0]);// Mod Yijun
            }
          else
            {
              phi2 = phi1 + std::atan2(rotation_matrix[0][1],rotation_matrix[0][0]);
            }

        }



      AssertThrow(!std::isnan(phi1), ExcMessage(" phi1 is nan. theta = " + std::to_string(theta) + ", rotation_matrix[2][2]= " + std::to_string(rotation_matrix[2][2])
                                                + ", acos(rotation_matrix[2][2]) = " + std::to_string(std::acos(rotation_matrix[2][2])) + ", acos(1.0) = " + std::to_string(std::acos(1.0))));
      AssertThrow(!std::isnan(theta), ExcMessage(" theta is nan."));
      AssertThrow(!std::isnan(phi2), ExcMessage(" phi2 is nan."));

      euler_angles[0] = wrap_angle(phi1 * rad_to_degree);
      euler_angles[1] = wrap_angle(theta * rad_to_degree);
      euler_angles[2] = wrap_angle(phi2 * rad_to_degree);
      //std::cout<<"phi1 is: "<<euler_angles[0]<<std::endl;
      //std::cout<<"theta is: "<<euler_angles[1]<<std::endl;
      //std::cout<<"phi2 is: "<<euler_angles[2]<<std::endl;


      AssertThrow(!std::isnan(euler_angles[0]), ExcMessage(" euler_angles[0] is nan."));
      AssertThrow(!std::isnan(euler_angles[1]), ExcMessage(" euler_angles[1] is nan."));
      AssertThrow(!std::isnan(euler_angles[2]), ExcMessage(" euler_angles[2] is nan."));

      return euler_angles;
    }

    template <int dim>
    Tensor<2,3>
    LPO<dim>::euler_angles_to_rotation_matrix(double phi1_d, double theta_d, double phi2_d) const
    {
      // ZXZ Euler angles
      const double phi1 = phi1_d *degree_to_rad;
      const double theta = theta_d *degree_to_rad;
      const double phi2 = phi2_d *degree_to_rad;
      Tensor<2,3> rot_matrix;


      rot_matrix[0][0] = cos(phi2)*cos(phi1) - cos(theta)*sin(phi1)*sin(phi2);
      rot_matrix[0][1] = -cos(phi2)*sin(phi1) - cos(theta)*cos(phi1)*sin(phi2);
      rot_matrix[0][2] = sin(phi2)*sin(theta); //Mod Agi 2022.06.23.

      rot_matrix[1][0] = sin(phi2)*cos(phi1) + cos(theta)*sin(phi1)*cos(phi2);
      rot_matrix[1][1] = -sin(phi2)*sin(phi1) + cos(theta)*cos(phi1)*cos(phi2);
      rot_matrix[1][2] = -cos(phi2)*sin(theta);//Mod Agi 2022.06.23.

      rot_matrix[2][0] = sin(theta)*sin(phi1);//Mod Agi 2022.06.23.
      rot_matrix[2][1] = sin(theta)*cos(phi1);//Mod Agi 2022.06.23.
      rot_matrix[2][2] = cos(theta);
      AssertThrow(rot_matrix[2][2] <= 1.0, ExcMessage("rot_matrix[2][2] > 1.0"));
      return rot_matrix;
    }



    template <int dim>
    void
    LPO<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }

    template <int dim>
    typename LPO<dim>::Output
    LPO<dim>::string_to_output_enum(std::string string)
    {
      //olivine volume fraction, olivine A matrix, olivine Euler angles, enstatite volume fraction, enstatite A matrix, enstatite Euler angles
      if (string == "volume fraction")
        return Output::VolumeFraction;
      if (string == "rotation matrix")
        return Output::RotationMatrix;
      if (string == "Euler angles")
        return Output::EulerAngles;
      return Output::not_found;
    }


    template <int dim>
    void
    LPO<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.enter_subsection("LPO");
          {
            prm.declare_entry ("Number of grains per praticle", "50",
                               Patterns::Integer (0),
                               "The number of grains of olivine and the number of grain of enstatite "
                               "each particle contains.");

            prm.enter_subsection("D-Rex 2004");
            {
              prm.declare_entry ("Minerals", "Olivine: Karato 2008, Enstatite",
                                 Patterns::List(Patterns::Anything()),
                                 "This determines what minerals and fabrics or fabric selectors are used used for the LPO calculation. "
                                 "The options are Olivine: A-fabric, Olivine: B-fabric, Olivine: C-fabric, Olivine: D-fabric, "
                                 "Olivine: E-fabric, Olivine: Karato 2008 or Enstatite. The "
                                 "Karato 2008 selector selects a fabric based on stress and water content as defined in "
                                 "figure 4 of the Karato 2008 review paper (doi: 10.1146/annurev.earth.36.031207.124120).");
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
        prm.enter_subsection("LPO");
        {
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          prm.declare_entry ("Random number seed", "1",
                             Patterns::Integer (0),
                             "The seed used to generate random numbers. This will make sure that "
                             "results are reproducable as long as the problem is run with the "
                             "same amount of MPI processes. It is implemented as final seed = "
                             "user seed + MPI Rank. ");

          prm.declare_entry ("Write in background thread", "false",
                             Patterns::Bool(),
                             "File operations can potentially take a long time, blocking the "
                             "progress of the rest of the model run. Setting this variable to "
                             "`true' moves this process into a background threads, while the "
                             "rest of the model continues.");

          prm.declare_entry ("Temporary output location", "",
                             Patterns::Anything(),
                             "On large clusters it can be advantageous to first write the "
                             "output to a temporary file on a local file system and later "
                             "move this file to a network file system. If this variable is "
                             "set to a non-empty string it will be interpreted as a "
                             "temporary storage location.");

          prm.declare_entry ("Write out raw lpo data",
                             "olivine volume fraction,olivine Euler angles,enstatite volume fraction,enstatite Euler angles",
                             Patterns::List(Patterns::Anything()),
                             "A list containing the what part of the particle lpo data needs "
                             "to be written out after the particle id. This writes out the raw "
                             "lpo data files for each MPI process. It can write out the following data: "
                             "olivine volume fraction, olivine A matrix, olivine Euler angles, "
                             "enstatite volume fraction, enstatite A matrix, enstatite Euler angles. \n"
                             "Note that the A matrix and Euler angles both contain the same "
                             "information, but in a different format. Euler angles are recommended "
                             "over the A matrix since they only require to write 3 values instead "
                             "of 9. If the list is empty, this file will not be written."
                             "Furthermore, the entries will be written out in the order given, "
                             "and if entries are entered muliple times, they will be written "
                             "out multiple times.");

          prm.declare_entry ("Write out draw volume weighted lpo data",
                             "olivine Euler angles,enstatite Euler angles",
                             Patterns::List(Patterns::Anything()),
                             "A list containing the what part of the random draw volume "
                             "weighted particle lpo data needs to be written out after "
                             "the particle id. after using a random draw volume weighting. "
                             "The random draw volume weigthing uses a uniform random distribution "
                             "This writes out the raw lpo data files for "
                             "each MPI process. It can write out the following data: "
                             "olivine volume fraction, olivine A matrix, olivine Euler angles, "
                             "enstatite volume fraction, enstatite A matrix, enstatite Euler angles. \n"
                             "Note that the A matrix and Euler angles both contain the same "
                             "information, but in a different format. Euler angles are recommended "
                             "over the A matrix since they only require to write 3 values instead "
                             "of 9. If the list is empty, this file will not be written. "
                             "Furthermore, the entries will be written out in the order given, "
                             "and if entries are entered muliple times, they will be written "
                             "out multiple times.");
          prm.declare_entry ("Compress lpo data files", "true",
                             Patterns::Bool(),
                             "Wether to compress the raw and weighted lpo data output files with zlib.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

    }


    template <int dim>
    void
    LPO<dim>::parse_parameters (ParameterHandler &prm)
    {
      end_time = prm.get_double ("End time");
      if (this->convert_output_to_years())
        end_time *= year_in_seconds;

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.enter_subsection("LPO");
          {
            n_grains = prm.get_integer("Number of grains per praticle");
            prm.enter_subsection("D-Rex 2004");
            {
              // Static variable of LPO has not been initialize yet, so we need to get it directly.
              const std::vector<std::string> temp_deformation_type_selector = dealii::Utilities::split_string_list(prm.get("Minerals")); // todo trim?
              n_minerals = temp_deformation_type_selector.size();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
        prm.enter_subsection("LPO");
        {
          output_interval = prm.get_double ("Time between data output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          random_number_seed = prm.get_integer ("Random number seed");

          //AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
          //            ExcMessage("Postprocessing nonlinear iterations in models with "
          //                       "particles is currently not supported."));

          aspect::Utilities::create_directory (this->get_output_directory() + "particle_LPO/",
                                               this->get_mpi_communicator(),
                                               true);

          write_in_background_thread = prm.get_bool("Write in background thread");
          temporary_output_location = prm.get("Temporary output location");

          if (temporary_output_location != "")
            {
              // Check if a command-processor is available by calling system() with a
              // null pointer. System is guaranteed to return non-zero if it finds
              // a terminal and zero if there is none (like on the compute nodes of
              // some cluster architectures, e.g. IBM BlueGene/Q)
              AssertThrow(system((char *)nullptr) != 0,
                          ExcMessage("Usage of a temporary storage location is only supported if "
                                     "there is a terminal available to move the files to their final location "
                                     "after writing. The system() command did not succeed in finding such a terminal."));
            }

          std::vector<std::string> write_raw_lpo_list = Utilities::split_string_list(prm.get("Write out raw lpo data"));
          write_raw_lpo.resize(write_raw_lpo_list.size());
          bool found_euler_angles = false;
          for (unsigned int i = 0; i < write_raw_lpo_list.size(); ++i)
            {
              std::vector<std::string> split_raw_lpo_instructions = Utilities::split_string_list(write_raw_lpo_list[i],':');

              AssertThrow(split_raw_lpo_instructions.size() == 2,
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option "
                                     + "because it should contain a mineral identification and a output specifier seprated by a colon (:). This entry "
                                     + "does not follow those rules."));

              // get mineral number
              std::vector<std::string> mineral_instructions = Utilities::split_string_list(split_raw_lpo_instructions[0],' ');
              AssertThrow(mineral_instructions.size() == 2,
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option "
                                     + "because the mineral identification part should contain two elements, the word mineral and a number."));

              AssertThrow(mineral_instructions[0] == "mineral",
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option "
                                     + "because the mineral identification part should start with the word mineral and it starts with \""
                                     + mineral_instructions[0] + "\"."));

              int mineral_number = Utilities::string_to_int(mineral_instructions[1]);

              AssertThrow((unsigned int) mineral_number < n_minerals,
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option "
                                     + "because the mineral number (" + std::to_string(mineral_number) + ") is larger than the number of minerals "
                                     + "provided in the LPO subsection (" + std::to_string(n_minerals) + ")."));

              // get mineral fabric/deformation type
              Output lpo_fabric_instruction = string_to_output_enum(split_raw_lpo_instructions[1]);

              AssertThrow(lpo_fabric_instruction != Output::not_found,
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option."))

              if (lpo_fabric_instruction == Output::EulerAngles)
                found_euler_angles = true;

              write_raw_lpo[i] = std::make_pair(mineral_number,lpo_fabric_instruction);
            }

          std::vector<std::string> write_draw_volume_weighted_lpo_list = Utilities::split_string_list(prm.get("Write out draw volume weighted lpo data"));
          write_draw_volume_weighted_lpo.resize(write_draw_volume_weighted_lpo_list.size());
          bool found_A_matrix = false;
          for (unsigned int i = 0; i < write_draw_volume_weighted_lpo_list.size(); ++i)
            {
              std::vector<std::string> split_draw_volume_weighted_lpo_instructions = Utilities::split_string_list(write_draw_volume_weighted_lpo_list[i],':');

              AssertThrow(split_draw_volume_weighted_lpo_instructions.size() == 2,
                          ExcMessage("Value \""+ write_draw_volume_weighted_lpo_list[i] +"\", set in \"Write out draw volume weighted lpo data\", is not a correct option "
                                     + "because it should contain a mineral identification and a output specifier seprated by a colon (:). This entry "
                                     + "does not follow those rules."));

              // get mineral number
              std::vector<std::string> mineral_instructions = Utilities::split_string_list(split_draw_volume_weighted_lpo_instructions[0],' ');
              AssertThrow(mineral_instructions.size() == 2,
                          ExcMessage("Value \""+ write_draw_volume_weighted_lpo_list[i] +"\", set in \"Write out draw volume weighted lpo data\", is not a correct option "
                                     + "because the mineral identification part should contain two elements, the word mineral and a number."));

              AssertThrow(mineral_instructions[0] == "mineral",
                          ExcMessage("Value \""+ write_draw_volume_weighted_lpo_list[i] +"\", set in \"Write out draw volume weighted lpo data\", is not a correct option "
                                     + "because the mineral identification part should start with the word mineral and it starts with \""
                                     + mineral_instructions[0] + "\"."));

              int mineral_number = Utilities::string_to_int(mineral_instructions[1]);
              AssertThrow((unsigned int) mineral_number < n_minerals,
                          ExcMessage("Value \""+ write_raw_lpo_list[i] +"\", set in \"Write out raw lpo data\", is not a correct option "
                                     + "because the mineral number (" + std::to_string(mineral_number) + ") is larger than the number of minerals "
                                     + "provided in the LPO subsection (" + std::to_string(n_minerals) + ")."));

              // get mineral fabric/deformation type
              Output lpo_fabric_instruction = string_to_output_enum(split_draw_volume_weighted_lpo_instructions[1]);

              AssertThrow(lpo_fabric_instruction != Output::not_found,
                          ExcMessage("Value \""+ write_draw_volume_weighted_lpo_list[i] +"\", set in \"Write out draw volume weighted lpo data\", is not a correct option."))

              if (lpo_fabric_instruction == Output::RotationMatrix)
                found_A_matrix = true;

              write_draw_volume_weighted_lpo[i] = std::make_pair(mineral_number,lpo_fabric_instruction);
            }

          /*std::vector<std::string> write_draw_volume_weighted_lpo_tmp = Utilities::split_string_list(prm.get("Write out draw volume weighted lpo data"));
          write_draw_volume_weighted_lpo.resize(write_draw_volume_weighted_lpo_tmp.size());
          bool found_A_matrix = false;
          for (unsigned int i = 0; i < write_draw_volume_weighted_lpo_tmp.size(); ++i)
            {
              write_draw_volume_weighted_lpo[i] = string_to_output_enum(write_draw_volume_weighted_lpo_tmp[i]);
              AssertThrow(write_draw_volume_weighted_lpo[i] != Output::not_found,
                          ExcMessage("Value \""+ write_draw_volume_weighted_lpo_tmp[i] +"\", set in \"Write out raw lpo data\", is not a correct option."));

              if (write_draw_volume_weighted_lpo[i] == Output::olivine_A_matrix || write_draw_volume_weighted_lpo[i] == Output::enstatite_A_matrix)
                found_A_matrix = true;
            }*/




          if (write_draw_volume_weighted_lpo_list.size() != 0 || found_euler_angles == true)
            compute_raw_euler_angles = true;
          else
            compute_raw_euler_angles = false;

          if (write_draw_volume_weighted_lpo_list.size() != 0 && found_A_matrix == true)
            compute_weighted_A_matrix = true;
          else
            compute_weighted_A_matrix = false;

          compress_lpo_data_files = prm.get_bool("Compress lpo data files");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(LPO,
                                  "lpo",
                                  "A Postprocessor that creates particles that follow the "
                                  "velocity field of the simulation. The particles can be generated "
                                  "and propagated in various ways and they can carry a number of "
                                  "constant or time-varying properties. The postprocessor can write "
                                  "output positions and properties of all particles at chosen intervals, "
                                  "although this is not mandatory. It also allows other parts of the "
                                  "code to query the particles for information.")
  }
}

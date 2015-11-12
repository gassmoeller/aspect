/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/output/hdf5.h>

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim>
      HDF5Output<dim>::HDF5Output()
        :
        file_index(0)
      {

#ifndef DEAL_II_WITH_HDF5
        AssertThrow (false,
                     ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                 "so HDF5 output is not possible. Please "
                                 "recompile deal.ii with HDF5 support turned on "
                                 "or select a different tracer output format."));
#endif

      }

      template <int dim>
      std::string
      HDF5Output<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                            const std::vector<std::pair<std::string, unsigned int> > &/*property_component_list*/,
                                            const double current_time)
      {
#ifdef DEAL_II_WITH_HDF5
        const std::string output_file_prefix = "particle-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix = this->get_output_directory() + output_file_prefix;
        const std::string h5_filename = output_path_prefix+".h5";

        // Create parallel file access
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        // Create property list for collective dataset write
        hid_t write_properties = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
        H5Pset_fapl_mpio(plist_id, this->get_mpi_communicator(), MPI_INFO_NULL);
        H5Pset_dxpl_mpio(write_properties, H5FD_MPIO_COLLECTIVE);
#endif

        // Create the file
        hid_t h5_file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        H5Pclose(plist_id);


        // Create the file dataspace descriptions
        types::particle_index n_local_particles = particles.size();
        types::particle_index n_global_particles = Utilities::MPI::sum(n_local_particles,this->get_mpi_communicator());

        hsize_t dims[2];
        dims[0] = n_global_particles;
        dims[1] = 3;
        hid_t one_dim_ds_id = H5Screate_simple(1, dims, NULL);
        hid_t dim_dataspace_id = H5Screate_simple(2, dims, NULL);

        // Create the datasets
#if H5Dcreate_vers == 1
        hid_t position_dataset = H5Dcreate(h5_file_id, "nodes", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT);
        hid_t particle_index_dataset = H5Dcreate(h5_file_id, "id", HDF5_PARTICLE_INDEX_TYPE, one_dim_ds_id, H5P_DEFAULT);
#else
        hid_t position_dataset = H5Dcreate(h5_file_id, "nodes", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        hid_t particle_index_dataset = H5Dcreate(h5_file_id, "id", HDF5_PARTICLE_INDEX_TYPE, one_dim_ds_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

        // Close the file dataspaces
        H5Sclose(dim_dataspace_id);
        H5Sclose(one_dim_ds_id);

        // Just in case we forget what they are
        hsize_t count[2];
        count[0] = 1;
        count[1] = 0;
        hid_t file_dataspace_id = H5Screate_simple (1, count, NULL);
#if H5Acreate_vers == 1
        hid_t pattr_id = H5Acreate(h5_file_id, "Ermahgerd! Pertecrs!", H5T_NATIVE_INT, file_dataspace_id, H5P_DEFAULT);
#else
        hid_t pattr_id = H5Acreate(h5_file_id, "Ermahgerd! Pertecrs!", H5T_NATIVE_INT, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
        H5Aclose(pattr_id);
        H5Sclose(file_dataspace_id);

        // Get the offset of the local particles among all processes
        types::particle_index local_particle_index_offset;
        MPI_Scan(&n_local_particles, &local_particle_index_offset, 1, ASPECT_TRACER_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());

        count[0] = n_local_particles;
        count[1] = 3;

        hsize_t offset[2];
        offset[0] = local_particle_index_offset - n_local_particles;
        offset[1] = 0;

        // Select the appropriate dataspace for this process
        hid_t one_dim_mem_ds_id = H5Screate_simple(1, count, NULL);
        hid_t dim_mem_ds_id = H5Screate_simple(2, count, NULL);
        hid_t pos_file_dataspace_id = H5Dget_space(position_dataset);
        hid_t pid_file_dataspace_id = H5Dget_space(particle_index_dataset);

        // And select the hyperslabs from each dataspace
        H5Sselect_hyperslab(pos_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Sselect_hyperslab(pid_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

        // Read the local particle data
        std::vector<double> position_data (3 * n_local_particles,0.0);
        std::vector<types::particle_index> index_data (n_local_particles);

        typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator it = particles.begin();
        for (unsigned int i = 0; it != particles.end(); ++i, ++it)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                position_data[i*3+d] = it->second.get_location()(d);
              }
            index_data[i] = it->second.get_id();
          }

        // Write particle data to the HDF5 file
        H5Dwrite(position_dataset, H5T_NATIVE_DOUBLE, dim_mem_ds_id, pos_file_dataspace_id, write_properties, &position_data[0]);
        H5Dwrite(particle_index_dataset, HDF5_PARTICLE_INDEX_TYPE, one_dim_mem_ds_id, pid_file_dataspace_id, write_properties, &index_data[0]);

        H5Pclose(write_properties);
        H5Sclose(one_dim_mem_ds_id);
        H5Sclose(dim_mem_ds_id);
        H5Sclose(pos_file_dataspace_id);
        H5Sclose(pid_file_dataspace_id);
        H5Dclose(position_dataset);
        H5Dclose(particle_index_dataset);
        H5Fclose(h5_file_id);

        // Record and output XDMF info on root process
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            const std::string local_h5_filename = output_file_prefix+".h5";
            XDMFEntry   entry(local_h5_filename, current_time, n_global_particles, 0, 3);
            DataOut<dim> data_out;
            const std::string xdmf_filename = (this->get_output_directory() + "particle.xdmf");

            entry.add_attribute("id", 1);

            xdmf_entries.push_back(entry);

            data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), this->get_mpi_communicator());
          }

        file_index++;

        return output_path_prefix;
#else
        (void) particles;
        (void) current_time;
        return "";
#endif
      }


      template <int dim>
      template <class Archive>
      void HDF5Output<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &file_index
           &xdmf_entries
        ;
      }

      template <int dim>
      void
      HDF5Output<dim>::save (std::ostringstream &os) const
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      template <int dim>
      void
      HDF5Output<dim>::load (std::istringstream &is)
      {
        aspect::iarchive ia (is);
        ia >> (*this);
      }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      ASPECT_REGISTER_PARTICLE_OUTPUT(HDF5Output,
                                      "hdf5",
                                      "This particle output plugin writes particle "
                                      "positions and properties into hdf5 files.")
    }
  }
}


/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/velocity_boundary_conditions/plume_box.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {

    namespace internal
    {
      template <int dim>
      BoxPlatesLookup<dim>::BoxPlatesLookup(const std::string &filename,
                                            const unsigned int components,
                                            const double time_scale_factor,
                                            const double velocity_scale_factor,
                                            const bool interpolate_velocity)
        :
        components(components),
        data(components),
        time_scale_factor(time_scale_factor),
        velocity_scale_factor(velocity_scale_factor),
        interpolate_velocity(interpolate_velocity)
      {
        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (Utilities::fexists(filename),
                     ExcMessage (std::string("Plate velocity file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") + filename));
        getline(in, temp); // eat first line


        double old_time (std::numeric_limits<double>::quiet_NaN());
        velocity_map velocity_slice;
        double start_time,end_time,vx,vy,omega;
        unsigned char plate_id;
        while (in >> start_time >> end_time >> plate_id >> vx >> vy >> omega)
          {
            // scale all properties
            vx *= velocity_scale_factor;
            vy *= velocity_scale_factor;
            // omega is in [velocity/km] scale it to [1/s]
            omega *= velocity_scale_factor / 1000;
            start_time *= time_scale_factor;
            end_time *= time_scale_factor;

            // If we have just read in the next time slice
            // save the last one and start a new one
            if (start_time > old_time + 1e-16)
              {
                velocity_values.push_back(std::make_pair(old_time,velocity_slice));
                velocity_slice = velocity_map();
              }

            // each time will contain a map from plate character to velocity
            plate_velocity velocity;
            switch (dim)
              {
                case 1:
                  velocity.first[0] = vx;
                  break;
                case 2:
                  velocity.first[0] = vx;
                  velocity.first[1] = vy;
                  break;
                default:
                  AssertThrow(false,ExcNotImplemented());
                  break;
              }

            velocity.second = omega;

            velocity_slice.insert(std::pair<unsigned char,plate_velocity >
                                  (plate_id,velocity));

            old_time = start_time;
          }

        //save the last slice
        velocity_values.push_back(std::make_pair(start_time,velocity_slice));
      }

      template <int dim>
      void
      BoxPlatesLookup<dim>::load_file(const std::string &filename,
                                      const double time)
      {
        const double time_until_end = velocity_values.back().first - time;
        unsigned int old_index,next_index;
        double velocity_time_weight;

        if (time_until_end >= velocity_values.back().first)
          {
            old_index = velocity_values.size()-1;
            next_index = velocity_values.size()-1;
            velocity_time_weight = 0.0;
          }
        else if (time_until_end <= velocity_values[0].first)
          {
            old_index = 0;
            next_index = 0;
            velocity_time_weight = 1.0;
          }
        else
          {
            for (unsigned int i = velocity_values.size() - 2; i >= 0; i--)
              if (time_until_end > velocity_values[i].first)
                {
                  old_index = i + 1;
                  next_index = i;
                  velocity_time_weight = (velocity_values[old_index].first - time_until_end) / (velocity_values[old_index].first - velocity_values[next_index].first);
                  break;
                }
          }

        const velocity_map old_map = velocity_values[old_index].second;
        const velocity_map next_map = velocity_values[next_index].second;

        Assert((0.0 <= velocity_time_weight) && (1.0 >= velocity_time_weight),
               ExcMessage ("Velocity time weight is wrong"));

        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (Utilities::fexists(filename),
                     ExcMessage (std::string("Ascii data file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open data file <"
                                             +
                                             filename
                                             +
                                             ">.")));

        // Read header lines and if necessary reinit tables
        while (in.peek() == '#')
          {
            std::string line;
            getline(in,line);
            std::stringstream linestream(line);
            std::string word;
            while (linestream >> word)
              if (word == "POINTS:")
                for (unsigned int i = 0; i < dim; i++)
                  {
                    unsigned int temp_index;
                    linestream >> temp_index;

                    if (table_points[i] == 0)
                      table_points[i] = temp_index;
                    else
                      AssertThrow (table_points[i] == temp_index,
                                   ExcMessage("The file grid must not change over model runtime. "
                                              "Either you prescribed a conflicting number of points in "
                                              "the input file, or the POINTS comment in your data files "
                                              "is changing between following files."));
                  }
          }

        for (unsigned int i = 0; i < dim; i++)
          AssertThrow (table_points[i] != 0,
                       ExcMessage("There was no POINTS: keyword in the data file, or it could not "
                                  "be parsed correctly. Ensure that at least the first data file contains "
                                  "an information about the grid size in the form '# POINTS: x y', where "
                                  "x y (only x in 2D) are the number of points in the coordinate directions"
                                  " of the boundary grid."));

        /**
         * Table for the new data. This peculiar reinit is necessary, because
         * there is no constructor for Table, which takes TableIndices as
         * argument.
         */
        Table<dim,double> data_table;
        data_table.TableBase<dim,double>::reinit(table_points);
        std::vector<Table<dim,double> > data_tables(components+dim,data_table);

        // Read data lines
        unsigned int line = 0;
        double temp_data;

        while (!in.eof())
          {
            Point<dim> position;
            for (unsigned int i = 0; i < dim; i++)
              {
                if (!(in >> position[i]))
                  break;
                data_tables[i](compute_table_indices(line)) = position[i];
              }

            char plate_id;
            if (!(in >> plate_id))
              break;

            Tensor<1,dim+1> old_velocity = old_map.find(plate_id)->second.first;
            const double old_omega = old_map.find(plate_id)->second.second;

            Tensor<1,dim+1> velocity;
            double omega;

            // It might happen that the current plate disappears in the next time step
            // In case the plate is not longer there, use the old velocity
            if (next_map.find(plate_id) != next_map.end())
              {
                velocity = next_map.find(plate_id)->second.first;
                omega    = next_map.find(plate_id)->second.second;
              }
            else
              {
                velocity = old_velocity;
                omega    = old_omega;
              }


            Tensor<1,dim+1> rotation_velocity;
            Tensor<1,dim+1> old_rotation_velocity;

            if (dim == 2)
              {
                rotation_velocity[0] = -omega*position[1];
                rotation_velocity[1] = omega*position[0];
                old_rotation_velocity[0] = -old_omega*position[1];
                old_rotation_velocity[1] = old_omega*position[0];
              }
            else if (dim == 1)
              {
                rotation_velocity[0] = -omega*position[0];
                old_rotation_velocity[0] = -old_omega*position[0];
              }

            velocity += rotation_velocity;
            old_velocity += old_rotation_velocity;

            Tensor<1,dim+1> surface_velocity;
            if (interpolate_velocity)
              {
                surface_velocity = velocity_time_weight * velocity
                                   + (1-velocity_time_weight) * old_velocity;
              }
            else
              surface_velocity = old_velocity;

            for (unsigned int i = 0; i < dim+1; i++)
              data_tables[dim+i](compute_table_indices(line)) = surface_velocity[i];

            line++;

            // TODO: add checks for coordinate ordering in data files
          }

        AssertThrow(line == data_table.n_elements(),
                    ExcMessage (std::string("Number of read in points does not match number of expected points. File corrupted?")));

        std_cxx11::array<unsigned int,dim> table_intervals;

        for (unsigned int i = 0; i < dim; i++)
          {
            table_intervals[i] = table_points[i] - 1;

            TableIndices<dim> idx;
            grid_extent[i].first = data_tables[i](idx);
            idx[i] = table_points[i] - 1;
            grid_extent[i].second = data_tables[i](idx);
          }

        for (unsigned int i = 0; i < components; i++)
          {
            if (data[i])
              delete data[i];
            data[i] = new Functions::InterpolatedUniformGridData<dim> (grid_extent,
                                                                       table_intervals,
                                                                       data_tables[dim+i]);
          }
      }


      template <int dim>
      TableIndices<dim>
      BoxPlatesLookup<dim>::compute_table_indices(const unsigned int line) const
      {
        TableIndices<dim> idx;
        idx[0] = line % table_points[0];
        if (dim >= 2)
          idx[1] = (line / table_points[0]) % table_points[1];
        if (dim == 3)
          idx[2] = line / (table_points[0] * table_points[1]);

        return idx;
      }

      template <int dim>
      double
      BoxPlatesLookup<dim>::get_data(const Point<dim> &position,
                                     const unsigned int component) const
      {
        return data[component]->value(position);
      }

    }


    template <int dim>
    PlumeBox<dim>::PlumeBox ()
      :
      current_file_number(0),
      first_data_file_model_time(0.0),
      first_data_file_number(0),
      decreasing_file_order(false),
      boundary_id(numbers::invalid_boundary_id),
      data_file_time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      scale_factor(1.0),
      lookup()
    {}


    template <int dim>
    void
    PlumeBox<dim>::initialize ()
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model())
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      // Initialize the bottom plume influx
      plume_lookup.reset(new BoundaryTemperature::internal::PlumeLookup<dim>(plume_data_directory+plume_file_name,
                                                                             this->get_pcout()));

      const std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >
      bvs = this->get_prescribed_velocity_boundary_conditions();
      for (typename std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >::const_iterator
           p = bvs.begin();
           p != bvs.end(); ++p)
        {
          if (p->second.get() == this)
            boundary_id = p->first;
        }
      AssertThrow(boundary_id != numbers::invalid_boundary_id,
                  ExcMessage("Did not find the boundary indicator for the prescribed data plugin."));

      // Set the first file number and load the first files
      current_file_number = first_data_file_number;

      const int next_file_number =
        (decreasing_file_order) ?
        current_file_number - 1
        :
        current_file_number + 1;

      if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) != "bottom")
        {
          surface_lookup.reset(new internal::BoxPlatesLookup<dim-1>(surface_data_directory+surface_velocity_file_name,
                                                                    dim,
                                                                    surface_time_scale_factor,
                                                                    surface_scale_factor,
                                                                    interpolate_velocity));

          old_surface_lookup.reset(new internal::BoxPlatesLookup<dim-1>(surface_data_directory+surface_velocity_file_name,
                                                                        dim,
                                                                        surface_time_scale_factor,
                                                                        surface_scale_factor,
                                                                        interpolate_velocity));

          const std::string surface_filename (create_surface_filename (current_file_number));
          this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                            << surface_filename << "." << std::endl << std::endl;

          if (Utilities::fexists(surface_filename))
            surface_lookup->load_file(surface_filename,0.0);
          else
            AssertThrow(false,
                        ExcMessage (std::string("Ascii data file <")
                                    +
                                    surface_filename
                                    +
                                    "> not found!"));
        }
      if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) != "top")
        {
          // Load side and bottom boundaries, but only if this is not the top boundary
          lookup.reset(new Utilities::AsciiDataLookup<dim-1> (dim,
                                                              scale_factor));

          old_lookup.reset(new Utilities::AsciiDataLookup<dim-1> (dim,
                                                                  scale_factor));

          const std::string filename (create_filename (current_file_number));
          this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                            << filename << "." << std::endl << std::endl;

          if (Utilities::fexists(filename))
            lookup->load_file(filename);
          else
            AssertThrow(false,
                        ExcMessage (std::string("Ascii data file <")
                                    +
                                    filename
                                    +
                                    "> not found!"));
        }

      // If the boundary condition is constant, switch
      // off time_dependence immediately. This also sets time_weight to 1.0.
      // If not, also load the second file for interpolation.
      if (create_filename (current_file_number) == create_filename (current_file_number+1))
        end_time_dependence ();
      else
        {
          const std::string filename (create_filename (next_file_number));
          const std::string surface_filename (create_surface_filename (next_file_number));

          if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "top")
            {
              this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                                << surface_filename << "." << std::endl << std::endl;
              if (Utilities::fexists(surface_filename))
                {
                  surface_lookup.swap(old_surface_lookup);
                  surface_lookup->load_file(surface_filename,
                                            std::abs(next_file_number-first_data_file_number)*data_file_time_step);
                }
              else
                end_time_dependence();
            }
          else if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "bottom")
            {
              this->get_pcout() << std::endl << "   Loading Ascii data data boundary file "
                                << filename << "." << std::endl << std::endl;
              if (Utilities::fexists(filename))
                {
                  lookup.swap(old_lookup);
                  lookup->load_file(filename);
                }
              else
                end_time_dependence();
            }
          else
            {
              this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                                << surface_filename << "." << std::endl << std::endl;
              this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                                << filename << "." << std::endl << std::endl;

              if (Utilities::fexists(filename) && Utilities::fexists(surface_filename))
                {
                  surface_lookup.swap(old_surface_lookup);
                  surface_lookup->load_file(surface_filename,
                                            std::abs(next_file_number-first_data_file_number)*data_file_time_step);
                  lookup.swap(old_lookup);
                  lookup->load_file(filename);
                }
              else
                end_time_dependence ();
            }
        }
    }

    template <int dim>
    std::string
    PlumeBox<dim>::create_filename (const int filenumber) const
    {
      std::string templ = data_directory+data_file_name;
      const int size = templ.length();
      char *filename = (char *) (malloc ((size + 10) * sizeof(char)));
      const std::string boundary_id_name = this->get_geometry_model().translate_id_to_symbol_name(boundary_id);
      snprintf (filename, size + 10, templ.c_str (), boundary_id_name.c_str(),filenumber);
      std::string str_filename (filename);
      free (filename);
      return str_filename;
    }

    template <int dim>
    std::string
    PlumeBox<dim>::create_surface_filename (const int filenumber) const
    {
      std::string templ = surface_data_directory+surface_id_file_names;
      const int size = templ.length();
      char *filename = (char *) (malloc ((size + 10) * sizeof(char)));
      snprintf (filename, size + 10, templ.c_str (),filenumber);
      std::string str_filename (filename);
      free (filename);
      return str_filename;
    }


    template <int dim>
    void
    PlumeBox<dim>::update ()
    {
      Interface<dim>::update ();

      // First update the plume position
      plume_position = plume_lookup->plume_position(this->get_time() - model_time_to_start_plume_tail);

      if (time_dependent && (this->get_time() - first_data_file_model_time >= 0.0))
        {
          // whether we need to update our data files. This looks so complicated
          // because we need to catch increasing and decreasing file orders and all
          // possible first_data_file_model_times and first_data_file_numbers.
          const bool need_update =
            static_cast<int> ((this->get_time() - first_data_file_model_time) / data_file_time_step)
            > std::abs(current_file_number - first_data_file_number);

          if (need_update)
            {
              // The last file, which was tried to be loaded was
              // number current_file_number +/- 1, because current_file_number
              // is the file older than the current model time
              const int old_file_number =
                (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

              //Calculate new file_number
              current_file_number =
                (decreasing_file_order) ?
                first_data_file_number
                - static_cast<unsigned int> ((this->get_time() - first_data_file_model_time) / data_file_time_step)
                :
                first_data_file_number
                + static_cast<unsigned int> ((this->get_time() - first_data_file_model_time) / data_file_time_step);

              const bool load_both_files = std::abs(current_file_number - old_file_number) >= 1;

              update_data(load_both_files);
            }

          time_weight = (this->get_time() - first_data_file_model_time) / data_file_time_step
                        - std::abs(current_file_number - first_data_file_number);

          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage (
                    "Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }


    template <int dim>
    void
    PlumeBox<dim>::update_data (const bool load_both_files)
    {

      const int next_file_number =
        (decreasing_file_order) ?
        current_file_number - 1
        :
        current_file_number + 1;

      // If the time step was large enough to move forward more
      // then one data file we need to load both current files
      // to stay accurate in interpolation
      if (load_both_files)
        {
          const std::string filename (create_filename (current_file_number));
          const std::string surface_filename (create_surface_filename (current_file_number));

          if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "top")
            {
              this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                                << surface_filename << "." << std::endl << std::endl;
              if (Utilities::fexists(surface_filename))
                {
                  surface_lookup.swap(old_surface_lookup);
                  surface_lookup->load_file(surface_filename,
                                            std::abs(next_file_number-first_data_file_number)*data_file_time_step);
                }
              else
                end_time_dependence();
            }
          else if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "bottom")
            {
              this->get_pcout() << std::endl << "   Loading Ascii data data boundary file "
                                << filename << "." << std::endl << std::endl;
              if (Utilities::fexists(filename))
                {
                  lookup.swap(old_lookup);
                  lookup->load_file(filename);
                }
              else
                end_time_dependence();
            }
          else
            {
              this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                                << surface_filename << "." << std::endl << std::endl;
              this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                                << filename << "." << std::endl << std::endl;

              if (Utilities::fexists(filename) && Utilities::fexists(surface_filename))
                {
                  surface_lookup.swap(old_surface_lookup);
                  surface_lookup->load_file(surface_filename,
                                            std::abs(next_file_number-first_data_file_number)*data_file_time_step);
                  lookup.swap(old_lookup);
                  lookup->load_file(filename);
                }
              else
                end_time_dependence ();
            }
        }

      // Now load the data file. This part is the main purpose of this function.
      const std::string filename (create_filename (next_file_number));
      const std::string surface_filename (create_surface_filename (next_file_number));

      if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "top")
        {
          this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                            << surface_filename << "." << std::endl << std::endl;
          if (Utilities::fexists(surface_filename))
            {
              surface_lookup.swap(old_surface_lookup);
              surface_lookup->load_file(surface_filename,
                                        std::abs(next_file_number-first_data_file_number)*data_file_time_step);
            }
          else
            end_time_dependence();
        }
      else if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "bottom")
        {
          this->get_pcout() << std::endl << "   Loading Ascii data data boundary file "
                            << filename << "." << std::endl << std::endl;
          if (Utilities::fexists(filename))
            {
              lookup.swap(old_lookup);
              lookup->load_file(filename);
            }
          else
            end_time_dependence();
        }
      else
        {
          this->get_pcout() << std::endl << "   Loading BoxPlates data boundary file "
                            << surface_filename << "." << std::endl << std::endl;
          this->get_pcout() << std::endl << "   Loading Ascii data boundary file "
                            << filename << "." << std::endl << std::endl;

          if (Utilities::fexists(filename) && Utilities::fexists(surface_filename))
            {
              surface_lookup.swap(old_surface_lookup);
              surface_lookup->load_file(surface_filename,
                                        std::abs(next_file_number-first_data_file_number)*data_file_time_step);
              lookup.swap(old_lookup);
              lookup->load_file(filename);
            }
          else
            end_time_dependence ();
        }
    }


    template <int dim>
    void
    PlumeBox<dim>::end_time_dependence ()
    {
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      // TODO: implement support for varying length of surface and side files
      time_dependent = false;
      // Give warning if first processor
      this->get_pcout() << std::endl
                        << "   Loading new data file did not succeed." << std::endl
                        << "   Assuming constant boundary conditions for rest of model run."
                        << std::endl << std::endl;
    }

    template <int dim>
    Tensor<1,dim>
    PlumeBox<dim>::
    get_velocity (const Point<dim> position) const
    {
      const std_cxx11::array<unsigned int,dim-1> boundary_dimensions =
        get_boundary_dimensions();

      Point<dim-1> data_position;
      for (unsigned int i = 0; i < dim-1; i++)
        data_position[i] = position[boundary_dimensions[i]];

      Tensor<1,dim> velocity,old_velocity;
      for (unsigned int i = 0; i < dim; i++)
        velocity[i] = lookup->get_data(data_position,i);

      if (!time_dependent)
        return velocity;

      for (unsigned int i = 0; i < dim; i++)
        old_velocity[i] = old_lookup->get_data(data_position,i);

      return time_weight * velocity + (1 - time_weight) * old_velocity;
    }

    template <int dim>
    Tensor<1,dim>
    PlumeBox<dim>::
    get_surface_velocity (const Point<dim> position) const
    {
      Point<dim-1> data_position;
      for (unsigned int i = 0; i < dim-1; i++)
        data_position[i] = position[i];

      Tensor<1,dim> surface_velocity,old_surface_velocity;
      for (unsigned int i = 0; i < dim; i++)
        surface_velocity[i] = surface_lookup->get_data(data_position,i);

      if (!time_dependent)
        return surface_velocity;

      for (unsigned int i = 0; i < dim; i++)
        old_surface_velocity[i] = old_surface_lookup->get_data(data_position,i);

      return time_weight * surface_velocity + (1 - time_weight) * old_surface_velocity;
    }

    template <int dim>
    std_cxx11::array<unsigned int,dim-1>
    PlumeBox<dim>::get_boundary_dimensions () const
    {
      std_cxx11::array<unsigned int,dim-1> boundary_dimensions;

      switch (dim)
        {
          case 2:
            if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
              }
            else
              AssertThrow(false,ExcNotImplemented());

            break;

          case 3:
            if ((boundary_id == 4) || (boundary_id == 5))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 1;
              }
            else if ((boundary_id == 0) || (boundary_id == 1))
              {
                boundary_dimensions[0] = 1;
                boundary_dimensions[1] = 2;
              }
            else if ((boundary_id == 2) || (boundary_id == 3))
              {
                boundary_dimensions[0] = 0;
                boundary_dimensions[1] = 2;
              }
            else
              AssertThrow(false,ExcNotImplemented());

            break;

          default:
            AssertThrow(false,ExcNotImplemented());
        }
      return boundary_dimensions;
    }

    template <int dim>
    Tensor<1,dim>
    PlumeBox<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      if (this->get_time() - first_data_file_model_time >= 0.0)
        {
          Tensor<1,dim> velocity;
          // If at the top, only use plate velocity
          if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "top")
            velocity = get_surface_velocity(position);
          // If at the bottom, add the plume contribution
          else if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "bottom")
            {
              velocity = get_velocity(position);

              double distance_head_to_boundary,current_head_radius(0);
              if ((head_radius > 0) && (head_velocity > 0))
                {
                  distance_head_to_boundary = fabs(head_velocity * (this->get_time() - model_time_to_start_plume_tail));

                  // If the plume is not yet there, perturbation will not be set
                  if (distance_head_to_boundary < head_radius)
                    {
                      current_head_radius = sqrt(head_radius * head_radius
                                                 - distance_head_to_boundary * distance_head_to_boundary);
                    }
                }

              //Normal plume tail if most of the plume head has passed
              if ((this->get_time() >= model_time_to_start_plume_tail)
                  && (current_head_radius < tail_radius))
                {
                  velocity[dim-1] += tail_velocity * std::exp(-std::pow((position-plume_position).norm()/tail_radius,2));
                }
              else
                {
                  if ((position-plume_position).norm() < current_head_radius)
                    {
                      velocity[dim-1] += head_velocity;
                      // const double head_amplitude = maximum_head_amplitude * std::exp(-std::pow(distance_head_to_boundary/maximum_head_radius,2));
                      // velocity[dim-1] += head_amplitude * std::exp(-std::pow(distance.norm()/head_radius,2));
                    }
                }
            }
          // At the sides interpolate between side and top velocity
          else
            {
              const double depth = this->get_geometry_model().depth(position);
              const double depth_weight = 0.5*(1.0 - std::tanh((depth - lithosphere_thickness)
                                                               / depth_interpolation_width));

              velocity = (1 - depth_weight) * get_velocity(position);
              velocity += depth_weight * get_surface_velocity(position);
            }

          return velocity;
        }
      else
        return Tensor<1,dim> ();
    }


    template <int dim>
    void
    PlumeBox<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary velocity model");
      {
        prm.enter_subsection ("Ascii data model");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/ascii-data/test/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Data file name",
                             "box_2d_%s.%d.csv",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Velocity file name).\\%s%d where \\\\%s is a string specifying "
                             "the boundary of the model according to the names of the boundary "
                             "indicators (of a box or a spherical shell).%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number. ");
          prm.declare_entry ("Data file time step", "1e6",
                             Patterns::Double (0),
                             "Time step between following velocity files. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years. "
                             "The default is one million, i.e., either one million seconds or one million years.");
          prm.declare_entry ("First data file model time", "0",
                             Patterns::Double (0),
                             "Time from which on the velocity file with number 'First velocity "
                             "file number' is used as boundary condition. Previous to this "
                             "time, a no-slip boundary condition is assumed. Depending on the setting "
                             "of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("First data file number", "0",
                             Patterns::Integer (),
                             "Number of the first velocity file to be loaded when the model time "
                             "is larger than 'First velocity file model time'.");
          prm.declare_entry ("Decreasing file order", "false",
                             Patterns::Bool (),
                             "In some cases the boundary files are not numbered in increasing "
                             "but in decreasing order (e.g. 'Ma BP'). If this flag is set to "
                             "'True' the plugin will first load the file with the number "
                             "'First velocity file number' and decrease the file number during "
                             "the model run.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction. Another way to use this factor is to "
                             "convert units of the input files. The unit is assumed to be"
                             "m/s or m/yr depending on the 'Use years in output instead of "
                             "seconds' flag. If you provide velocities in cm/yr set this "
                             "factor to 0.01.");
        }
        prm.leave_subsection();

        prm.enter_subsection ("Box plates model");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/tristan_plume/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Id file names", "plates_%d",
                             Patterns::Anything (),
                             "The file name of the id data. Provide file in format: "
                             "(Velocity file name).\\%d.gpml where \\%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number.");
          prm.declare_entry ("Velocity file name", "velocities",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Velocity file name).\\%d.gpml where \\%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number.");
          prm.declare_entry ("Lithosphere thickness", "2e5",
                             Patterns::Double (0),
                             "The velocity at the sides is interpolated between a plate velocity "
                             "at the surface and ascii data velocities below the lithosphere. "
                             "This value is the depth of the transition between surface and side "
                             "velocities.");
          prm.declare_entry ("Depth interpolation width", "5e4",
                             Patterns::Double (0),
                             "The width of the interpolation zone described for 'lithosphere "
                             "thickness'.");
          prm.declare_entry ("Time scale factor", "1e6",
                             Patterns::Double (0),
                             "Determines the factor applied to the times in the velocity file. The unit is assumed to be"
                             "s or yr depending on the 'Use years in output instead of "
                             "seconds' flag. If you provide time in Ma set this "
                             "factor to 1e6.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction.");
          prm.declare_entry ("Interpolate plate velocity", "true",
                             Patterns::Bool (),
                             "A boolean value, which determines whether to interpolate "
                             "between old and new plate velocities. If false, only the "
                             "old one will be used for a piecewise-constant velocity.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Plume");
      {
        prm.declare_entry ("Data directory",
                           "$ASPECT_SOURCE_DIR/data/boundary-temperature/plume/",
                           Patterns::DirectoryName (),
                           "The name of a directory that contains the model data. This path "
                           "may either be absolute (if starting with a '/') or relative to "
                           "the current directory. The path may also include the special "
                           "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the 'data/' subdirectory of ASPECT. ");
        prm.declare_entry ("Plume position file name", "Tristan.sur",
                           Patterns::Anything (),
                           "The file name of the plume position data.");
        prm.declare_entry ("Inflow velocity", "0",
                           Patterns::Double (),
                           "Magnitude of the velocity inflow. Units: m/s or m/yr.");
        prm.declare_entry ("Radius", "0",
                           Patterns::Double (),
                           "Radius of the anomaly. Units: m.");
        prm.declare_entry ("Head radius", "0",
                           Patterns::Double (),
                           "Radius of the plume head velocity anomaly. Units: m.");
        prm.declare_entry ("Head velocity", "0",
                           Patterns::Double (),
                           "Magnitude of the plume head velocity inflow. Units: m/s or m/yr.");
        prm.declare_entry ("Model time to start plume tail", "0",
                           Patterns::Double (),
                           "Time before the start of the plume position data at which "
                           "the head starts to flow into the model. Units: years or "
                           "seconds.");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    PlumeBox<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Ascii data model");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory    = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = data_directory.find (subst_text),  position!=std::string::npos)
              data_directory.replace (data_directory.begin()+position,
                                      data_directory.begin()+position+subst_text.size(),
                                      ASPECT_SOURCE_DIR);
          }

          data_file_name    = prm.get ("Data file name");

          scale_factor      = prm.get_double ("Scale factor");

          data_file_time_step             = prm.get_double ("Data file time step");
          first_data_file_model_time      = prm.get_double ("First data file model time");
          first_data_file_number          = prm.get_double ("First data file number");
          decreasing_file_order           = prm.get_bool   ("Decreasing file order");

          if (this->convert_output_to_years() == true)
            {
              data_file_time_step        *= year_in_seconds;
              first_data_file_model_time *= year_in_seconds;
              scale_factor               /= year_in_seconds;
            }
        }
        prm.leave_subsection();

        prm.enter_subsection("Box plates model");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          surface_data_directory        = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = surface_data_directory.find (subst_text),  position!=std::string::npos)
              surface_data_directory.replace (surface_data_directory.begin()+position,
                                              surface_data_directory.begin()+position+subst_text.size(),
                                              ASPECT_SOURCE_DIR);
          }

          surface_velocity_file_name    = prm.get ("Velocity file name");
          surface_id_file_names         = prm.get ("Id file names");

          lithosphere_thickness = prm.get_double ("Lithosphere thickness");
          depth_interpolation_width     = prm.get_double ("Depth interpolation width");
          surface_time_scale_factor     = prm.get_double ("Time scale factor");
          surface_scale_factor  = prm.get_double ("Scale factor");
          interpolate_velocity  = prm.get_bool ("Interpolate plate velocity");


          if (this->convert_output_to_years() == true)
            {
              surface_time_scale_factor *= year_in_seconds;
              surface_scale_factor      /= year_in_seconds;
            }
        }
        prm.leave_subsection();

      }
      prm.leave_subsection();

      prm.enter_subsection("Plume");
      {
        // Get the path to the data files. If it contains a reference
        // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
        // as a #define
        plume_data_directory        = prm.get ("Data directory");
        {
          const std::string      subst_text = "$ASPECT_SOURCE_DIR";
          std::string::size_type position;
          while (position = plume_data_directory.find (subst_text),  position!=std::string::npos)
            plume_data_directory.replace (plume_data_directory.begin()+position,
                                          plume_data_directory.begin()+position+subst_text.size(),
                                          ASPECT_SOURCE_DIR);
        }

        plume_file_name    = prm.get ("Plume position file name");
        tail_velocity = prm.get_double ("Inflow velocity");
        tail_radius = prm.get_double ("Radius");

        head_radius = prm.get_double("Head radius");
        head_velocity = prm.get_double("Head velocity");
        model_time_to_start_plume_tail = prm.get_double ("Model time to start plume tail");

        if (this->convert_output_to_years() == true)
          {
            tail_velocity /= year_in_seconds;
            head_velocity /= year_in_seconds;
            model_time_to_start_plume_tail *= year_in_seconds;
          }
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(PlumeBox,
                                                 "plume box",
                                                 "This is a velocity plugin that is a combination "
                                                 "of the Ascii data, Box plates and Plume plugins. "
                                                 "It is specifically designed for modelling a box "
                                                 "of the mantle with prescribed velocities at every "
                                                 "boundary with a plume influx at the bottom, a plate "
                                                 "velocity at the top, side and bottom velocities from "
                                                 "a different geodynamic model and interpolation "
                                                 "between these constraints.")
  }
}

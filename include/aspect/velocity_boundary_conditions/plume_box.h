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


#ifndef __aspect__velocity_boundary_conditions_plume_box_h
#define __aspect__velocity_boundary_conditions_plume_box_h

#include <aspect/velocity_boundary_conditions/interface.h>

// Additional lookup classes are within these
#include <aspect/boundary_temperature/plume.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    namespace internal
    {
      /**
       * AsciiDataLookup reads in files containing input data
         in ascii format. Note the required format of the
         input data: The first lines may contain any number of comments
         if they begin with '#', but one of these lines needs to
         contain the number of grid points in each dimension as
         for example '# POINTS: 3 3'. The order of the columns
         has to be 'coordinates data' with @p dim coordinate columns
         and @p components data columns. Note that the data in the input
         files need to be sorted in a specific order:
         the first coordinate needs to ascend first,
         followed by the second and so on in order to
         assign the correct data to the prescribed coordinates.
       */
      template <int dim>
      class BoxPlatesLookup
      {
        public:
          BoxPlatesLookup(const std::string &filename,
                          const unsigned int components,
                          const double time_scale_factor,
                          const double velocity_scale_factor,
                          const bool interpolate_velocity);

          /**
           * Loads a data text file. Throws an exception if the file does not exist,
           * if the data file format is incorrect or if the file grid changes over model runtime.
           */
          void
          load_file(const std::string &filename,
                    const double time);

          /**
           * Returns the computed velocity
           * in cartesian coordinates.
           *
           * @param position The current position to compute the data (velocity, temperature, etc.)
           * @param component Number of the component that is requested
           */
          double
          get_data(const Point<dim> &position,
                   const unsigned int component) const;

        private:
          /**
           * The number of data components read in (=columns in the data file).
           */
          const unsigned int components;

          /**
           * Interpolation functions to access the data.
           */
          std::vector<Functions::InterpolatedUniformGridData<dim> *> data;

          /**
           * Model size
           */
          std_cxx11::array<std::pair<double,double>,dim> grid_extent;

          /**
           * Number of points in the data grid.
           */
          TableIndices<dim> table_points;

          /**
           * Scales the data times by a scalar factor. Can be
           * used to transform the unit of the data.
           */
          const double time_scale_factor;

          /**
           * Scales the data boundary velocity by a scalar factor. Can be
           * used to transform the unit of the data.
           */
          const double velocity_scale_factor;

          /**
           * Determines whether to interpolate between old and new plate
           * velocities. If false, only the old one will be used.
           */
          const bool interpolate_velocity;

          /**
           * Any plate velocity consists of a velocity and the rotation
           * of the associated plate in the used map projection
           */
          typedef std::pair<Tensor<1,dim+1>,double> plate_velocity;

          typedef std::map<unsigned char,plate_velocity > velocity_map;

          /**
           * Table for the plate velocities of all times.
           */
          std::vector<std::pair<double, velocity_map> > velocity_values;

          /**
           * Computes the table indices of each entry in the input data file.
           * The index depends on dim, grid_dim and the number of components.
           */
          TableIndices<dim>
          compute_table_indices(const unsigned int line) const;
      };
    }



    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a AsciiData input file.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class PlumeBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        PlumeBox ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next data files if
         * necessary and outputs a warning if the end of the set of data
         * files is reached.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from the text files.
         */
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * A variable that stores the currently used data file of a
         * series. It gets updated if necessary by set_current_time.
         */
        int  current_file_number;

        /**
         * Time from which on the data file with number 'First data
         * file number' is used as boundary condition. Previous to this
         * time, a no-slip boundary condition is assumed. Depending on the setting
         * of the global 'Use years in output instead of seconds' flag
         * in the input file, this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

        /**
         * Number of the first data file to be loaded when the model time
         * is larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing
         * but in decreasing order (e.g. 'Ma BP'). If this flag is set to
         * 'True' the plugin will first load the file with the number
         * 'First data file number' and decrease the file number during
         * the model run.
         */
        bool decreasing_file_order;

        /**
         * Directory in which the data files are present.
         */
        std::string data_directory;

        /**
         * Filename of data file. The file names can contain
         * the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Directory in which the surface velocity files are present.
         */
        std::string surface_data_directory;

        /**
         * Filename of velocity file.
         */
        std::string surface_velocity_file_name;

        /**
         * First part of filename of id files. The files have to have
         * the pattern velocity_file_name.n.gpml where n is the number of the
         * current timestep (starts from 0).
         */
        std::string surface_id_file_names;

        /**
         * Distance between two surface grid points in x direction.
         */
        double x_step;

        /**
         * Distance between two surface grid points in y direction.
         */
        double y_step;

        /**
         * Boundary id for which this plugin is created. Is initialized in
         * initialize() and used to determine the appropriate file name.
         */
        types::boundary_id boundary_id;

        /**
         * Time in model units (depends on other model inputs) between two
         * data files.
         */
        double data_file_time_step;

        /**
         * Number of grid points in data file
         */
        std_cxx1x::array<unsigned int,3> data_points;

        /**
         * Weight between data file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more data files to suppress attempts to
         * read in new files.
         */
        bool time_dependent;

        /**
         * Scale the boundary condition by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double scale_factor;

        /**
         * Determines the width of the velocity interpolation zone around the
         * current point. Currently equals the distance between evaluation
         * point and velocity data point that is still included in the
         * interpolation. The weighting of the points currently only accounts
         * for the surface area a single data point is covering ('moving
         * window' interpolation without distance weighting).
         */
        double interpolation_width;

        /**
         * The velocity at the sides is interpolated between a plate velocity
         * at the surface and ascii data velocities below the lithosphere.
         * This value is the depth of the transition between surface and side
         * velocities.
         */
        double lithosphere_thickness;

        /**
         * The width of the interpolation zone described for 'lithosphere
         * thickness'.
         */
        double depth_interpolation_width;

        /**
         * Scale the time steps of the surface by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double surface_time_scale_factor;

        /**
         * Scale the surface velocity by a scalar factor. Can be
         * used to transform the unit of the velocities (if they are not
         * specified in the default unit (m/s or m/yr depending on the
         * "Use years in output instead of seconds" parameter).
         */
        double surface_scale_factor;

        /**
         * Determines whether to interpolate between old and new plate
         * velocities. If false, only the old one will be used.
         */
        bool interpolate_velocity;

        /**
         * Pointer to an object that reads and processes data we get from
         * text files.
         */
        std_cxx11::shared_ptr<Utilities::AsciiDataLookup<dim-1> > lookup;
        std_cxx11::shared_ptr<Utilities::AsciiDataLookup<dim-1> > old_lookup;


        /**
         * Pointer to an object that reads and processes data we get from
         * text files.
         */
        std_cxx11::shared_ptr<internal::BoxPlatesLookup<dim-1> > surface_lookup;
        std_cxx11::shared_ptr<internal::BoxPlatesLookup<dim-1> > old_surface_lookup;


        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files.
         */
        std_cxx1x::shared_ptr<BoundaryTemperature::internal::PlumeLookup<dim> > plume_lookup;

        /**
         * Current plume position. Used to avoid looking up the plume position
         * for every quadrature point, since this is only time-dependent.
         */
        Point<dim> plume_position;

        /**
         * Filename of plume file.
         */
        std::string plume_file_name;

        /**
         * Directory in which the plume file is present.
         */
        std::string plume_data_directory;

        /**
         * Magnitude of the plume velocity anomaly
         */
        double tail_velocity;

        /**
         * Radius of the plume velocity anomaly
         */
        double tail_radius;

        /**
         * Radius of the plume head temperature anomaly
         */
        double head_radius;

        /**
         * Velocity of the plume head inflow
         */
        double head_velocity;

        /**
         * Model time at which the plume tail will start to move according to
         * the position data file. This is equivalent to the difference between
         * model start time (in the past) and the first data time of the plume
         * positions file.
         */
        double model_time_to_start_plume_tail;

        /**
         * Handles the update of the data in lookup.
         */
        void
        update_data (const bool load_both_files);

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence ();

        /**
         * Create a filename out of the name template for the side and bottom
         * ascii data model.
         */
        std::string
        create_filename (const int filenumber) const;

        /**
         * Create a filename out of the name template for the surface
         * box plates model.
         */
        std::string
        create_surface_filename (const int filenumber) const;

        Tensor<1,dim>
        get_velocity (const Point<dim> position) const;

        Tensor<1,dim>
        get_surface_velocity (const Point<dim> position) const;

        /**
         * Determines which of the dimensions of the position is used to find
         * the data point in the data grid. E.g. the left boundary of a box model extents in
         * the y and z direction (position[1] and position[2]), therefore the function
         * would return [1,2] for dim==3 or [1] for dim==2. We are lucky that these indices are
         * identical for the box and the spherical shell (if we use spherical coordinates for the
         * spherical shell), therefore we do not need to distinguish between them. For the initial
         * condition this function is trivial, because the position in the data grid is the same as
         * the actual position (the function returns [0,1,2] or [0,1]), but for the boundary
         * conditions it matters.
         */
        std_cxx11::array<unsigned int,dim-1>
        get_boundary_dimensions () const;
    };
  }
}


#endif

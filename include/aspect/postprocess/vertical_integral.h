/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_vertical_integral_h
#define __aspect__postprocess_vertical_integral_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A class derived from CellDataVectorCreator that takes an output
     * vector and computes a variable that represents the vertical integral
     * of a given compositional field. This quantity only makes sense at the
     * surface of the domain. Thus, the value is set to zero in all the
     * cells inside of the domain.
     *
     * The member functions are all implementations of those declared in the
     * base class. See there for their meaning.
     */
    template <int dim>
    class VerticalIntegral
      : public Interface<dim>,
        public SimulatorAccess<dim>
    {
      public:
        VerticalIntegral();

        /**
         * Evaluate the solution for the vertical integral of a given
         * compositional field.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


        /**
         * Save the state of this object.
         */
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object.
         */
        virtual
        void load (const std::map<std::string, std::string> &status_strings);

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * A parameter that we read from the input file that denotes, which
         * compositional field to integrate.
         */
        std::string name_of_compositional_field;

        /**
         * The format in which to produce graphical output. This also
         * determines the extension of the file name to which to write.
         */
        DataOutBase::OutputFormat output_format;

        /**
         * A parameter that can be used to exclude the upper part
         * of the model from integration. All cells with a smaller
         * depth are ignored.
         */
        double minimum_depth;

        /**
         * A parameter that can be used to exclude the lower part
         * of the model from integration. All cells with a larger
         * depth are ignored. The default value of 0 will be
         * replaced by the maximum depth of the model.
         */
        double maximum_depth;

        /**
         * Interval between the generation of graphical output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time;

        /**
         * Consecutively counted number indicating the how-manyth time we will
         * create output the next time we get to it.
         */
        unsigned int output_file_number;

        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * A list of pairs (time, pvtu_filename) that have so far been written
         * and that we will pass to DataOutInterface::write_pvd_record to
         * create a master file that can make the association between
         * simulation time and corresponding file name (this is done because
         * there is no way to store the simulation time inside the .pvtu or
         * .vtu files).
         */
        std::vector<std::pair<double,std::string> > times_and_vtu_names;

        /**
         * A list of list of filenames, sorted by timestep, that correspond to
         * what has been created as output. This is used to create a master
         * .visit file for the entire simulation.
         */
        std::vector<std::vector<std::string> > output_file_names_by_timestep;
    };
  }
}

#endif

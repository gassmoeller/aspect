/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_entropy_model_h
#define _aspect_material_model_entropy_model_h

#include <aspect/material_model/utilities.h>

#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace MaterialUtilities
    {
      namespace Lookup
      {
        /**
         * An implementation of the above base class that reads in files created
         * by the Perplex software.
         */
        class PerplexEntropyReader : public MaterialLookup
        {
          public:
            PerplexEntropyReader(const std::string &filename,
                                 const MPI_Comm &comm);

            double
            temperature(const double temperature,
                        const double pressure) const;

          private:
            dealii::Table<2,double> temperature_values;
        };
      }
    }

    /**
     * A model that is formulated in terms of pressure and entropy
     * to allow more realistic thermodynamic computations.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class EntropyModel: public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        virtual
        void
        initialize ();

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
         * @}
         */

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


      private:
        /**
         * Constant reference viscosity.
         */
        double reference_eta;

        /**
         * The value for thermal conductivity. This model only
         * implements a constant thermal conductivity for the whole domain.
         */
        double thermal_conductivity_value;

        /**
         * Information about the location of data files.
         */
        std::string data_directory;
        std::string material_file_name;

        /**
         * List of pointers to objects that read and process data we get from
         * Perplex files.
         */
        std::unique_ptr<MaterialUtilities::Lookup::PerplexEntropyReader> material_lookup;
    };
  }
}

#endif

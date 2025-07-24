/*
  Copyright (C) 2024 - by the authors of the ASPECT code.

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

#include <aspect/material_model/thermal_conductivity/hofmeister_branlund_2015.h>

// Helper functions in anonymous namespace to compute thermal conductivities using the Hofmeister and Branlund (2015) formulations
namespace 
{
  double compute_lattice_thermal_conductivity_hbr2015(const double latTC_room, const double T_dep_factor, const double T_coefficient, const double T_model)
  {
   double lattice_thermal_conductivity = latTC_room + (T_dep_factor / (T_model + (T_coefficient*T_model) ));
   return lattice_thermal_conductivity;
  }
}


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {   
      // Main function: 
      template <int dim>
      void
      hofmeister_branlund_2015<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                            MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for oceanic crust (pyroxene, plagioclase)
        constexpr int oceanicrus_index = 0;
        const double oceanicrus_hbr15_latTC_room =  1.217;     // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double oceanicrus_hbr15_latTC_Tdep =  407.34;    // T-dependent factor (Tdep) in [W m^-1]
        const double oceanicrus_hbr15_latTC_Tcoe =  8.0555e-5; // temperature coefficient (Tcoe) in [/]
          
        unsigned int litholopar_index = oceanicrus_index+1; // Number of lithologies

        // Define room temperature [K] 
        // const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {

          // Preallocate a vector for storing thermal conductivities of lithologies
          std::vector<double> hbr15_litholog_latTcond(litholopar_index, 0.0); // Lattice thermal conductivity
          std::vector<double> hbr15_litholog_totTcond(litholopar_index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of lithologies
          std::vector<std::vector<double>> hbr15_all_litholog_Tconds(litholopar_index, std::vector<double>(3, 0.0));

          // Convert pressure unit from [Pa] to [MPa]
          // double P_MPa = in.pressure[i]/1e6;
          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];
          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];

          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given lihtologies
          {
            case oceanicrus_index: // Oceanic Crust  
            { 
             double oceanicrus_hbr15_latTcon = compute_lattice_thermal_conductivity_hbr2015(
             oceanicrus_hbr15_latTC_room,
             oceanicrus_hbr15_latTC_Tdep,
             oceanicrus_hbr15_latTC_Tcoe,
             T_mod); 
             double oceanicrus_hbr15_totTcon = oceanicrus_hbr15_latTcon; 
             // Store the thermal conductivities in the vector
             hbr15_litholog_latTcond[oceanicrus_index] = oceanicrus_hbr15_latTcon;
             hbr15_litholog_totTcond[oceanicrus_index] = oceanicrus_hbr15_totTcon;
             break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < litholopar_index; ++row)
          {
            hbr15_all_litholog_Tconds[row][0] = hbr15_litholog_latTcond[row]; // Column 0: Lattice conductivities
            hbr15_all_litholog_Tconds[row][1] = hbr15_litholog_totTcond[row]; // Column 1: Total conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_hbr15_Tcond = std::pow(hbr15_all_litholog_Tconds[mID][1], min_frac);

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_hbr15_Tcond;

        }
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      #define INSTANTIATE(dim) \
      template class hofmeister_branlund_2015<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

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

#include <aspect/material_model/thermal_conductivity/nondimensional.h>

// Helper functions in anonymous namespace to compute thermal conductivities using the Marzotto et al. (2025) formulations
namespace 
{
  // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_lattice_thermal_conductivity_nondimn(double a0, double b1, double ymin, double ymax, double P_log, double T_mod, double T_room, double n_exp)
  { 
    double zsimpl = a0 + b1 * P_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = ymin + (ymax - ymin) * yprime;
    double lattice_thermal_conductivity = std::exp(yreals);
    return lattice_thermal_conductivity * std::pow((T_room / T_mod), n_exp);
  }
  
  // Compute the radiative thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_radiative_thermal_conductivity_nondimn(double c0, double d1, double jmin, double jmax, double T_log)
  {
    double zsimpl = c0 + d1 * T_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = jmin + (jmax - jmin) * yprime;
    return std::exp(yreals);
  }
     
  double compute_total_thermal_conductivity_nondimn(double lattice_thermal_conductivity, double radiative_thermal_conductivity)
  {
    double total_thermal_conductivity = lattice_thermal_conductivity + radiative_thermal_conductivity;
    return total_thermal_conductivity;
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
      nondimensional<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Define coefficients for lattice thermal conductivity of the mantle

        // Coefficients for NON-DIMENSIONAL dry olivine 
        // retreived from fitting TDTR dataset of
        // [Chang et al., 2017, PNAS, vol 114, p. 4078-4081]
        // https://doi.org/10.1073/pnas.1616216114
        // mineral composition [Mg1.8 Fe0.2 SiO4]    
        constexpr int nondimoliv_index = 0;
        const double nondimoliv_latTC_a0 =    1.6023;
        const double nondimoliv_latTC_b1 =    2.1469;
        const double nondimoliv_latTC_ymin =  1.28093384543429;
        const double nondimoliv_latTC_ymax =  2.60726381956037;  
        const double nondimoliv_Tdep_n_exp =  0.5;  

        unsigned int mineralpar_index = nondimoliv_index+1; // Number of minerals

        // Coefficients for NON-DIMENSIONAL dry olivine
        // retreived from fitting dataset of
        // [Marzotto et al. 2025, Nature Communication, 16, 6058]
        // https://doi.org/10.1038/s41467-025-61148-8
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double nondimoliv_radTC_c0 =   -5.6162;
        const double nondimoliv_radTC_d1 =    1.8839;
        const double nondimoliv_radTC_jmin = -23.02585093;
        const double nondimoliv_radTC_jmax =  1.289885976;

        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Define room temperature [K] 
        const double T_room = 298.15; 

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> nondim_minerals_latTcond(mineralpar_index, 0.0); // Lattice thermal conductivity
          std::vector<double> nondim_minerals_radTcond(mineralpar_index, 0.0); // Radiative thermal conductivity
          std::vector<double> nondim_minerals_totTcond(mineralpar_index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> nondim_all_minerals_Tconds(mineralpar_index, std::vector<double>(3, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;

          // Compute natural logarithm of pressure and temperature 
          double P_log = std::log(P_GPa);
          double T_log = std::log(in.temperature[i]);

          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];

          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];

          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
          {
            case nondimoliv_index: // Dry Olivine
            {      
              double nondimoliv_latTCon = compute_lattice_thermal_conductivity_nondimn(
              nondimoliv_latTC_a0, 
              nondimoliv_latTC_b1, 
              nondimoliv_latTC_ymin, 
              nondimoliv_latTC_ymax,
              P_log, 
              T_mod, 
              T_room, 
              nondimoliv_Tdep_n_exp); 
              double nondimoliv_radTCon = compute_radiative_thermal_conductivity_nondimn(
              nondimoliv_radTC_c0, 
              nondimoliv_radTC_d1, 
              nondimoliv_radTC_jmin, 
              nondimoliv_radTC_jmax, 
              T_log); 
              double nondimoliv_TotTCon = compute_total_thermal_conductivity_nondimn(
              nondimoliv_latTCon, 
              nondimoliv_radTCon); 
              // Store the thermal conductivities in the vector
              nondim_minerals_latTcond[nondimoliv_index] = nondimoliv_latTCon;
              nondim_minerals_radTcond[nondimoliv_index] = nondimoliv_radTCon;
              nondim_minerals_totTcond[nondimoliv_index] = nondimoliv_TotTCon;
              break;
            }

          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            nondim_all_minerals_Tconds[row][0] = nondim_minerals_latTcond[row]; // Column 0: Lattice conductivities
            nondim_all_minerals_Tconds[row][1] = nondim_minerals_radTcond[row]; // Column 1: Radiative conductivities
            nondim_all_minerals_Tconds[row][2] = nondim_minerals_totTcond[row]; // Column 2: Total conductivities
          }


          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_nondi_Tcond = std::pow(nondim_all_minerals_Tconds[mID][2], min_frac);
          std::vector<double> nondimensional_Tcond(5, aggrock_testcase_nondi_Tcond);

          // Test Case
          out.thermal_conductivities = nondimensional_Tcond;

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
      template class nondimensional<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

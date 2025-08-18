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

#include <aspect/material_model/thermal_conductivity/nondimensional_Tcond.h>

// Helper functions in anonymous namespace to compute thermal conductivities using the Marzotto et al. (2025) formulations
namespace 
{
  // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_lattice_thermal_conductivity_nondimn(double a0, double b1, double ymin, double ymax, double P_log, double T_mod, double T_nondim, double n_exp)
  { 
    double zsimpl = a0 + b1 * P_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = ymin + (ymax - ymin) * yprime;
    double lattice_thermal_conductivity = std::exp(yreals);
    return lattice_thermal_conductivity * std::pow((T_nondim / T_mod), n_exp);
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
      nondimensional_Tcond<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
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

        // Define coefficients for radiative thermal conductivity of different minerals

        // Coefficients for NON-DIMENSIONAL dry olivine
        // retreived from fitting dataset of
        // [Marzotto et al. 2025, Nature Communication, 16, 6058]
        // https://doi.org/10.1038/s41467-025-61148-8
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double nondimoliv_radTC_c0 =    5.6162;
        const double nondimoliv_radTC_d1 =    1.8839;
        const double nondimoliv_radTC_jmin = -23.0258509299405;
        const double nondimoliv_radTC_jmax =  1.28988597569074;

        // Preallocate a vector for mineral fractions of different rocks and a vector for mineral indices

        // Uniform Lithology - Upper Mantle  (100% olivine) 
        std::vector<double> minfract_unifolit_UpMa = {1.00};
        std::vector<unsigned int> minindex_unifolit_UpMa = {nondimoliv_index};
        // Uniform Lithology - Upper Mantle Transition Zone (100% olivine)
        std::vector<double> minfract_unifolit_UMTZ = {1.00};
        std::vector<unsigned int> minindex_unifolit_UMTZ = {nondimoliv_index};
        // Uniform Lithology - Lower Mantle Transition Zone (100% olivine)
        std::vector<double> minfract_unifolit_LMTZ = {1.00};
        std::vector<unsigned int> minindex_unifolit_LMTZ = {nondimoliv_index};
        // Uniform Lithology - Lower Mantle (100% olivine)
        std::vector<double> minfract_unifolit_LoMa = {1.00};
        std::vector<unsigned int> minindex_unifolit_LoMa = {nondimoliv_index};

        // All Minerals (test)
        std::vector<double> minfract_allminerals_test(mineralpar_index);
        std::vector<unsigned int> minindex_allminerals_test(mineralpar_index);
        for (unsigned int i = 0; i < mineralpar_index; ++i)
        {
          minindex_allminerals_test[i] = i;
          minfract_allminerals_test[i] = 1.00;
        }

        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Check if the sum of Rock Mineral Fraction is equal to 1
        double sum_min_fract_unifolit_UpMa = std::accumulate(minfract_unifolit_UpMa.begin(), minfract_unifolit_UpMa.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_unifolit_UpMa - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_unifolit_UM must be equal to 1."));
        double sum_min_fract_unifolit_UMTZ = std::accumulate(minfract_unifolit_UMTZ.begin(), minfract_unifolit_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_unifolit_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_unifolit_UMTZ must be equal to 1."));
        double sum_min_fract_unifolit_LMTZ = std::accumulate(minfract_unifolit_LMTZ.begin(), minfract_unifolit_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_unifolit_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_unifolit_LMTZ must be equal to 1."));
        double sum_min_fract_unifolit_LoMa = std::accumulate(minfract_unifolit_LoMa.begin(), minfract_unifolit_LoMa.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_unifolit_LoMa - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_unifolit_LoMa must be equal to 1."));

        // Define room temperature [K] 
        const double T_room = 298.15; 

        // Define maximum temperature [K] and pressure temperature [Pa]
        const double T_max = 4000;
        const double P_max = 132.0222017e9; 

        // Define minimum temperature [K] and pressure temperature [Pa]
        const double T_min = 273.15;
        const double P_min = 1e-10; 

        // Define pressure at the top and at the bottom of each mantle layer [Pa]
        const double P_UpMa_top = P_min;         // Upper Mantle top pressure in [Pa]
        const double P_UpMa_bot = 13.59280101e9; // Upper Mantle bottom pressure in [Pa]
        const double P_UMTZ_top = 13.59280101e9; // Upper Mantle Transition Zone top pressure in [Pa]
        const double P_UMTZ_bot = 17.69264984e9; // Upper Mantle Transition Zone bottom pressure in [Pa]
        const double P_LMTZ_top = 17.69264984e9; // Lower Mantle Transition Zone top pressure in [Pa]
        const double P_LMTZ_bot = 23.11221520e9; // Lower Mantle Transition Zone bottom pressure in [Pa]
        const double P_LoMa_top = 23.11221520e9; // Lower Mantle top pressure in [Pa]
        const double P_LoMa_bot = P_max;         // Lower Mantle bottom pressure in [Pa]

        // Define pressure ratio [/]
        const double P_ratio_UpMa_top = P_UpMa_top / P_max;
        const double P_ratio_UpMa_bot = P_UpMa_bot / P_max;
        const double P_ratio_UMTZ_top = P_UMTZ_top / P_max;
        const double P_ratio_UMTZ_bot = P_UMTZ_bot / P_max;
        const double P_ratio_LMTZ_top = P_LMTZ_top / P_max;
        const double P_ratio_LMTZ_bot = P_LMTZ_bot / P_max;
        const double P_ratio_LoMa_top = P_LoMa_top / P_max;
        const double P_ratio_LoMa_bot = P_LoMa_bot / P_max;

        // Define nondimensional temperature ratio [/] 
        const double T_nondim = T_room / T_max;

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          double current_temperature = in.temperature[i];
          if (in.temperature[i] <= 0.0) // Avoid log(0) or negative temperature
          {
            current_temperature = T_min; 
          }

          double current_pressure = in.pressure[i];
          if (in.pressure[i] <= 0.0) // Avoid log(0) or negative pressure
          {
            current_pressure = P_min; 
          }

          // Compute pressure and temperature ratios
          double P_ratio = current_pressure / P_max;
          double T_ratio = current_temperature / T_max;

          // Compute natural logarithm of pressure and temperature 
          double P_log = std::log(P_ratio);
          double T_log = std::log(T_ratio);

          // Take lithology of the model
          double lithology = in.composition[0][i];

          std::vector<double> mineral_fraction;    // Mineral fractions for the current lithology
          std::vector<unsigned int> mineral_index; // Mineral indexes for the current lithology

          if (P_ratio >= P_ratio_UpMa_top && P_ratio <= P_ratio_UpMa_bot) // Upper Mantle
          {
            if (lithology == 0) // Uniform Lithology
            {
              mineral_fraction = minfract_unifolit_UpMa; 
              mineral_index = minindex_unifolit_UpMa;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_fraction = {1.00}; 
              mineral_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the upper mantle."));
            }
          }
          else if (P_ratio > P_ratio_UMTZ_top && P_ratio <= P_ratio_UMTZ_bot) // Upper Transition Zone
          {
            if (lithology == 0) // Uniform Lithology
            {
              mineral_fraction = minfract_unifolit_UMTZ; 
              mineral_index = minindex_unifolit_UMTZ;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_fraction = {1.00}; 
              mineral_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the upper mantle transition zone."));
            }
          }
          else if (P_ratio > P_ratio_LMTZ_top && P_ratio <= P_ratio_LMTZ_bot) // Lower Transition Zone
          {
            if (lithology == 0) // Uniform Lithology
            {
              mineral_fraction = minfract_unifolit_LMTZ; 
              mineral_index = minindex_unifolit_LMTZ;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_fraction = {1.00}; 
              mineral_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the lower mantle transition zone."));
            }
          }
          else if (P_ratio > P_ratio_LoMa_top && P_ratio <= P_ratio_LoMa_bot) // Lower Mantle
          {
            if (lithology == 0) // Uniform Lithology
            {
              mineral_fraction = minfract_unifolit_LoMa; 
              mineral_index = minindex_unifolit_LoMa;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_fraction = {1.00}; 
              mineral_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the lower mantle."));
            }
          }
          else
          {
           AssertThrow(false, dealii::ExcMessage("Invalid pressure range for the mantle."));
          }

          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> nondim_minerals_latTcond(mineral_fraction.size(), 0.0); // Lattice thermal conductivity
          std::vector<double> nondim_minerals_radTcond(mineral_fraction.size(), 0.0); // Radiative thermal conductivity
          std::vector<double> nondim_minerals_totTcond(mineral_fraction.size(), 0.0); // Total thermal conductivity

          // Preallocate total thermal conductivity of the aggregate rock
          double nondim_aggregate_rock_totTcond = 1;

          for (size_t col = 0; col < mineral_fraction.size(); ++col)
          {

           unsigned int mID = mineral_index[col];

           switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
           {
             case nondimoliv_index: // Dry Olivine
             {  
               double nondimoliv_latTcond = compute_lattice_thermal_conductivity_nondimn(
               nondimoliv_latTC_a0, 
               nondimoliv_latTC_b1, 
               nondimoliv_latTC_ymin, 
               nondimoliv_latTC_ymax,
               P_log, 
               T_ratio, 
               T_nondim, 
               nondimoliv_Tdep_n_exp); 
               double nondimoliv_radTcond = compute_radiative_thermal_conductivity_nondimn(
               nondimoliv_radTC_c0, 
               nondimoliv_radTC_d1, 
               nondimoliv_radTC_jmin, 
               nondimoliv_radTC_jmax, 
               T_log); 
               double nondimoliv_totTcond = compute_total_thermal_conductivity_nondimn(
               nondimoliv_latTcond, 
               nondimoliv_radTcond); 
               // Store the thermal conductivities in the vector
               nondim_minerals_latTcond[col] = nondimoliv_latTcond;
               nondim_minerals_radTcond[col] = nondimoliv_radTcond;
               nondim_minerals_totTcond[col] = nondimoliv_totTcond;
               break;
              }
            }

           // Thermal conductivity of the aggregate rock is computed as the
           // geometric mean of the total thermal conductivities of the minerals weighted by their fraction
           nondim_aggregate_rock_totTcond = nondim_aggregate_rock_totTcond * std::pow(nondim_minerals_totTcond[col], mineral_fraction[col]);
          
          }

          if (lithology != 99)
          {
             out.thermal_conductivities[i] = nondim_aggregate_rock_totTcond;
          }
          else if (lithology == 99)
             out.thermal_conductivities[i] = nondim_minerals_totTcond[0];
          else
          {
             AssertThrow(false, dealii::ExcMessage("Invalid lithology for the mantle."));
          }
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
      template class nondimensional_Tcond<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

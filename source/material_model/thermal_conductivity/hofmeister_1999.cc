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

#include <aspect/material_model/thermal_conductivity/hofmeister_1999.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Hofmeister (1999) formulations
namespace 
{
  double compute_lattice_thermal_conductivity_hof1999(const double latTC_room, const double n_Texp, const double gamma, const double alpha, const double K0, const double K_prime, const double T_room, const double T_model, const double P_model)
  {
   double factor_1 = latTC_room * std::pow((T_room / T_model), n_Texp);
   double factor_2 = std::exp(-((4*gamma)+(1.0/3.0)) * (alpha * (T_model - T_room)));
   double factor_3 = (1+((K_prime*P_model)/K0));
   double lattice_thermal_conductivity = factor_1 * factor_2 * factor_3;
   return lattice_thermal_conductivity;
  }

  double compute_radiative_thermal_conductivity_hof1999(const double a0, const double b1, const double c2, const double d3, const double T_model)
  {
   double radiative_thermal_conductivity = a0 + (b1 * T_model) + (c2 * std::pow(T_model, 2)) + (d3 * std::pow(T_model, 3));
   return radiative_thermal_conductivity;
  }

  double compute_total_thermal_conductivity_hof1999(const double lattice_thermal_conductivity, const double radiative_thermal_conductivity)
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
      hofmeister_1999<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        
        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int olivinedry_index = 0;
        const double olivinedry_hof99_latTC_room =  4.7;        // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double olivinedry_hof99_latTC_Texp =  0.3;        // temperature exponent (n_Texp)
        const double olivinedry_hof99_latTC_grun =  1.28;       // Grueneisen parameter(gamma)
        const double olivinedry_hof99_latTC_alph =  2.7500e-5;  // thermal expansion coefficient (alpha)
        const double olivinedry_hof99_latTC_Kzer =  128.1;      // bulk modulus (K0)
        const double olivinedry_hof99_latTC_Kpri =  4.6;        // pressure derivative of bulk modulus (K_prime)
        const double olivinedry_hof99_radTC_coa0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (a0)
        const double olivinedry_hof99_radTC_cob1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (b1)
        const double olivinedry_hof99_radTC_coc2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (c2)
        const double olivinedry_hof99_radTC_cod3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (d3)
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int wadsleydry_index = 1;
        const double wadsleydry_hof99_latTC_room =  7.7;        // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double wadsleydry_hof99_latTC_Texp =  0.3;        // temperature exponent (n_Texp)
        const double wadsleydry_hof99_latTC_grun =  1.00;       // Grueneisen parameter (gamma)
        const double wadsleydry_hof99_latTC_alph =  2.2821e-5;  // thermal expansion coefficient (alpha)
        const double wadsleydry_hof99_latTC_Kzer =  172;        // bulk modulus (K0)
        const double wadsleydry_hof99_latTC_Kpri =  4.8;        // pressure derivative of bulk modulus (K_prime)
        const double wadsleydry_hof99_radTC_coa0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (a0)
        const double wadsleydry_hof99_radTC_cob1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (b1)
        const double wadsleydry_hof99_radTC_coc2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (c2)
        const double wadsleydry_hof99_radTC_cod3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (d3)

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int ringwoodry_index = 2;
        const double ringwoodry_hof99_latTC_room =  7.7;        // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double ringwoodry_hof99_latTC_Texp =  0.3;        // temperature exponent (n_Texp)
        const double ringwoodry_hof99_latTC_grun =  1.25;       // Grueneisen parameter (gamma)
        const double ringwoodry_hof99_latTC_alph =  2.0535e-5;  // thermal expansion coefficient (alpha)
        const double ringwoodry_hof99_latTC_Kzer =  183;        // bulk modulus (K0)
        const double ringwoodry_hof99_latTC_Kpri =  5.2;        // pressure derivative of bulk modulus (K_prime)
        const double ringwoodry_hof99_radTC_coa0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (a0)
        const double ringwoodry_hof99_radTC_cob1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (b1)
        const double ringwoodry_hof99_radTC_coc2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (c2)
        const double ringwoodry_hof99_radTC_cod3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (d3)

        // Coefficients for Fe-bridgmanite (10%)
        // mineral chemical formula [Fe0.1Mg0.9SiO3]
        constexpr int brigma90Mg_index = 3;
        const double brigma90Mg_hof99_latTC_room =  7.7;        // lattice thermal conductivity at room P,T conditions (latTC_room)
        const double brigma90Mg_hof99_latTC_Texp =  0.3;        // temperature exponent (n_Texp)
        const double brigma90Mg_hof99_latTC_grun =  1.25;       // Grueneisen parameter (gamma)
        const double brigma90Mg_hof99_latTC_alph =  2.0535e-5;  // thermal expansion coefficient (alpha)
        const double brigma90Mg_hof99_latTC_Kzer =  183;        // bulk modulus (K0)
        const double brigma90Mg_hof99_latTC_Kpri =  5.2;        // pressure derivative of bulk modulus (K_prime)
        const double brigma90Mg_hof99_radTC_coa0 =  1.7530e-2;  // coefficient for radiative thermal conductivity T^0 (a0)
        const double brigma90Mg_hof99_radTC_cob1 = -1.0365e-4;  // coefficient for radiative thermal conductivity T^1 (b1)
        const double brigma90Mg_hof99_radTC_coc2 =  2.2451e-7;  // coefficient for radiative thermal conductivity T^2 (c2)
        const double brigma90Mg_hof99_radTC_cod3 = -3.4071e-11; // coefficient for radiative thermal conductivity T^3 (d3)

        unsigned int mineralpar_index = brigma90Mg_index+1; // Number of minerals

        // Dunite Upper Mantle (100% olivine)
        std::vector<double> minfract_hof99_duniteOl_UM = {1.00}; 
        std::vector<unsigned int> minindex_hof99_duniteOl_UM = {olivinedry_index};
        // Dunite Upper Mantle Transition Zone (100% wadsleyite)
        std::vector<double> minfract_hof99_duniteOl_UMTZ = {1.00};
        std::vector<unsigned int> minindex_hof99_duniteOl_UMTZ = {wadsleydry_index};
        // Dunite Lower Mantle Transition Zone (100% ringwoodite)
        std::vector<double> minfract_hof99_duniteOl_LMTZ = {1.00};
        std::vector<unsigned int> minindex_hof99_duniteOl_LMTZ = {ringwoodry_index};
        // Dunite Lower Mantle (100% bridgmanite)
        std::vector<double> minfract_hof99_duniteOl_LM = {1.00};
        std::vector<unsigned int> minindex_hof99_duniteOl_LM = {brigma90Mg_index};

        // Define room temperature [K] 
        const double T_room = 298.15; 
        
        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {

         // Preallocate a vector for storing thermal conductivities of minerals
         std::vector<double> hof99_minerals_latTcond(mineralpar_index, 0.0); // Lattice thermal conductivity
         std::vector<double> hof99_minerals_radTcond(mineralpar_index, 0.0); // Radiative thermal conductivity
         std::vector<double> hof99_minerals_totTcond(mineralpar_index, 0.0); // Total thermal conductivity
         // Preallocate a matrix for storing thermal conductivities of minerals
         std::vector<std::vector<double>> hof99_all_minerals_Tconds(mineralpar_index, std::vector<double>(3, 0.0));

         // Convert pressure unit from [Pa] to [GPa]
         double P_GPa = in.pressure[i]/1e9;
         // Take the temperature field of the model [K]
         double T_mod = in.temperature[i];
         // Take the mineral fraction of the model
         double min_frac = in.composition[0][i];
 
         unsigned int mID = in.Mineral_ID;

         switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
         {
           case olivinedry_index: // Dry Olivine
           { 
             double olivinedry_hof99_latTcon = compute_lattice_thermal_conductivity_hof1999(
             olivinedry_hof99_latTC_room, 
             olivinedry_hof99_latTC_Texp, 
             olivinedry_hof99_latTC_grun, 
             olivinedry_hof99_latTC_alph,
             olivinedry_hof99_latTC_Kzer, 
             olivinedry_hof99_latTC_Kpri, 
             T_room, 
             T_mod, 
             P_GPa);   
             double olivinedry_hof99_radTcon = compute_radiative_thermal_conductivity_hof1999(
             olivinedry_hof99_radTC_coa0, 
             olivinedry_hof99_radTC_cob1, 
             olivinedry_hof99_radTC_coc2, 
             olivinedry_hof99_radTC_cod3, 
             T_mod); 
             double olivinedry_hof99_totTcon = compute_total_thermal_conductivity_hof1999(
             olivinedry_hof99_latTcon, 
             olivinedry_hof99_radTcon); 
             // Store the thermal conductivities in the vector
             hof99_minerals_latTcond[olivinedry_index] = olivinedry_hof99_latTcon;
             hof99_minerals_radTcond[olivinedry_index] = olivinedry_hof99_radTcon;
             hof99_minerals_totTcond[olivinedry_index] = olivinedry_hof99_totTcon;
             break;
           }
           case wadsleydry_index: // Dry Wadsleyite 
           { 
             double wadsleydry_hof99_latTcon = compute_lattice_thermal_conductivity_hof1999(
             wadsleydry_hof99_latTC_room, 
             wadsleydry_hof99_latTC_Texp, 
             wadsleydry_hof99_latTC_grun, 
             wadsleydry_hof99_latTC_alph,
             wadsleydry_hof99_latTC_Kzer, 
             wadsleydry_hof99_latTC_Kpri, 
             T_room, 
             T_mod, 
             P_GPa); 
             double wadsleydry_hof99_radTcon = compute_radiative_thermal_conductivity_hof1999(
             wadsleydry_hof99_radTC_coa0,
             wadsleydry_hof99_radTC_cob1, 
             wadsleydry_hof99_radTC_coc2, 
             wadsleydry_hof99_radTC_cod3, 
             T_mod); 
             double wadsleydry_hof99_totTcon = compute_total_thermal_conductivity_hof1999(
             wadsleydry_hof99_latTcon, 
             wadsleydry_hof99_radTcon);
             // Store the thermal conductivities in the vector
             hof99_minerals_latTcond[wadsleydry_index] = wadsleydry_hof99_latTcon;
             hof99_minerals_radTcond[wadsleydry_index] = wadsleydry_hof99_radTcon;
             hof99_minerals_totTcond[wadsleydry_index] = wadsleydry_hof99_totTcon;
             break;
           }
           case ringwoodry_index: // Dry Ringwoodite
           { 
             double ringwoodry_hof99_latTcon = compute_lattice_thermal_conductivity_hof1999(
             ringwoodry_hof99_latTC_room, 
             ringwoodry_hof99_latTC_Texp, 
             ringwoodry_hof99_latTC_grun, 
             ringwoodry_hof99_latTC_alph,
             ringwoodry_hof99_latTC_Kzer, 
             ringwoodry_hof99_latTC_Kpri, 
             T_room, 
             T_mod, 
             P_GPa); 
             double ringwoodry_hof99_radTcon = compute_radiative_thermal_conductivity_hof1999(
             ringwoodry_hof99_radTC_coa0, 
             ringwoodry_hof99_radTC_cob1, 
             ringwoodry_hof99_radTC_coc2, 
             ringwoodry_hof99_radTC_cod3, 
             T_mod); 
             double ringwoodry_hof99_totTcon = compute_total_thermal_conductivity_hof1999(
             ringwoodry_hof99_latTcon, 
             ringwoodry_hof99_radTcon);
             // Store the thermal conductivities in the vector
             hof99_minerals_latTcond[ringwoodry_index] = ringwoodry_hof99_latTcon;
             hof99_minerals_radTcond[ringwoodry_index] = ringwoodry_hof99_radTcon;
             hof99_minerals_totTcond[ringwoodry_index] = ringwoodry_hof99_totTcon;
             break;
           }
         } 

         // Fill the matrix column by column
         for (unsigned int row = 0; row < mineralpar_index; ++row)
         {
           hof99_all_minerals_Tconds[row][0] = hof99_minerals_latTcond[row]; // Column 0: Lattice conductivities
           hof99_all_minerals_Tconds[row][1] = hof99_minerals_radTcond[row]; // Column 1: Radiative conductivities
           hof99_all_minerals_Tconds[row][2] = hof99_minerals_totTcond[row]; // Column 2: Total conductivities
         }

         // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
         double aggrock_testcase_hof99_Tcond = std::pow(hof99_all_minerals_Tconds[mID][2], min_frac); 

         // Test Case
         out.thermal_conductivities[i] = aggrock_testcase_hof99_Tcond;

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
      template class hofmeister_1999<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

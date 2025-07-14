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

#include <aspect/material_model/thermal_conductivity/anderson_1987.h>

// Helper functions in anonymous namespace to compute thermal conductivities using Anderson (1987) formulation
namespace 
{
  // Compute the lattice thermal conductivity using density 
  double compute_lattice_thermal_conductivity_and1987(double latTC_room, double densi_room, double densi_model, double k_Pexp)
  { 
    return latTC_room * std::pow((densi_model / densi_room), k_Pexp);
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
      anderson_1987<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                  MaterialModel::MaterialModelOutputs<dim> &out) const
     {
        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Coefficients for dry olivine 
        // mineral composition [(Mg1.8Fe0.2)SiO4]    
        constexpr int olivinedry_index = 0;
        const double olivinedry_and87_latTC_room =  3.6;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double olivinedry_and87_densi_room =  3329; // density at room P,T conditions [kg m^-3]
        const double olivinedry_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)
        
        // Coefficients for dry wadsleyite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int wadsleydry_index = 1;
        const double wadsleydry_and87_latTC_room =  5.9;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double wadsleydry_and87_densi_room =  3626; // density at room P,T conditions [kg m^-3]
        const double wadsleydry_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for dry ringwoodite
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int ringwoodry_index = 2;
        const double ringwoodry_and87_latTC_room =  4.9;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double ringwoodry_and87_densi_room =  3842; // density at room P,T conditions [kg m^-3]
        const double ringwoodry_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for Fe-bridgmanite (10%)
        // mineral composition [Fe0.1Mg0.9SiO3]
        constexpr int brigma90Mg_index = 3;
        const double brigma90Mg_and87_latTC_room =  3.8;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double brigma90Mg_and87_densi_room =  4262; // density at room P,T conditions [kg m^-3]
        const double brigma90Mg_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for orthopyroxene (enstatite)
        // mineral composition [Mg2Si2O6]
        constexpr int opxenstati_index = 4;
        const double opxenstati_and87_latTC_room =  5.8;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double opxenstati_and87_densi_room =  3304; // density at room P,T conditions [kg m^-3]
        const double opxenstati_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for clinopyroxene (diopside)
        // mineral composition [CaMgSi2O6]
        constexpr int cpxdiopsid_index = 5;
        const double cpxdiopsid_and87_latTC_room =  6.0;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double cpxdiopsid_and87_densi_room =  3288; // density at room P,T conditions [kg m^-3]
        const double cpxdiopsid_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for garnet (pyrope)
        // mineral composition [Mg3Al2Si3O12]
        constexpr int grtpyropes_index = 6;
        const double grtpyropes_and87_latTC_room =  4.4;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double grtpyropes_and87_densi_room =  3568; // density at room P,T conditions [kg m^-3]
        const double grtpyropes_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for garnet (majorite)
        // mineral composition [Mg3(MgSi)(SiO4)3]
        constexpr int grtmajorit_index = 7;
        const double grtmajorit_and87_latTC_room =  9.7;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double grtmajorit_and87_densi_room =  3477;  // density at room P,T conditions [kg m^-3]
        const double grtmajorit_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for davemaoite 
        // mineral composition [CaSiO3]
        constexpr int davemaoite_index = 8;
        const double davemaoite_and87_latTC_room =  10.9; // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double davemaoite_and87_densi_room =  4184; // density at room P,T conditions [kg m^-3]
        const double davemaoite_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for ferropericlase (10% Iron)
        // mineral composition [Mg0.90Fe0.10O] 
        constexpr int ferroper10_index = 9;
        const double ferroper10_and87_latTC_room =  3.5;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double ferroper10_and87_densi_room =  4006; // density at room P,T conditions [kg m^-3]
        const double ferroper10_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        unsigned int mineralpar_index = ferroper10_index+1; // Number of minerals

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> and87_minerals_latTcond(mineralpar_index, 0.0); // Lattice thermal conductivity
          std::vector<double> and87_minerals_totTcond(mineralpar_index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> and87_all_minerals_Tconds(mineralpar_index, std::vector<double>(3, 0.0));

          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];

          double densi_model = in.density[i];

          unsigned int mID = in.Mineral_ID;

          switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
          {
           case olivinedry_index: // Dry Olivine
            {      
             double olivinedry_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
             olivinedry_and87_latTC_room,
             olivinedry_and87_densi_room,
             densi_model,
             olivinedry_and87_k_Pexponen); 
             double olivinedry_and87_totTcon = olivinedry_and87_latTcon;
             // Store the thermal conductivities in the vector
             and87_minerals_latTcond[olivinedry_index] = olivinedry_and87_latTcon;
             and87_minerals_totTcond[olivinedry_index] = olivinedry_and87_totTcon;
             break;
            }
            case wadsleydry_index: // Dry Wadsleyite 
            {      
             double wadsleydry_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
             wadsleydry_and87_latTC_room,
             wadsleydry_and87_densi_room,
             densi_model,
             wadsleydry_and87_k_Pexponen); 
             double wadsleydry_and87_totTcon = wadsleydry_and87_latTcon;
             // Store the thermal conductivities in the vector
             and87_minerals_latTcond[wadsleydry_index] = wadsleydry_and87_latTcon;
             and87_minerals_totTcond[wadsleydry_index] = wadsleydry_and87_totTcon;
             break;
            }
            case ringwoodry_index: // Dry Ringwoodite
            {      
             double ringwoodry_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
             ringwoodry_and87_latTC_room,
             ringwoodry_and87_densi_room,
             densi_model,
             ringwoodry_and87_k_Pexponen); 
             double ringwoodry_and87_totTcon = ringwoodry_and87_latTcon;
             // Store the thermal conductivities in the vector
             and87_minerals_latTcond[ringwoodry_index] = ringwoodry_and87_latTcon;
             and87_minerals_totTcond[ringwoodry_index] = ringwoodry_and87_totTcon;
             break;
            }
            case brigma90Mg_index: // Fe-Bridgmanite (10%)
            {      
             double brigma90Mg_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
             brigma90Mg_and87_latTC_room,
             brigma90Mg_and87_densi_room,
             densi_model,
             brigma90Mg_and87_k_Pexponen); 
             double brigma90Mg_and87_totTcon = brigma90Mg_and87_latTcon;
             // Store the thermal conductivities in the vector
             and87_minerals_latTcond[brigma90Mg_index] = brigma90Mg_and87_latTcon;
             and87_minerals_totTcond[brigma90Mg_index] = brigma90Mg_and87_totTcon;
             break;
            }
            case opxenstati_index: // Orthopyroxene (Enstatite)
            {      
             double opxenstati_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
             opxenstati_and87_latTC_room,
             opxenstati_and87_densi_room,
             densi_model,
             opxenstati_and87_k_Pexponen); 
             double opxenstati_and87_totTcon = opxenstati_and87_latTcon;
             // Store the thermal conductivities in the vector
             and87_minerals_latTcond[opxenstati_index] = opxenstati_and87_latTcon;
             and87_minerals_totTcond[opxenstati_index] = opxenstati_and87_totTcon;
             break;
            }
            case cpxdiopsid_index: // Clinopyroxene (Diopside)
            { 
              double cpxdiopsid_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
              cpxdiopsid_and87_latTC_room,
              cpxdiopsid_and87_densi_room,
              densi_model,
              cpxdiopsid_and87_k_Pexponen); 
              double cpxdiopsid_and87_totTcon = cpxdiopsid_and87_latTcon;
              // Store the thermal conductivities in the vector
              and87_minerals_latTcond[cpxdiopsid_index] = cpxdiopsid_and87_latTcon;
              and87_minerals_totTcond[cpxdiopsid_index] = cpxdiopsid_and87_totTcon;
              break;
            }
            case grtpyropes_index: // Garnet (Pyrope)
            { 
              double grtpyropes_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
              grtpyropes_and87_latTC_room,
              grtpyropes_and87_densi_room,
              densi_model,
              grtpyropes_and87_k_Pexponen); 
              double grtpyropes_and87_totTcon = grtpyropes_and87_latTcon;
              // Store the thermal conductivities in the vector
              and87_minerals_latTcond[grtpyropes_index] = grtpyropes_and87_latTcon;
              and87_minerals_totTcond[grtpyropes_index] = grtpyropes_and87_totTcon;
              break;
            }
            case grtmajorit_index: // Garnet (Majorite)
            { 
              double grtmajorit_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
              grtmajorit_and87_latTC_room,
              grtmajorit_and87_densi_room,
              densi_model,
              grtmajorit_and87_k_Pexponen); 
              double grtmajorit_and87_totTcon = grtmajorit_and87_latTcon;
              // Store the thermal conductivities in the vector
              and87_minerals_latTcond[grtmajorit_index] = grtmajorit_and87_latTcon;
              and87_minerals_totTcond[grtmajorit_index] = grtmajorit_and87_totTcon;
              break;
            }
            case ferroper10_index: // Ferropericlase (Mg90Fe10O)
            { 
              double ferroper10_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
              ferroper10_and87_latTC_room,
              ferroper10_and87_densi_room,
              densi_model,
              ferroper10_and87_k_Pexponen); 
              double ferroper10_and87_totTcon = ferroper10_and87_latTcon;
              // Store the thermal conductivities in the vector
              and87_minerals_latTcond[ferroper10_index] = ferroper10_and87_latTcon;
              and87_minerals_totTcond[ferroper10_index] = ferroper10_and87_totTcon;
              break;
            }
            case davemaoite_index: // Davemaoite
            {
              double davemaoite_and87_latTcon = compute_lattice_thermal_conductivity_and1987(
              davemaoite_and87_latTC_room,
              davemaoite_and87_densi_room,
              densi_model,
              davemaoite_and87_k_Pexponen); 
              double davemaoite_and87_totTcon = davemaoite_and87_latTcon;
              // Store the thermal conductivities in the vector
              and87_minerals_latTcond[davemaoite_index] = davemaoite_and87_latTcon;
              and87_minerals_totTcond[davemaoite_index] = davemaoite_and87_totTcon;
              break;
            }
          }

          // Fill the matrix column by column
          for (unsigned int row = 0; row < mineralpar_index; ++row)
          {
            and87_all_minerals_Tconds[row][0] = and87_minerals_latTcond[row]; // Column 0: Lattice conductivities
            and87_all_minerals_Tconds[row][1] = and87_minerals_totTcond[row]; // Column 1: Total conductivities
          }

          // Aggregate rock thermal conductivity: geometric mean of the total thermal conductivities  of the minerals weighted by their fraction
          double aggrock_testcase_and87_Tcond = std::pow(and87_all_minerals_Tconds[mID][1], min_frac);

          // Test Case
          out.thermal_conductivities[i] = aggrock_testcase_and87_Tcond;
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
      template class anderson_1987<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

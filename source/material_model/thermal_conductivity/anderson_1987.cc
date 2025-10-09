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
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

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

        // Coefficients for ferropericlase (10% Iron)
        // mineral composition [Mg0.90Fe0.10O] 
        constexpr int ferroper10_index = 8;
        const double ferroper10_and87_latTC_room =  10.9;  // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double ferroper10_and87_densi_room =  4184; // density at room P,T conditions [kg m^-3]
        const double ferroper10_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // Coefficients for davemaoite 
        // mineral composition [CaSiO3]
        constexpr int davemaoite_index = 9;
        const double davemaoite_and87_latTC_room =  3.5; // lattice thermal conductivity at room P,T conditions (latTC_room) [W m^-1 K^-1]
        const double davemaoite_and87_densi_room =  4006; // density at room P,T conditions [kg m^-3]
        const double davemaoite_and87_k_Pexponen =  1.0;  // exponent for P-dependent thermal conductivity (k_Pexp)

        // unsigned int mineralpar_index = davemaoite_index+1; // Number of minerals

        // Preallocate a vector for mineral fractions of different rocks and a vector for mineral indices

        // pyrolite Upper Mantle (58% olivine, 13% pyrope, 17% ensatite, 12% diopside)
        std::vector<double> minfract_and87_pyrolite_UM = {0.58,0.13,0.17,0.12}; 
        std::vector<unsigned int> minindex_and87_pyrolite_UM = {olivinedry_index, grtpyropes_index, opxenstati_index, cpxdiopsid_index}; 
        // pyrolite Upper Mantle Transition Zone (58% wadsleyite, 28% majorite, 14% diopside)
        std::vector<double> minfract_and87_pyrolite_UMTZ = {0.58,0.28,0.14}; 
        std::vector<unsigned int> minindex_and87_pyrolite_UMTZ = {wadsleydry_index, grtmajorit_index, cpxdiopsid_index}; 
        // pyrolite Lower Mantle Transition Zone (58% ringwoodite, 42% majorite)
        std::vector<double> minfract_and87_pyrolite_LMTZ = {0.58,0.42}; 
        std::vector<unsigned int> minindex_and87_pyrolite_LMTZ = {ringwoodry_index, grtmajorit_index}; 
        // pyrolite Lower Mantle (80% bridgmanite, 14% ferropericlase, 6% davemaoite)
        std::vector<double> minfract_and87_pyrolite_LM = {0.80,0.14,0.06}; 
        std::vector<unsigned int> minindex_and87_pyrolite_LM = {brigma90Mg_index, ferroper10_index, davemaoite_index}; 

        // harzburgite Upper Mantle (80% olivine, 20% ensatite)
        std::vector<double> minfract_and87_harzburg_UM = {0.80,0.20}; 
        std::vector<unsigned int> minindex_and87_harzburg_UM = {olivinedry_index, opxenstati_index};
        // harzburgite Upper Mantle Transition Zone (80% wadsleyite, 13% diopside, 7% majorite)
        std::vector<double> minfract_and87_harzburg_UMTZ = {0.80,0.13,0.07}; 
        std::vector<unsigned int> minindex_and87_harzburg_UMTZ = {wadsleydry_index, cpxdiopsid_index, grtmajorit_index};
        // harzburgite Lower Mantle Transition Zone (80% ringwoodite, 20% majorite)
        std::vector<double> minfract_and87_harzburg_LMTZ = {0.80,0.20}; 
        std::vector<unsigned int> minindex_and87_harzburg_LMTZ = {ringwoodry_index, grtmajorit_index};
        // harzburgite Lower Mantle (76% bridgmanite, 24% ferropericlase)
        std::vector<double> minfract_and87_harzburg_LM = {0.76,0.24}; 
        std::vector<unsigned int> minindex_and87_harzburg_LM = {brigma90Mg_index, ferroper10_index};

        // Meta-basaltic crust MORB Upper Mantle (80% diopside, 20% pyrope)
        std::vector<double> minfract_and87_metaMORB_UM = {0.80,0.20}; 
        std::vector<unsigned int> minindex_and87_metaMORB_UM = {cpxdiopsid_index, grtpyropes_index};
        // Meta-basaltic crust MORB Upper Mantle Transition Zone (50% majorite, 50% diopside)
        std::vector<double> minfract_and87_metaMORB_UMTZ = {0.50,0.50}; 
        std::vector<unsigned int> minindex_and87_metaMORB_UMTZ = {grtmajorit_index, cpxdiopsid_index};
        // Meta-basaltic crust MORB Lower Mantle Transition Zone (100% majorite)
        std::vector<double> minfract_and87_metaMORB_LMTZ = {1.00}; 
        std::vector<unsigned int> minindex_and87_metaMORB_LMTZ = {grtmajorit_index};
        // Meta-basaltic crust MORB Lower Mantle (80% bridgmanite, 20% davemaoite)
        std::vector<double> minfract_and87_metaMORB_LM = {0.80,0.20}; 
        std::vector<unsigned int> minindex_and87_metaMORB_LM = {brigma90Mg_index, davemaoite_index};

        // Dunite Upper Mantle (100% olivine)
        std::vector<double> minfract_and87_duniteOl_UM = {1.00}; 
        std::vector<unsigned int> minindex_and87_duniteOl_UM = {olivinedry_index};
        // Dunite Upper Mantle Transition Zone (100% wadsleyite)
        std::vector<double> minfract_and87_duniteOl_UMTZ = {1.00};
        std::vector<unsigned int> minindex_and87_duniteOl_UMTZ = {wadsleydry_index};
        // Dunite Lower Mantle Transition Zone (100% ringwoodite)
        std::vector<double> minfract_and87_duniteOl_LMTZ = {1.00};
        std::vector<unsigned int> minindex_and87_duniteOl_LMTZ = {ringwoodry_index};
        // Dunite Lower Mantle (100% bridgmanite)
        std::vector<double> minfract_and87_duniteOl_LM = {1.00};
        std::vector<unsigned int> minindex_and87_duniteOl_LM = {brigma90Mg_index};

         // Check if the sum of Rock Mineral Fraction is equal to 1
        double sum_min_fract_and87_pyrolite_UM = std::accumulate(minfract_and87_pyrolite_UM.begin(), minfract_and87_pyrolite_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_pyrolite_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_pyrolite_UM must be equal to 1."));
        double sum_min_fract_and87_pyrolite_UMTZ = std::accumulate(minfract_and87_pyrolite_UMTZ.begin(), minfract_and87_pyrolite_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_pyrolite_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_pyrolite_UMTZ must be equal to 1."));
        double sum_min_fract_and87_pyrolite_LMTZ = std::accumulate(minfract_and87_pyrolite_LMTZ.begin(), minfract_and87_pyrolite_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_pyrolite_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_pyrolite_LMTZ must be equal to 1."));
        double sum_min_fract_and87_pyrolite_LM = std::accumulate(minfract_and87_pyrolite_LM.begin(), minfract_and87_pyrolite_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_pyrolite_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_pyrolite_LM must be equal to 1."));

        double sum_min_fract_and87_harzburg_UM = std::accumulate(minfract_and87_harzburg_UM.begin(), minfract_and87_harzburg_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_harzburg_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_harzburg_UM must be equal to 1."));
        double sum_min_fract_and87_harzburg_UMTZ = std::accumulate(minfract_and87_harzburg_UMTZ.begin(), minfract_and87_harzburg_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_harzburg_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_harzburg_UMTZ must be equal to 1."));
        double sum_min_fract_and87_harzburg_LMTZ = std::accumulate(minfract_and87_harzburg_LMTZ.begin(), minfract_and87_harzburg_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_harzburg_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_harzburg_LMTZ must be equal to 1."));
        double sum_min_fract_and87_harzburg_LM = std::accumulate(minfract_and87_harzburg_LM.begin(), minfract_and87_harzburg_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_harzburg_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_harzburg_LM must be equal to 1."));

        double sum_min_fract_and87_metaMORB_UM = std::accumulate(minfract_and87_metaMORB_UM.begin(), minfract_and87_metaMORB_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_metaMORB_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_metaMORB_UM must be equal to 1."));
        double sum_min_fract_and87_metaMORB_UMTZ = std::accumulate(minfract_and87_metaMORB_UMTZ.begin(), minfract_and87_metaMORB_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_metaMORB_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_metaMORB_UMTZ must be equal to 1."));
        double sum_min_fract_and87_metaMORB_LMTZ = std::accumulate(minfract_and87_metaMORB_LMTZ.begin(), minfract_and87_metaMORB_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_metaMORB_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_metaMORB_LMTZ must be equal to 1."));
        double sum_min_fract_and87_metaMORB_LM = std::accumulate(minfract_and87_metaMORB_LM.begin(), minfract_and87_metaMORB_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_metaMORB_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_metaMORB_LM must be equal to 1."));

        double sum_min_fract_and87_duniteOl_UM = std::accumulate(minfract_and87_duniteOl_UM.begin(), minfract_and87_duniteOl_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_duniteOl_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_duniteOl_UM must be equal to 1."));
        double sum_min_fract_and87_duniteOl_UMTZ = std::accumulate(minfract_and87_duniteOl_UMTZ.begin(), minfract_and87_duniteOl_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_duniteOl_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_duniteOl_UMTZ must be equal to 1."));
        double sum_min_fract_and87_duniteOl_LMTZ = std::accumulate(minfract_and87_duniteOl_LMTZ.begin(), minfract_and87_duniteOl_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_duniteOl_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_duniteOl_LMTZ must be equal to 1."));
        double sum_min_fract_and87_duniteOl_LM = std::accumulate(minfract_and87_duniteOl_LM.begin(), minfract_and87_duniteOl_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_and87_duniteOl_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_and87_duniteOl_LM must be equal to 1."));

        double density_UpMa_bot = 3725; // Upper Mantle bottom density in kg/m^3
        double density_UMTZ_top = 3725; // Upper Transition Zone top density in kg/m^3
        double density_UMTZ_bot = 3975; // Upper Transition Zone bottom density in kg/m^3
        double density_LMTZ_top = 3975; // Lower Transition Zone top density in kg/m^3
        double density_LMTZ_bot = 4000; // Lower Transition Zone bottom density in kg/m^3

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          double densi_model = in.density[i];

          double lithology = 0.0;
          // if there is a compositional field, use the first one as indicator for lithology
          if (in.composition[i].size() > 0)
          {
            lithology = in.composition[i][0];
          }

          std::vector<double> mineral_and87_fraction;    // Mineral fractions for the current lithology
          std::vector<unsigned int> mineral_and87_index; // Mineral indexes for the current lithology

          if (densi_model < density_UpMa_bot) // upper mantle
          {
            if (lithology == 0) // pyrolite
            {
              mineral_and87_fraction = minfract_and87_pyrolite_UM; 
              mineral_and87_index = minindex_and87_pyrolite_UM;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_and87_fraction = minfract_and87_harzburg_UM; 
              mineral_and87_index = minindex_and87_harzburg_UM;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_and87_fraction = minfract_and87_metaMORB_UM; 
              mineral_and87_index = minindex_and87_metaMORB_UM;
            }
            else if (lithology == 3) // dunite
            {
              mineral_and87_fraction = minfract_and87_duniteOl_UM; 
              mineral_and87_index = minindex_and87_duniteOl_UM;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_and87_fraction = {1.00}; 
              mineral_and87_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the upper mantle."));
            }
          }
          else if (densi_model >= density_UMTZ_top && densi_model <= density_UMTZ_bot) // upper transition zone
          {
            if (lithology == 0) // pyrolite
            {
              mineral_and87_fraction = minfract_and87_pyrolite_UMTZ; 
              mineral_and87_index = minindex_and87_pyrolite_UMTZ;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_and87_fraction = minfract_and87_harzburg_UMTZ; 
              mineral_and87_index = minindex_and87_harzburg_UMTZ;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_and87_fraction = minfract_and87_metaMORB_UMTZ; 
              mineral_and87_index = minindex_and87_metaMORB_UMTZ;
            }
            else if (lithology == 3) // dunite
            {
              mineral_and87_fraction = minfract_and87_duniteOl_UMTZ; 
              mineral_and87_index = minindex_and87_duniteOl_UMTZ;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_and87_fraction = {1.00};  
              mineral_and87_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the upper mantle transition zone."));
            }
          }
          else if (densi_model > density_LMTZ_top && densi_model <= density_LMTZ_bot) // lower transition zone
          {
            if (lithology == 0) // pyrolite
            {
              mineral_and87_fraction = minfract_and87_pyrolite_LMTZ; 
              mineral_and87_index = minindex_and87_pyrolite_LMTZ;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_and87_fraction = minfract_and87_harzburg_LMTZ; 
              mineral_and87_index = minindex_and87_harzburg_LMTZ;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_and87_fraction = minfract_and87_metaMORB_LMTZ; 
              mineral_and87_index = minindex_and87_metaMORB_LMTZ;
            }
            else if (lithology == 3) // dunite
            {
              mineral_and87_fraction = minfract_and87_duniteOl_LMTZ; 
              mineral_and87_index = minindex_and87_duniteOl_LMTZ;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_and87_fraction = {1.00}; 
              mineral_and87_index = {in.Mineral_ID};
            }
            else
            {
              AssertThrow(false, dealii::ExcMessage("Invalid lithology for the lower mantle transition zone."));
            }
          }
          else if (densi_model > density_LMTZ_bot) // lower mantle
          {
            if (lithology == 0) // pyrolite
            {
              mineral_and87_fraction = minfract_and87_pyrolite_LM; 
              mineral_and87_index = minindex_and87_pyrolite_LM;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_and87_fraction = minfract_and87_harzburg_LM; 
              mineral_and87_index = minindex_and87_harzburg_LM;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_and87_fraction = minfract_and87_metaMORB_LM; 
              mineral_and87_index = minindex_and87_metaMORB_LM;
            }
            else if (lithology == 3) // dunite
            {
              mineral_and87_fraction = minfract_and87_duniteOl_LM; 
              mineral_and87_index = minindex_and87_duniteOl_LM;
            }
            else if (lithology == 99) // test (all minerals)
            {
              mineral_and87_fraction = {1.00}; 
              mineral_and87_index = {in.Mineral_ID};
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
          std::vector<double> and87_minerals_latTcond(mineral_and87_fraction.size(), 0.0); // Lattice thermal conductivity
          std::vector<double> and87_minerals_totTcond(mineral_and87_fraction.size(), 0.0); // Total thermal conductivity

          // Preallocate total thermal conductivity of the aggregate rock
          double and87_aggregate_rock_totTcond = 1;

          for (size_t col = 0; col < mineral_and87_fraction.size(); ++col)
          {

           unsigned int mID = mineral_and87_index[col];

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
               and87_minerals_latTcond[col] = olivinedry_and87_latTcon;
               and87_minerals_totTcond[col] = olivinedry_and87_totTcon;
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
               and87_minerals_latTcond[col] = wadsleydry_and87_latTcon;
               and87_minerals_totTcond[col] = wadsleydry_and87_totTcon;
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
               and87_minerals_latTcond[col] = ringwoodry_and87_latTcon;
               and87_minerals_totTcond[col] = ringwoodry_and87_totTcon;
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
               and87_minerals_latTcond[col] = brigma90Mg_and87_latTcon;
               and87_minerals_totTcond[col] = brigma90Mg_and87_totTcon;
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
               and87_minerals_latTcond[col] = opxenstati_and87_latTcon;
               and87_minerals_totTcond[col] = opxenstati_and87_totTcon;
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
               and87_minerals_latTcond[col] = cpxdiopsid_and87_latTcon;
               and87_minerals_totTcond[col] = cpxdiopsid_and87_totTcon;
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
               and87_minerals_latTcond[col] = grtpyropes_and87_latTcon;
               and87_minerals_totTcond[col] = grtpyropes_and87_totTcon;
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
               and87_minerals_latTcond[col] = grtmajorit_and87_latTcon;
               and87_minerals_totTcond[col] = grtmajorit_and87_totTcon;
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
               and87_minerals_latTcond[col] = ferroper10_and87_latTcon;
               and87_minerals_totTcond[col] = ferroper10_and87_totTcon;
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
               and87_minerals_latTcond[col] = davemaoite_and87_latTcon;
               and87_minerals_totTcond[col] = davemaoite_and87_totTcon;
               break;
             }
            }

           // Thermal conductivity of the aggregate rock is computed as the
           // geometric mean of the total thermal conductivities of the minerals weighted by their fraction
           and87_aggregate_rock_totTcond = and87_aggregate_rock_totTcond * std::pow(and87_minerals_totTcond[col], mineral_and87_fraction[col]);

          }

          if (lithology != 99)
          {
             out.thermal_conductivities[i] = and87_aggregate_rock_totTcond;
          }
          else if (lithology == 99)
             out.thermal_conductivities[i] = and87_minerals_totTcond[0];
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
      template class anderson_1987<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

      #undef INSTANTIATE
    }
  }
}

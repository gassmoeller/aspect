/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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

#include "common.h"
#include <aspect/utilities.h>
#include <aspect/material_model/thermal_conductivity/anderson_1987.h>
#include <aspect/material_model/thermal_conductivity/gerya_2021.h>
#include <aspect/material_model/thermal_conductivity/grose_afonso_2019.h>
#include <aspect/material_model/thermal_conductivity/hofmeister_1999.h>
#include <aspect/material_model/thermal_conductivity/hofmeister_2005.h>
#include <aspect/material_model/thermal_conductivity/hofmeister_branlund_2015.h>
#include <aspect/material_model/thermal_conductivity/marzotto_2025.h>
#include <aspect/material_model/thermal_conductivity/nondimensional_Tcond.h>
#include <aspect/material_model/thermal_conductivity/stackhouse_2015.h>
#include <aspect/material_model/thermal_conductivity/tosi_2016.h>
#include <aspect/material_model/thermal_conductivity/xu_2004.h>

TEST_CASE("Utilities::weighted_p_norm_average")
{
  std::vector<double> weights = {1,1,2,2,3,3};
  std::vector<double> values = {6,5,4,3,2,1};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};
  std::vector<double> expected = {1., 1.59237, 1.6974 , 1.98895, 2.38899, 2.83333, 3.24037, 3.41824, 3.57872, 3.85347, 6. };

  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      INFO("check i=" << i << ": ");
      REQUIRE(aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) == Approx(expected[i]));
    }

}

TEST_CASE("Utilities:: P,T dependent thermal conductivity Marzotto et al., 2025")
{
  aspect::MaterialModel::ThermalConductivity::marzotto_2025<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1700, 1800, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 5e9, 16e9, 22e9, 100e9};
  in.pressure = pressures;

  // Assigning a lithology index to in.composition 
  //(0: pyrolite, 1: harzburgite, 2: meta-MORB, 3: dunite, 99: test)
  std::vector<std::vector<double>> lithologies = 
  {
    {0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1},
    {2, 2, 2, 2, 2},
    {3, 3, 3, 3, 3},
    {99, 99, 99, 99, 99}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K] of minerals
  constexpr int olivinedry_exptID = 0;
  std::vector<double> olivinedry_expt_totTcon(temperatures.size());
  constexpr int wadsleydry_exptID = 1;
  std::vector<double> wadsleydry_expt_totTcon(temperatures.size());
  constexpr int ringwoodry_exptID = 2;
  std::vector<double> ringwoodry_expt_totTcon(temperatures.size());
  constexpr int brigm100Mg_exptID = 3;
  std::vector<double> brigm100Mg_expt_totTcon(temperatures.size());
  constexpr int brigma97Mg_exptID = 4;
  std::vector<double> brigma97Mg_expt_totTcon(temperatures.size());
  constexpr int brigma90Mg_exptID = 5;
  std::vector<double> brigma90Mg_expt_totTcon(temperatures.size());
  constexpr int brigmaAlMg_exptID = 6;
  std::vector<double> brigmaAlMg_expt_totTcon(temperatures.size());
  constexpr int brigmaFeAl_exptID = 7;
  std::vector<double> brigmaFeAl_expt_totTcon(temperatures.size());
  constexpr int opxenstati_exptID = 8;
  std::vector<double> opxenstati_expt_totTcon(temperatures.size());
  constexpr int cpxdiopsid_exptID = 9;
  std::vector<double> cpxdiopsid_expt_totTcon(temperatures.size());
  constexpr int grtpyropes_exptID = 10;
  std::vector<double> grtpyropes_expt_totTcon(temperatures.size());
  constexpr int grtgrossul_exptID = 11;
  std::vector<double> grtgrossul_expt_totTcon(temperatures.size());
  constexpr int grtalmandi_exptID = 12;
  std::vector<double> grtalmandi_expt_totTcon(temperatures.size());
  constexpr int grtmajorit_exptID = 13;
  std::vector<double> grtmajorit_expt_totTcon(temperatures.size());
  constexpr int quartzpure_exptID = 14;
  std::vector<double> quartzpure_expt_totTcon(temperatures.size());
  constexpr int coesitSiO2_exptID = 15;
  std::vector<double> coesitSiO2_expt_totTcon(temperatures.size());
  constexpr int stishovite_exptID = 16;
  std::vector<double> stishovite_expt_totTcon(temperatures.size());
  constexpr int stisho05Al_exptID = 17;
  std::vector<double> stisho05Al_expt_totTcon(temperatures.size());
  constexpr int antigor010_exptID = 18;
  std::vector<double> antigor010_expt_totTcon(temperatures.size());
  constexpr int antigor001_exptID = 19;
  std::vector<double> antigor001_expt_totTcon(temperatures.size());
  constexpr int phaseDFeAl_exptID = 20;
  std::vector<double> phaseDFeAl_expt_totTcon(temperatures.size());
  constexpr int phaseD02Al_exptID = 21;
  std::vector<double> phaseD02Al_expt_totTcon(temperatures.size());
  constexpr int ferroper08_exptID = 22;
  std::vector<double> ferroper08_expt_totTcon(temperatures.size());
  constexpr int ferroper10_exptID = 23;
  std::vector<double> ferroper10_expt_totTcon(temperatures.size());
  constexpr int ferroper20_exptID = 24;
  std::vector<double> ferroper20_expt_totTcon(temperatures.size());
  constexpr int ferroper56_exptID = 25;
  std::vector<double> ferroper56_expt_totTcon(temperatures.size());
  constexpr int davemaoite_exptID = 26;
  std::vector<double> davemaoite_expt_totTcon(temperatures.size());
  constexpr int newhexAlph_exptID = 27;
  std::vector<double> newhexAlph_expt_totTcon(temperatures.size());
  constexpr int akimotoite_exptID = 28;
  std::vector<double> akimotoite_expt_totTcon(temperatures.size());

  // Preallocate the expected total thermal conductivities (k) in [W/m/K] of rocks
  constexpr int pyrolite_exptID = 0;
  std::vector<double> pyrolite_expt_totTcon(temperatures.size());
  constexpr int harzburg_exptID = 1;
  std::vector<double> harzburg_expt_totTcon(temperatures.size());
  constexpr int metaMORB_exptID = 2;
  std::vector<double> metaMORB_expt_totTcon(temperatures.size());
  constexpr int duniteOl_exptID = 3;
  std::vector<double> duniteOl_expt_totTcon(temperatures.size());
  // constexpr int testmine_exptID = 4;
  // std::vector<double> testmine_expt_totTcon(temperatures.size());

  unsigned int mineralpar_index = akimotoite_exptID+1; // Number of minerals
  unsigned int rockspar_index = duniteOl_exptID+1; // Number of rocks

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_minerals_latTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_minerals_radTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_minerals_totTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Total thermal conductivity

  // Preallocate matrixes for storing thermal conductivities of rocks
  // std::vector<std::vector<double>> expt_rocks_latTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Lattice thermal conductivity
  // std::vector<std::vector<double>> expt_rocks_radTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_rocks_totTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Total thermal conductivity

  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_latTcon = {3.588882835, 2.435539710, 4.726870069, 4.996780305, 4.257676002};
  std::vector<double> olivinedry_expt_radTcon = {0.00138288, 2.23152, 2.34978, 2.45486, 3.12273};
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> wadsleydry_expt_latTcon = {5.883641867, 3.313822545, 3.396102843, 3.349763783, 2.773684075};
  std::vector<double> wadsleydry_expt_radTcon = {1.1834e-9, 1.51515, 1.71067, 1.87995, 2.71457};
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> ringwoodry_expt_latTcon = {8.672339322, 4.374552957, 5.684237930, 6.320588146, 13.717405178};
  std::vector<double> ringwoodry_expt_radTcon = {5.52667e-10, 0.74367, 0.85004, 0.94227, 1.38818};
  // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigm100Mg_expt_latTcon = {10.695037531, 2.351863333, 3.145064580, 3.482255641, 5.690029119};
  std::vector<double> brigm100Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigma97Mg_expt_latTcon = {5.737509108, 2.574655637, 3.294977563, 3.635423378, 6.952419599};
  std::vector<double> brigma97Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigma90Mg_expt_latTcon = {3.791217949, 3.174583362, 3.861598953, 4.224518327, 8.920763745};
  std::vector<double> brigma90Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> brigmaAlMg_expt_latTcon = {6.304027763, 2.666956718, 3.582814931, 4.018317477, 7.886262751};
  std::vector<double> brigmaAlMg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigmaFeAl_expt_latTcon = {3.999620927, 2.112809941, 2.758962273, 3.081860257, 6.171981516};
  std::vector<double> brigmaFeAl_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> opxenstati_expt_latTcon = {5.799503276, 3.247391474, 3.392606680, 3.306155169, 2.566413033};
  std::vector<double> opxenstati_expt_radTcon = {5.45056e-5, 2.93307, 3.08161, 3.20925, 3.90772};
  // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> cpxdiopsid_expt_latTcon = {5.992731146, 3.234777477, 4.072304813, 4.127063506, 3.416805199};
  std::vector<double> cpxdiopsid_expt_radTcon = {2.64987e-5, 2.99817, 3.14622, 3.27243, 3.94085};
  // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> grtpyropes_expt_latTcon = {4.388274092, 2.716973203, 4.408565964, 4.691962737, 4.224051220};
  std::vector<double> grtpyropes_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> grtgrossul_expt_latTcon = {4.088378059, 2.329479982, 3.974765217, 4.338634179, 4.014125491};
  std::vector<double> grtgrossul_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> grtalmandi_expt_latTcon = {3.391236937, 2.235613596, 4.074648592, 4.416571698, 4.075429705};
  std::vector<double> grtalmandi_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
  std::vector<double> grtmajorit_expt_latTcon = {9.739829248, 4.711434237, 5.768834705, 5.833221258, 4.762499667};
  std::vector<double> grtmajorit_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> quartzpure_expt_latTcon = {9.532431770, 2.656822576, 2.647526672, 2.503780368, 1.493230149};
  std::vector<double> quartzpure_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> coesitSiO2_expt_latTcon = {7.211963110, 1.317885657, 1.243056895, 1.178804741, 0.849825723};
  std::vector<double> coesitSiO2_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> stishovite_expt_latTcon = {67.686173356, 29.308701033, 28.378068245, 27.409461557, 29.604731196};
  std::vector<double> stishovite_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> stisho05Al_expt_latTcon = {24.185714678, 10.678393756,11.357373198, 11.723480486, 15.115699610};
  std::vector<double> stisho05Al_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> antigor010_expt_latTcon = {4.555887376, 2.486216120, 3.908524818, 4.127554249, 3.574463165};
  std::vector<double> antigor010_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> antigor001_expt_latTcon = {1.066695732, 1.049619491, 1.787724820, 1.821134848, 1.485778351};
  std::vector<double> antigor001_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> phaseDFeAl_expt_latTcon = {2.593251373, 1.360157872, 1.918148669, 2.017880533, 6.126485225};
  std::vector<double> phaseDFeAl_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> phaseD02Al_expt_latTcon = {3.606657774, 1.702912781, 2.691714734, 3.471527545, 8.617758947};
  std::vector<double> phaseD02Al_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper08_expt_latTcon = {5.084250682, 2.31862542, 3.26811261, 4.140228226, 14.39940943};
  std::vector<double> ferroper08_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper10_expt_latTcon = {4.486103543, 2.041799599, 2.822614132, 3.535729625, 12.65637458};
  std::vector<double> ferroper10_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg80Fe20O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper20_expt_latTcon = {3.486472241, 3.569902514, 4.043245052, 4.297325425, 7.573244781};
  std::vector<double> ferroper20_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper56_expt_latTcon = {2.691665917, 1.59851609, 2.774573165, 3.378749002, 7.040686707};
  std::vector<double> ferroper56_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> davemaoite_expt_latTcon = {10.86634023, 5.311473757, 7.244743745, 8.242377128, 13.48442897};
  std::vector<double> davemaoite_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> newhexAlph_expt_latTcon = {10.59581461, 4.588122584, 4.45336532, 4.351526862, 12.15183764};
  std::vector<double> newhexAlph_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> akimotoite_expt_latTcon = {10.695037531, 2.351863333, 3.145064580, 3.482255641, 5.690029119};
  std::vector<double> akimotoite_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    olivinedry_expt_totTcon[row] = olivinedry_expt_latTcon[row] + olivinedry_expt_radTcon[row];
    expt_minerals_latTcond[row][olivinedry_exptID] = olivinedry_expt_latTcon[row];
    expt_minerals_radTcond[row][olivinedry_exptID] = olivinedry_expt_radTcon[row];  
    expt_minerals_totTcond[row][olivinedry_exptID] = olivinedry_expt_totTcon[row];
    AssertThrow(olivinedry_expt_totTcon[row] > 0, dealii::ExcMessage("olivinedry is < 0 ; Base for pow must be > 0"));
    // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    wadsleydry_expt_totTcon[row] = wadsleydry_expt_latTcon[row] + wadsleydry_expt_radTcon[row];
    expt_minerals_latTcond[row][wadsleydry_exptID] = wadsleydry_expt_latTcon[row];
    expt_minerals_radTcond[row][wadsleydry_exptID] = wadsleydry_expt_radTcon[row];   
    expt_minerals_totTcond[row][wadsleydry_exptID] = wadsleydry_expt_totTcon[row];
    AssertThrow(wadsleydry_expt_totTcon[row] > 0, dealii::ExcMessage("wadsleydry is < 0 ; Base for pow must be > 0"));
    // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    ringwoodry_expt_totTcon[row] = ringwoodry_expt_latTcon[row] + ringwoodry_expt_radTcon[row];
    expt_minerals_latTcond[row][ringwoodry_exptID] = ringwoodry_expt_latTcon[row];
    expt_minerals_radTcond[row][ringwoodry_exptID] = ringwoodry_expt_radTcon[row]; 
    expt_minerals_totTcond[row][ringwoodry_exptID] = ringwoodry_expt_totTcon[row];
    AssertThrow(ringwoodry_expt_totTcon[row] > 0, dealii::ExcMessage("ringwoodry is < 0 ; Base for pow must be > 0"));
    // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigm100Mg_expt_totTcon[row] = brigm100Mg_expt_latTcon[row] + brigm100Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_radTcon[row]; 
    expt_minerals_totTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_totTcon[row];
    AssertThrow(brigm100Mg_expt_totTcon[row] > 0, dealii::ExcMessage("brigm100Mg is < 0 ; Base for pow must be > 0"));
    // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigma97Mg_expt_totTcon[row] = brigma97Mg_expt_latTcon[row] + brigma97Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_radTcon[row];
    expt_minerals_totTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_totTcon[row];
    AssertThrow(brigma97Mg_expt_totTcon[row] > 0, dealii::ExcMessage("brigma97Mg is < 0 ; Base for pow must be > 0"));
    // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigma90Mg_expt_totTcon[row] = brigma90Mg_expt_latTcon[row] + brigma90Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_radTcon[row];
    expt_minerals_totTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_totTcon[row];
    AssertThrow(brigma90Mg_expt_totTcon[row] > 0, dealii::ExcMessage("brigma90Mg is < 0 ; Base for pow must be > 0"));
    // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    brigmaAlMg_expt_totTcon[row] = brigmaAlMg_expt_latTcon[row] + brigmaAlMg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_radTcon[row]; 
    expt_minerals_totTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_totTcon[row];
    AssertThrow(brigmaAlMg_expt_totTcon[row] > 0, dealii::ExcMessage("brigmaAlMg is < 0 ; Base for pow must be > 0"));
    // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigmaFeAl_expt_totTcon[row] = brigmaFeAl_expt_latTcon[row] + brigmaFeAl_expt_radTcon[row];
    expt_minerals_latTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_latTcon[row];
    expt_minerals_radTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_radTcon[row];
    expt_minerals_totTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_totTcon[row];
    AssertThrow(brigmaFeAl_expt_totTcon[row] > 0, dealii::ExcMessage("brigmaFeAl is < 0 ; Base for pow must be > 0"));
    // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    opxenstati_expt_totTcon[row] = opxenstati_expt_latTcon[row] + opxenstati_expt_radTcon[row];
    expt_minerals_latTcond[row][opxenstati_exptID] = opxenstati_expt_latTcon[row];
    expt_minerals_radTcond[row][opxenstati_exptID] = opxenstati_expt_radTcon[row];
    expt_minerals_totTcond[row][opxenstati_exptID] = opxenstati_expt_totTcon[row];
    AssertThrow(opxenstati_expt_totTcon[row] > 0, dealii::ExcMessage("opxenstati is < 0 ; Base for pow must be > 0"));
    // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    cpxdiopsid_expt_totTcon[row] = cpxdiopsid_expt_latTcon[row] + cpxdiopsid_expt_radTcon[row];
    expt_minerals_latTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_latTcon[row];
    expt_minerals_radTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_radTcon[row];
    expt_minerals_totTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_totTcon[row];
    AssertThrow(cpxdiopsid_expt_totTcon[row] > 0, dealii::ExcMessage("cpxdiopsid is < 0 ; Base for pow must be > 0"));
    // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    grtpyropes_expt_totTcon[row] = grtpyropes_expt_latTcon[row] + grtpyropes_expt_radTcon[row];
    expt_minerals_latTcond[row][grtpyropes_exptID] = grtpyropes_expt_latTcon[row];
    expt_minerals_radTcond[row][grtpyropes_exptID] = grtpyropes_expt_radTcon[row];
    expt_minerals_totTcond[row][grtpyropes_exptID] = grtpyropes_expt_totTcon[row];
    AssertThrow(grtpyropes_expt_totTcon[row] > 0, dealii::ExcMessage("grtpyropes is < 0 ; Base for pow must be > 0"));
    // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    grtgrossul_expt_totTcon[row] = grtgrossul_expt_latTcon[row] + grtgrossul_expt_radTcon[row];
    expt_minerals_latTcond[row][grtgrossul_exptID] = grtgrossul_expt_latTcon[row];
    expt_minerals_radTcond[row][grtgrossul_exptID] = grtgrossul_expt_radTcon[row];
    expt_minerals_totTcond[row][grtgrossul_exptID] = grtgrossul_expt_totTcon[row];
    AssertThrow(grtgrossul_expt_totTcon[row] > 0, dealii::ExcMessage("grtgrossul is < 0 ; Base for pow must be > 0"));
    // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    grtalmandi_expt_totTcon[row] = grtalmandi_expt_latTcon[row] + grtalmandi_expt_radTcon[row];
    expt_minerals_latTcond[row][grtalmandi_exptID] = grtalmandi_expt_latTcon[row];
    expt_minerals_radTcond[row][grtalmandi_exptID] = grtalmandi_expt_radTcon[row];
    expt_minerals_totTcond[row][grtalmandi_exptID] = grtalmandi_expt_totTcon[row];
    AssertThrow(grtalmandi_expt_totTcon[row] > 0, dealii::ExcMessage("grtalmandi is < 0 ; Base for pow must be > 0"));
    // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
    grtmajorit_expt_totTcon[row] = grtmajorit_expt_latTcon[row] + grtmajorit_expt_radTcon[row];
    expt_minerals_latTcond[row][grtmajorit_exptID] = grtmajorit_expt_latTcon[row];
    expt_minerals_radTcond[row][grtmajorit_exptID] = grtmajorit_expt_radTcon[row];   
    expt_minerals_totTcond[row][grtmajorit_exptID] = grtmajorit_expt_totTcon[row];
    AssertThrow(grtmajorit_expt_totTcon[row] > 0, dealii::ExcMessage("grtmajorit is < 0 ; Base for pow must be > 0"));
    // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    quartzpure_expt_totTcon[row] = quartzpure_expt_latTcon[row] + quartzpure_expt_radTcon[row];
    expt_minerals_latTcond[row][quartzpure_exptID] = quartzpure_expt_latTcon[row];
    expt_minerals_radTcond[row][quartzpure_exptID] = quartzpure_expt_radTcon[row];
    expt_minerals_totTcond[row][quartzpure_exptID] = quartzpure_expt_totTcon[row];
    AssertThrow(quartzpure_expt_totTcon[row] > 0, dealii::ExcMessage("quartzpure is < 0 ; Base for pow must be > 0"));
    // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    coesitSiO2_expt_totTcon[row] = coesitSiO2_expt_latTcon[row] + coesitSiO2_expt_radTcon[row];
    expt_minerals_latTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_latTcon[row];
    expt_minerals_radTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_radTcon[row];
    expt_minerals_totTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_totTcon[row];
    AssertThrow(coesitSiO2_expt_totTcon[row] > 0, dealii::ExcMessage("coesitSiO2 is < 0 ; Base for pow must be > 0"));
    // stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    stishovite_expt_totTcon[row] = stishovite_expt_latTcon[row] + stishovite_expt_radTcon[row];
    expt_minerals_latTcond[row][stishovite_exptID] = stishovite_expt_latTcon[row];
    expt_minerals_radTcond[row][stishovite_exptID] = stishovite_expt_radTcon[row];    
    expt_minerals_totTcond[row][stishovite_exptID] = stishovite_expt_totTcon[row];
    AssertThrow(stishovite_expt_totTcon[row] > 0, dealii::ExcMessage("stishovite is < 0 ; Base for pow must be > 0"));
    // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    stisho05Al_expt_totTcon[row] = stisho05Al_expt_latTcon[row] + stisho05Al_expt_radTcon[row];
    expt_minerals_latTcond[row][stisho05Al_exptID] = stisho05Al_expt_latTcon[row];
    expt_minerals_radTcond[row][stisho05Al_exptID] = stisho05Al_expt_radTcon[row];
    expt_minerals_totTcond[row][stisho05Al_exptID] = stisho05Al_expt_totTcon[row];
    AssertThrow(stisho05Al_expt_totTcon[row] > 0, dealii::ExcMessage("stisho05Al is < 0 ; Base for pow must be > 0"));
    // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    antigor010_expt_totTcon[row] = antigor010_expt_latTcon[row] + antigor010_expt_radTcon[row];
    expt_minerals_latTcond[row][antigor010_exptID] = antigor010_expt_latTcon[row];
    expt_minerals_radTcond[row][antigor010_exptID] = antigor010_expt_radTcon[row];
    expt_minerals_totTcond[row][antigor010_exptID] = antigor010_expt_totTcon[row];
    AssertThrow(antigor010_expt_totTcon[row] > 0, dealii::ExcMessage("antigor010 is < 0 ; Base for pow must be > 0"));
    // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    antigor001_expt_totTcon[row] = antigor001_expt_latTcon[row] + antigor001_expt_radTcon[row]; 
    expt_minerals_latTcond[row][antigor001_exptID] = antigor001_expt_latTcon[row];
    expt_minerals_radTcond[row][antigor001_exptID] = antigor001_expt_radTcon[row];
    expt_minerals_totTcond[row][antigor001_exptID] = antigor001_expt_totTcon[row];   
    AssertThrow(antigor001_expt_totTcon[row] > 0, dealii::ExcMessage("antigor001 is < 0 ; Base for pow must be > 0"));
    // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    phaseDFeAl_expt_totTcon[row] = phaseDFeAl_expt_latTcon[row] + phaseDFeAl_expt_radTcon[row];
    expt_minerals_latTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_latTcon[row];
    expt_minerals_radTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_radTcon[row];
    expt_minerals_totTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_totTcon[row]; 
    AssertThrow(phaseDFeAl_expt_totTcon[row] > 0, dealii::ExcMessage("phaseDFeAl is < 0 ; Base for pow must be > 0"));
    // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    phaseD02Al_expt_totTcon[row] = phaseD02Al_expt_latTcon[row] + phaseD02Al_expt_radTcon[row];
    expt_minerals_latTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_latTcon[row];
    expt_minerals_radTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_radTcon[row];
    expt_minerals_totTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_totTcon[row];
    AssertThrow(phaseD02Al_expt_totTcon[row] > 0, dealii::ExcMessage("phaseD02Al is < 0 ; Base for pow must be > 0"));
    // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper08_expt_totTcon[row] = ferroper08_expt_latTcon[row] + ferroper08_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper08_exptID] = ferroper08_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper08_exptID] = ferroper08_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper08_exptID] = ferroper08_expt_totTcon[row];
    AssertThrow(ferroper08_expt_totTcon[row] > 0, dealii::ExcMessage("ferroper08 is < 0 ; Base for pow must be > 0"));
    // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper10_expt_totTcon[row] = ferroper10_expt_latTcon[row] + ferroper10_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper10_exptID] = ferroper10_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper10_exptID] = ferroper10_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper10_exptID] = ferroper10_expt_totTcon[row];
    AssertThrow(ferroper10_expt_totTcon[row] > 0, dealii::ExcMessage("ferroper10 is < 0 ; Base for pow must be > 0"));
    // Ferropericlase (Mg80Fe20O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper20_expt_totTcon[row] = ferroper20_expt_latTcon[row] + ferroper20_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper20_exptID] = ferroper20_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper20_exptID] = ferroper20_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper20_exptID] = ferroper20_expt_totTcon[row];
    AssertThrow(ferroper20_expt_totTcon[row] > 0, dealii::ExcMessage("ferroper20 is < 0 ; Base for pow must be > 0"));
    // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper56_expt_totTcon[row] = ferroper56_expt_latTcon[row] + ferroper56_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper56_exptID] = ferroper56_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper56_exptID] = ferroper56_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper56_exptID] = ferroper56_expt_totTcon[row];
    AssertThrow(ferroper56_expt_totTcon[row] > 0, dealii::ExcMessage("ferroper56 is < 0 ; Base for pow must be > 0"));
    // davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    davemaoite_expt_totTcon[row] = davemaoite_expt_latTcon[row] + davemaoite_expt_radTcon[row]; 
    expt_minerals_latTcond[row][davemaoite_exptID] = davemaoite_expt_latTcon[row];
    expt_minerals_radTcond[row][davemaoite_exptID] = davemaoite_expt_radTcon[row];
    expt_minerals_totTcond[row][davemaoite_exptID] = davemaoite_expt_totTcon[row];
    AssertThrow(davemaoite_expt_totTcon[row] > 0, dealii::ExcMessage("davemaoite is < 0 ; Base for pow must be > 0"));
    // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    newhexAlph_expt_totTcon[row] = newhexAlph_expt_latTcon[row] + newhexAlph_expt_radTcon[row];
    expt_minerals_latTcond[row][newhexAlph_exptID] = newhexAlph_expt_latTcon[row];
    expt_minerals_radTcond[row][newhexAlph_exptID] = newhexAlph_expt_radTcon[row];
    expt_minerals_totTcond[row][newhexAlph_exptID] = newhexAlph_expt_totTcon[row];
    AssertThrow(newhexAlph_expt_totTcon[row] > 0, dealii::ExcMessage("newhexAlph is < 0 ; Base for pow must be > 0"));
    // akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    akimotoite_expt_totTcon[row] = akimotoite_expt_latTcon[row] + akimotoite_expt_radTcon[row];
    expt_minerals_latTcond[row][akimotoite_exptID] = akimotoite_expt_latTcon[row];
    expt_minerals_radTcond[row][akimotoite_exptID] = akimotoite_expt_radTcon[row]; 
    expt_minerals_totTcond[row][akimotoite_exptID] = akimotoite_expt_totTcon[row];
    AssertThrow(akimotoite_expt_totTcon[row] > 0, dealii::ExcMessage("akimotoite is < 0 ; Base for pow must be > 0"));
  }

  // Preallocate a vector for mineral fractions of different rocks and a vector for mineral indices
  // Pyrolite 
  std::vector<double> minfract_expt_pyrolite_UpMa = {0.58, 0.13, 0.17, 0.12};  // Upper Mantle (58% olivine, 13% pyrope, 17% ensatite, 12% diopside)
  std::vector<double> minfract_expt_pyrolite_UMTZ = {0.58, 0.28, 0.14};        // Upper Mantle Transition Zone (58% wadsleyite, 28% majorite, 14% diopside)  
  std::vector<double> minfract_expt_pyrolite_LMTZ = {0.58, 0.42};              // Lower Mantle Transition Zone (58% ringwoodite, 42% majorite)
  std::vector<double> minfract_expt_pyrolite_LoMa = {0.80, 0.14, 0.06};        // Lower Mantle (80% bridgmanite, 14% ferropericlase, 6% davemaoite)
  // Harzburgite 
  std::vector<double> minfract_expt_harzburg_UpMa = {0.80, 0.20};              // Upper Mantle (80% olivine, 20% ensatite)
  std::vector<double> minfract_expt_harzburg_UMTZ = {0.80, 0.13, 0.07};        // Upper Mantle Transition Zone (80% wadsleyite, 13% diopside, 7% majorite)
  std::vector<double> minfract_expt_harzburg_LMTZ = {0.80, 0.20};              // Lower Mantle Transition Zone (80% ringwoodite, 20% majorite)
  std::vector<double> minfract_expt_harzburg_LoMa = {0.76, 0.24};              // Lower Mantle (76% bridgmanite, 24% ferropericlase)
  // Meta-basalts (MORB)
  std::vector<double> minfract_expt_metaMORB_UpMa = {0.80, 0.20};              // Upper Mantle (80% diopside, 20% pyrope)
  std::vector<double> minfract_expt_metaMORB_UMTZ = {0.50, 0.04 ,0.46};        // Upper Mantle Transition Zone (50% majorite, 4% stishovite, 46% diopside)
  std::vector<double> minfract_expt_metaMORB_LMTZ = {0.92, 0.08};              // Lower Mantle Transition Zone (92% majorite, 8% stishovite)
  std::vector<double> minfract_expt_metaMORB_LoMa = {0.35, 0.28, 0.19, 0.18};  // Lower Mantle (35% bridgmanite, 28% davemaoite, 19% Fe-NAL, 18% stishovite)
  // Dunite (> 90% olivine)
  std::vector<double> minfract_expt_duniteOl_UpMa = {1.00};                    // Upper Mantle (100% olivine)
  std::vector<double> minfract_expt_duniteOl_UMTZ = {1.00};                    // Upper Mantle Transition Zone (100% wadsleyite)
  std::vector<double> minfract_expt_duniteOl_LMTZ = {1.00};                    // Lower Mantle Transition Zone (100% ringwoodite)
  std::vector<double> minfract_expt_duniteOl_LoMa = {1.00};                    // Lower Mantle (100% bridgmanite)

  // Check if the sum of Rock Mineral Fraction is equal to 1
  double sum_expt_fract_pyrolite_UpMa = std::accumulate(minfract_expt_pyrolite_UpMa.begin(), minfract_expt_pyrolite_UpMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_pyrolite_UpMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_pyrolite_UpMa must be equal to 1."));
  double sum_expt_fract_pyrolite_UMTZ = std::accumulate(minfract_expt_pyrolite_UMTZ.begin(), minfract_expt_pyrolite_UMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_pyrolite_UMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_pyrolite_UMTZ must be equal to 1."));
  double sum_expt_fract_pyrolite_LMTZ = std::accumulate(minfract_expt_pyrolite_LMTZ.begin(), minfract_expt_pyrolite_LMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_pyrolite_LMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_pyrolite_LMTZ must be equal to 1."));
  double sum_expt_fract_pyrolite_LoMa = std::accumulate(minfract_expt_pyrolite_LoMa.begin(), minfract_expt_pyrolite_LoMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_pyrolite_LoMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_pyrolite_LoMa must be equal to 1."));

  double sum_expt_fract_harzburg_UpMa = std::accumulate(minfract_expt_harzburg_UpMa.begin(), minfract_expt_harzburg_UpMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_harzburg_UpMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_harzburg_UpMa must be equal to 1."));
  double sum_expt_fract_harzburg_UMTZ = std::accumulate(minfract_expt_harzburg_UMTZ.begin(), minfract_expt_harzburg_UMTZ.end(), 0.0); 
  AssertThrow(std::abs(sum_expt_fract_harzburg_UMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_harzburg_UMTZ must be equal to 1."));
  double sum_expt_fract_harzburg_LMTZ = std::accumulate(minfract_expt_harzburg_LMTZ.begin(), minfract_expt_harzburg_LMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_harzburg_LMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_harzburg_LMTZ must be equal to 1."));
  double sum_expt_fract_harzburg_LoMa = std::accumulate(minfract_expt_harzburg_LoMa.begin(), minfract_expt_harzburg_LoMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_harzburg_LoMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_harzburg_LoMa must be equal to 1."));

  double sum_expt_fract_metaMORB_UpMa = std::accumulate(minfract_expt_metaMORB_UpMa.begin(), minfract_expt_metaMORB_UpMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_metaMORB_UpMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_metaMORB_UpMa must be equal to 1."));
  double sum_expt_fract_metaMORB_UMTZ = std::accumulate(minfract_expt_metaMORB_UMTZ.begin(), minfract_expt_metaMORB_UMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_metaMORB_UMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_metaMORB_UMTZ must be equal to 1."));
  double sum_expt_fract_metaMORB_LMTZ = std::accumulate(minfract_expt_metaMORB_LMTZ.begin(), minfract_expt_metaMORB_LMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_metaMORB_LMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_metaMORB_LMTZ must be equal to 1."));
  double sum_expt_fract_metaMORB_LoMa = std::accumulate(minfract_expt_metaMORB_LoMa.begin(), minfract_expt_metaMORB_LoMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_metaMORB_LoMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_metaMORB_LoMa must be equal to 1."));

  double sum_expt_fract_duniteOl_UpMa = std::accumulate(minfract_expt_duniteOl_UpMa.begin(), minfract_expt_duniteOl_UpMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_duniteOl_UpMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_duniteOl_UpMa must be equal to 1."));
  double sum_expt_fract_duniteOl_UMTZ = std::accumulate(minfract_expt_duniteOl_UMTZ.begin(), minfract_expt_duniteOl_UMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_duniteOl_UMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_duniteOl_UMTZ must be equal to 1."));
  double sum_expt_fract_duniteOl_LMTZ = std::accumulate(minfract_expt_duniteOl_LMTZ.begin(), minfract_expt_duniteOl_LMTZ.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_duniteOl_LMTZ - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_duniteOl_LMTZ must be equal to 1."));
  double sum_expt_fract_duniteOl_LoMa = std::accumulate(minfract_expt_duniteOl_LoMa.begin(), minfract_expt_duniteOl_LoMa.end(), 0.0);
  AssertThrow(std::abs(sum_expt_fract_duniteOl_LoMa - 1.0) < 1e-6, dealii::ExcMessage("Error: The sum of minfract_expt_duniteOl_LoMa must be equal to 1."));

  // Pyrolite
  std::vector<double> aggrock_pyrolite_UpMa_TCond(temperatures.size());
  std::vector<double> aggrock_pyrolite_UMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_pyrolite_LMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_pyrolite_LoMa_TCond(temperatures.size());
  // Harzburgite
  std::vector<double> aggrock_harzburg_UpMa_TCond(temperatures.size());
  std::vector<double> aggrock_harzburg_UMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_harzburg_LMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_harzburg_LoMa_TCond(temperatures.size());
  // Meta-basalts (MORB)
  std::vector<double> aggrock_metaMORB_UpMa_TCond(temperatures.size());
  std::vector<double> aggrock_metaMORB_UMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_metaMORB_LMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_metaMORB_LoMa_TCond(temperatures.size());
  // Dunite (> 90% olivine)
  std::vector<double> aggrock_duniteOl_UpMa_TCond(temperatures.size());
  std::vector<double> aggrock_duniteOl_UMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_duniteOl_LMTZ_TCond(temperatures.size());
  std::vector<double> aggrock_duniteOl_LoMa_TCond(temperatures.size());

  // Compute P,T-dependent thermal conductivities of aggregate rocks 
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    // Pyrolite
    aggrock_pyrolite_UpMa_TCond[row] = std::pow(olivinedry_expt_totTcon[row], minfract_expt_pyrolite_UpMa[0])*
                                       std::pow(grtpyropes_expt_totTcon[row], minfract_expt_pyrolite_UpMa[1])*
                                       std::pow(opxenstati_expt_totTcon[row], minfract_expt_pyrolite_UpMa[2])*
                                       std::pow(cpxdiopsid_expt_totTcon[row], minfract_expt_pyrolite_UpMa[3]);  
    aggrock_pyrolite_UMTZ_TCond[row] = std::pow(wadsleydry_expt_totTcon[row], minfract_expt_pyrolite_UMTZ[0])*
                                       std::pow(grtmajorit_expt_totTcon[row], minfract_expt_pyrolite_UMTZ[1])*
                                       std::pow(cpxdiopsid_expt_totTcon[row], minfract_expt_pyrolite_UMTZ[2]);  
    aggrock_pyrolite_LMTZ_TCond[row] = std::pow(ringwoodry_expt_totTcon[row], minfract_expt_pyrolite_LMTZ[0])*
                                       std::pow(grtmajorit_expt_totTcon[row], minfract_expt_pyrolite_LMTZ[1]); 
    aggrock_pyrolite_LoMa_TCond[row] = std::pow(brigmaAlMg_expt_totTcon[row], minfract_expt_pyrolite_LoMa[0])*
                                       std::pow(ferroper10_expt_totTcon[row], minfract_expt_pyrolite_LoMa[1])*
                                       std::pow(davemaoite_expt_totTcon[row], minfract_expt_pyrolite_LoMa[2]);
    AssertThrow(aggrock_pyrolite_UpMa_TCond[row] > 0, dealii::ExcMessage("aggrock_pyrolite_UpMa_TCond is <= 0")); 
    AssertThrow(aggrock_pyrolite_UMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_pyrolite_UMTZ_TCond is <= 0"));
    AssertThrow(aggrock_pyrolite_LMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_pyrolite_LMTZ_TCond is <= 0"));
    AssertThrow(aggrock_pyrolite_LoMa_TCond[row] > 0, dealii::ExcMessage("aggrock_pyrolite_LoMa_TCond is <= 0"));
    // Harzburgite
    aggrock_harzburg_UpMa_TCond[row] = std::pow(olivinedry_expt_totTcon[row], minfract_expt_harzburg_UpMa[0])*
                                       std::pow(opxenstati_expt_totTcon[row], minfract_expt_harzburg_UpMa[1]);
    aggrock_harzburg_UMTZ_TCond[row] = std::pow(wadsleydry_expt_totTcon[row], minfract_expt_harzburg_UMTZ[0])*
                                       std::pow(cpxdiopsid_expt_totTcon[row], minfract_expt_harzburg_UMTZ[1])*
                                       std::pow(grtmajorit_expt_totTcon[row], minfract_expt_harzburg_UMTZ[2]);
    aggrock_harzburg_LMTZ_TCond[row] = std::pow(ringwoodry_expt_totTcon[row], minfract_expt_harzburg_LMTZ[0])*
                                       std::pow(grtmajorit_expt_totTcon[row], minfract_expt_harzburg_LMTZ[1]);
    aggrock_harzburg_LoMa_TCond[row] = std::pow(brigmaAlMg_expt_totTcon[row], minfract_expt_harzburg_LoMa[0])*
                                       std::pow(ferroper10_expt_totTcon[row], minfract_expt_harzburg_LoMa[1]);
    AssertThrow(aggrock_harzburg_UpMa_TCond[row] > 0, dealii::ExcMessage("aggrock_harzburg_UpMa_TCond is <= 0")); 
    AssertThrow(aggrock_harzburg_UMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_harzburg_UMTZ_TCond is <= 0"));
    AssertThrow(aggrock_harzburg_LMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_harzburg_LMTZ_TCond is <= 0"));
    AssertThrow(aggrock_harzburg_LoMa_TCond[row] > 0, dealii::ExcMessage("aggrock_harzburg_LoMa_TCond is <= 0"));
    // Meta-basalts (MORB)
    aggrock_metaMORB_UpMa_TCond[row] = std::pow(cpxdiopsid_expt_totTcon[row], minfract_expt_metaMORB_UpMa[0])*
                                       std::pow(grtpyropes_expt_totTcon[row], minfract_expt_metaMORB_UpMa[1]);
    aggrock_metaMORB_UMTZ_TCond[row] = std::pow(grtmajorit_expt_totTcon[row], minfract_expt_metaMORB_UMTZ[0])*
                                       std::pow(stisho05Al_expt_totTcon[row], minfract_expt_metaMORB_UMTZ[1])*
                                       std::pow(cpxdiopsid_expt_totTcon[row], minfract_expt_metaMORB_UMTZ[2]);
    aggrock_metaMORB_LMTZ_TCond[row] = std::pow(grtmajorit_expt_totTcon[row], minfract_expt_metaMORB_LMTZ[0])*
                                       std::pow(stisho05Al_expt_totTcon[row], minfract_expt_metaMORB_LMTZ[1]);
    aggrock_metaMORB_LoMa_TCond[row] = std::pow(brigmaFeAl_expt_totTcon[row], minfract_expt_metaMORB_LoMa[0])*
                                       std::pow(davemaoite_expt_totTcon[row], minfract_expt_metaMORB_LoMa[1])*
                                       std::pow(newhexAlph_expt_totTcon[row], minfract_expt_metaMORB_LoMa[2])*
                                       std::pow(stisho05Al_expt_totTcon[row], minfract_expt_metaMORB_LoMa[3]);
    AssertThrow(aggrock_metaMORB_UpMa_TCond[row] > 0, dealii::ExcMessage("aggrock_metaMORB_UpMa_TCond is <= 0")); 
    AssertThrow(aggrock_metaMORB_UMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_metaMORB_UMTZ_TCond is <= 0"));
    AssertThrow(aggrock_metaMORB_LMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_metaMORB_LMTZ_TCond is <= 0"));
    AssertThrow(aggrock_metaMORB_LoMa_TCond[row] > 0, dealii::ExcMessage("aggrock_metaMORB_LoMa_TCond is <= 0"));
    // Dunite (> 90% olivine)
    aggrock_duniteOl_UpMa_TCond[row] = std::pow(olivinedry_expt_totTcon[row], minfract_expt_duniteOl_UpMa[0]);
    aggrock_duniteOl_UMTZ_TCond[row] = std::pow(wadsleydry_expt_totTcon[row], minfract_expt_duniteOl_UMTZ[0]);
    aggrock_duniteOl_LMTZ_TCond[row] = std::pow(ringwoodry_expt_totTcon[row], minfract_expt_duniteOl_LMTZ[0]);
    aggrock_duniteOl_LoMa_TCond[row] = std::pow(brigma90Mg_expt_totTcon[row], minfract_expt_duniteOl_LoMa[0]);
    AssertThrow(aggrock_duniteOl_UpMa_TCond[row] > 0, dealii::ExcMessage("aggrock_duniteOl_UpMa_TCond is <= 0")); 
    AssertThrow(aggrock_duniteOl_UMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_duniteOl_UMTZ_TCond is <= 0"));
    AssertThrow(aggrock_duniteOl_LMTZ_TCond[row] > 0, dealii::ExcMessage("aggrock_duniteOl_LMTZ_TCond is <= 0"));
    AssertThrow(aggrock_duniteOl_LoMa_TCond[row] > 0, dealii::ExcMessage("aggrock_duniteOl_LoMa_TCond is <= 0"));
  }
  
  double pressure_UpMa_top = 1e5;             // Upper Mantle top pressure in Pa
  double pressure_UpMa_bot = 13.59280101*1e9; // Upper Mantle bottom pressure in Pa
  double pressure_UMTZ_top = 13.59280101*1e9; // Upper Transition Zone top pressure in Pa
  double pressure_UMTZ_bot = 17.69264984*1e9; // Upper Transition Zone bottom pressure in Pa
  double pressure_LMTZ_top = 17.69264984*1e9; // Lower Transition Zone top pressure in Pa
  double pressure_LMTZ_bot = 23.1122152*1e9;  // Lower Transition Zone bottom pressure in Pa
  double pressure_LoMa_top = 23.1122152*1e9;  // Lower Mantle top pressure in Pa
  double pressure_LoMa_bot = 136*1e9;         // Lower Mantle bottom pressure in Pa


  // expected thermal conductivities (k) in [W/m/K] 
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    double current_pressure = in.pressure[row];

    if (current_pressure >= pressure_UpMa_top && current_pressure <= pressure_UpMa_bot) // Upper Mantle
    {     
       expt_rocks_totTcond[row][pyrolite_exptID] = aggrock_pyrolite_UpMa_TCond[row]; // Pyrolite
       expt_rocks_totTcond[row][harzburg_exptID] = aggrock_harzburg_UpMa_TCond[row]; // Harzburgite
       expt_rocks_totTcond[row][metaMORB_exptID] = aggrock_metaMORB_UpMa_TCond[row]; // Meta-basalts (MORB)
       expt_rocks_totTcond[row][duniteOl_exptID] = aggrock_duniteOl_UpMa_TCond[row]; // Dunite (> 90% olivine)
    }
    else if (current_pressure > pressure_UMTZ_top && current_pressure <= pressure_UMTZ_bot) // Upper Transition Zone
    {
       expt_rocks_totTcond[row][pyrolite_exptID] = aggrock_pyrolite_UMTZ_TCond[row]; // Pyrolite
       expt_rocks_totTcond[row][harzburg_exptID] = aggrock_harzburg_UMTZ_TCond[row]; // Harzburgite
       expt_rocks_totTcond[row][metaMORB_exptID] = aggrock_metaMORB_UMTZ_TCond[row]; // Meta-basalts (MORB)
       expt_rocks_totTcond[row][duniteOl_exptID] = aggrock_duniteOl_UMTZ_TCond[row]; // Dunite (> 90% olivine)
    }
    else if (current_pressure > pressure_LMTZ_top && current_pressure <= pressure_LMTZ_bot) // Lower Transition Zone
    {
       expt_rocks_totTcond[row][pyrolite_exptID] = aggrock_pyrolite_LMTZ_TCond[row]; // Pyrolite
       expt_rocks_totTcond[row][harzburg_exptID] = aggrock_harzburg_LMTZ_TCond[row]; // Harzburgite
       expt_rocks_totTcond[row][metaMORB_exptID] = aggrock_metaMORB_LMTZ_TCond[row]; // Meta-basalts (MORB)
       expt_rocks_totTcond[row][duniteOl_exptID] = aggrock_duniteOl_LMTZ_TCond[row]; // Dunite (> 90% olivine)
    }
    else if (current_pressure > pressure_LoMa_top && current_pressure <= pressure_LoMa_bot) // lower mantle
    {    
       expt_rocks_totTcond[row][pyrolite_exptID] = aggrock_pyrolite_LoMa_TCond[row]; // Pyrolite
       expt_rocks_totTcond[row][harzburg_exptID] = aggrock_harzburg_LoMa_TCond[row]; // Harzburgite
       expt_rocks_totTcond[row][metaMORB_exptID] = aggrock_metaMORB_LoMa_TCond[row]; // Meta-basalts (MORB)
       expt_rocks_totTcond[row][duniteOl_exptID] = aggrock_duniteOl_LoMa_TCond[row]; // Dunite (> 90% olivine)
    }
    else
    {
       AssertThrow(false, dealii::ExcMessage("Invalid pressure range for the mantle."));
    }
  }

  // Loop over all lithologies
  for (unsigned int lID = 0; lID < lithologies.size(); ++lID)
  {

    std::vector<double> current_lithology = lithologies[lID]; // Set the current lithology

    if (current_lithology[0] == 99)
    {

      INFO("Checking thermal conductivity (k) for different minerals as a function of temperature (T) and pressure (P)");

      // Loop over all mID values
      for (unsigned int mID = 0; mID < mineralpar_index; ++mID)
      {
        in.Mineral_ID = mID; // Set the current mID

        // Loop over the different combinations of pressures (P) and temperatures (T)
        for (size_t row = 0; row < temperatures.size(); ++row)
        {

         in.composition[row] = current_lithology;  // Assign the current lithology as model in.composition input

         model.evaluate(in, out);  // Call the function to compute the thermal conductivities

         switch (mID) // Compare the computed thermal conductivity with the expected value
         {
           case olivinedry_exptID: // olivinedry
           {
              INFO("Mineral: " << in.Mineral_ID << " (Dry Olivine) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("olivinedry expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("olivinedry computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case wadsleydry_exptID: // wadsleydry
            {
              INFO("Mineral: " << in.Mineral_ID << " (Dry Wadsleyite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("wadsleydry expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("wadsleydry computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case ringwoodry_exptID: // ringwoodry
            {
              INFO("Mineral: " << in.Mineral_ID << " (Dry Ringwoodite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("ringwoodry expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("ringwoodry computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case brigm100Mg_exptID: // brigm100Mg
            {
              INFO("Mineral: " << in.Mineral_ID << " (Bridgmanite 0% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("brigm100Mg expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("brigm100Mg computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case brigma97Mg_exptID: // brigma97Mg
            {
              INFO("Mineral: " << in.Mineral_ID << " (Bridgmanite 3% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("brigma97Mg expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("brigma97Mg computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case brigma90Mg_exptID: // brigma90Mg
            {
              INFO("Mineral: " << in.Mineral_ID << " (Bridgmanite 10% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("brigma90Mg expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("brigma90Mg computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case brigmaAlMg_exptID: // brigmaAlMg
            {
              INFO("Mineral: " << in.Mineral_ID << " (Bridgmanite Al) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("brigmaAlMg expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("brigmaAlMg computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case brigmaFeAl_exptID: // brigmaFeAl
            {
              INFO("Mineral: " << in.Mineral_ID << " (Bridgmanite Fe,Al) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("brigmaFeAl expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("brigmaFeAl computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case opxenstati_exptID: // opxenstati
            {
              INFO("Mineral: " << in.Mineral_ID << " (Enstatite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("opxenstati expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("opxenstati computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case cpxdiopsid_exptID: // cpxdiopsid
            {
              INFO("Mineral: " << in.Mineral_ID << " (Diopside) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("cpxdiopsid expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("cpxdiopsid computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case grtpyropes_exptID: // grtpyropes
            {
              INFO("Mineral: " << in.Mineral_ID << " (Pyrope) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("grtpyropes expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("grtpyropes computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case grtgrossul_exptID: // grtgrossul
            {
              INFO("Mineral: " << in.Mineral_ID << " (Grossular) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("grtgrossul expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("grtgrossul computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case grtalmandi_exptID: // grtalmandi
            {
              INFO("Mineral: " << in.Mineral_ID << " (Almandine) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("grtalmandi expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("grtalmandi computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case grtmajorit_exptID: // grtmajorit
            {
              INFO("Mineral: " << in.Mineral_ID << " (Majorite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("grtmajorit expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("grtmajorit computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case quartzpure_exptID: // quartzpure
            {
              INFO("Mineral: " << in.Mineral_ID << " (Quartz) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("quartzpure expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("quartzpure computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case coesitSiO2_exptID: // coesitSiO2
            {
              INFO("Mineral: " << in.Mineral_ID << " (Coesite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("coesitSiO2 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("coesitSiO2 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case stishovite_exptID: // stishovite
            {
              INFO("Mineral: " << in.Mineral_ID << " (Stishovite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("stishovite expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("stishovite computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case stisho05Al_exptID: // stisho05Al
            {
              INFO("Mineral: " << in.Mineral_ID << " (Stishovite 5% Al) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("stisho05Al expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("stisho05Al computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case antigor010_exptID: // antigor010
            {
              INFO("Mineral: " << in.Mineral_ID << " (Antigorite [010] direction) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("antigor010 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("antigor010 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case antigor001_exptID: // antigor001
            {
              INFO("Mineral: " << in.Mineral_ID << " (Antigorite [001] direction) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("antigor001 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("antigor001 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case phaseDFeAl_exptID: // phaseDFeAl
            {
              INFO("Mineral: " << in.Mineral_ID << " (Phase-D Fe,Al) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("phaseDFeAl expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("phaseDFeAl computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case phaseD02Al_exptID: // phaseD02Al
            {
              INFO("Mineral: " << in.Mineral_ID << " (Phase-D 2% Al) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("phaseD02Al expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("phaseD02Al computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case ferroper08_exptID: // ferroper08
            {
              INFO("Mineral: " << in.Mineral_ID << " (Ferropericlase 8% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("ferroper08 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("ferroper08 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case ferroper10_exptID: // ferroper10
            {
              INFO("Mineral: " << in.Mineral_ID << " (Ferropericlase 10% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("ferroper10 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("ferroper10 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case ferroper20_exptID: // ferroper20
            {
              INFO("Mineral: " << in.Mineral_ID << " (Ferropericlase 20% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("ferroper20 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("ferroper20 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case ferroper56_exptID: // ferroper56
            {
              INFO("Mineral: " << in.Mineral_ID << " (Ferropericlase 56% Fe) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("ferroper56 expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("ferroper56 computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case davemaoite_exptID: // davemaoite
            {
              INFO("Mineral: " << in.Mineral_ID << " (Davemaoite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("davemaoite expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("davemaoite computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case newhexAlph_exptID: // newhexAlph
            {
              INFO("Mineral: " << in.Mineral_ID << " (New Hexagonal Al-Phase) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("newhexAlph expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("newhexAlph computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
            case akimotoite_exptID: // akimotoite
            {
              INFO("Mineral: " << in.Mineral_ID << " (Akimotoite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
              INFO("akimotoite expected k= " << expt_minerals_totTcond[row][mID] << "[W/m/K]");
              INFO("akimotoite computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
              REQUIRE(out.thermal_conductivities[row] == Approx(expt_minerals_totTcond[row][mID]));
              break;
            }
          } 
        }
      }
    }
    else if (current_lithology[0] != 99) // To test the thermal conductivity of a all minerals
    {

     INFO("Checking thermal conductivity (k) for different lithologies as a function of temperature (T) and pressure (P)");

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t row = 0; row < temperatures.size(); ++row)
     {
       in.composition[row] = current_lithology;  // Assign the current lithology as model in.composition input

       model.evaluate(in, out);  // Call the function to compute the thermal conductivities

       switch (lID) // Compare the computed thermal conductivity with the expected value
       {
         case pyrolite_exptID: // Pyrolite
         {
           INFO("Lithology: " << in.composition[row][lID] << " (Pyrolite) ; Conditions: T = " << in.temperature[row] << "[K] ; P = " << (in.pressure[row]/1e9) << "[GPa]");
           INFO("Pyrolite expected k = " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
           INFO("Pyrolite computed k = " << out.thermal_conductivities[row] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
           break;
          }
         case harzburg_exptID: // Harzburgite
         {
           INFO("Lithology: " << in.composition[row][lID] << " (Harzburgite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
           INFO("Harzburgite expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
           INFO("Harzburgite computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
           break;
          }
         case metaMORB_exptID: // Metabasalt (MORB)
         {
           INFO("Lithology: " << in.composition[row][lID] << " (Metabasalt) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
           INFO("Metabasalt expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
           INFO("Metabasalt computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
           break;
          }
         case duniteOl_exptID: // Dunite (100% olivine)
         {
           INFO("Lithology: " << in.composition[row][lID] << " (Dunite) ; Conditions: T= " << in.temperature[row] << "[K] ; P= " << (in.pressure[row]/1e9) << "[GPa]");
           INFO("Dunite expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
           INFO("Dunite computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
           break;
          }
        }
      }
    }
    else
    {
       AssertThrow(false, dealii::ExcMessage("Invalid lithology for the mantle."));
    }
  }
}

TEST_CASE("Utilities:: Thermal Conductivity Hofmeister 1999")
{
  aspect::MaterialModel::ThermalConductivity::hofmeister_1999<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1700, 1800, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 1e9, 5e9, 10e9, 100e9};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };
    
  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int olivinedry_hof99_ID = 0;
  std::vector<double> olivinedry_expt_hof99_totTcon(temperatures.size());
  constexpr int wadsleydry_hof99_ID = 1;
  std::vector<double> wadsleydry_expt_hof99_totTcon(temperatures.size());
  constexpr int ringwoodry_hof99_ID = 2;
  std::vector<double> ringwoodry_expt_hof99_totTcon(temperatures.size());

  unsigned int hof99_index = ringwoodry_hof99_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_hof99_latTcond(hof99_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_hof99_radTcond(hof99_index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_hof99_totTcond(hof99_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_hof99_latTcon = {4.69000, 2.41949, 2.66506, 2.97359, 7.19819};
  std::vector<double> olivinedry_expt_hof99_radTcon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  expt_hof99_latTcond[olivinedry_hof99_ID] = olivinedry_expt_hof99_latTcon;
  expt_hof99_radTcond[olivinedry_hof99_ID] = olivinedry_expt_hof99_radTcon;
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> wadsleydry_expt_hof99_latTcon = {7.68434, 4.20366, 4.53114, 4.95031, 11.17803};
  std::vector<double> wadsleydry_expt_hof99_radTcon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  expt_hof99_latTcond[wadsleydry_hof99_ID] = wadsleydry_expt_hof99_latTcon;
  expt_hof99_radTcond[wadsleydry_hof99_ID] = wadsleydry_expt_hof99_radTcon;
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> ringwoodry_expt_hof99_latTcon = {7.68419, 4.14794, 4.47408, 4.89128, 11.00724};
  std::vector<double> ringwoodry_expt_hof99_radTcon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  expt_hof99_latTcond[ringwoodry_hof99_ID] = ringwoodry_expt_hof99_latTcon;
  expt_hof99_radTcond[ringwoodry_hof99_ID] = ringwoodry_expt_hof99_radTcon;

 // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    olivinedry_expt_hof99_totTcon[row] = olivinedry_expt_hof99_latTcon[row] + olivinedry_expt_hof99_radTcon[row];
    expt_hof99_totTcond[olivinedry_hof99_ID] = olivinedry_expt_hof99_totTcon;
    wadsleydry_expt_hof99_totTcon[row] = wadsleydry_expt_hof99_latTcon[row] + wadsleydry_expt_hof99_radTcon[row];
    expt_hof99_totTcond[wadsleydry_hof99_ID] = wadsleydry_expt_hof99_totTcon;
    ringwoodry_expt_hof99_totTcon[row] = ringwoodry_expt_hof99_latTcon[row] + ringwoodry_expt_hof99_radTcon[row];
    expt_hof99_totTcond[ringwoodry_hof99_ID] = ringwoodry_expt_hof99_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < hof99_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_hof99_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_hof99_Tcond[row][col] = std::pow(expt_hof99_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_hof99_Tcond;

   INFO("Checking hof99 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case olivinedry_hof99_ID: // olivinedry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case wadsleydry_hof99_ID: // wadsleydry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("wadsleydry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("wadsleydry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ringwoodry_hof99_ID: // ringwoodry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ringwoodry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ringwoodry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Anderson 1987")
{
  aspect::MaterialModel::ThermalConductivity::anderson_1987<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.densities in [Kg m^-3]
  std::vector<double> density = {3000, 3300, 3600, 3900, 4200};
  in.density = density;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
      {1.00, 1.00, 1.00, 1.00, 1.00},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.50, 0.50, 0.50, 0.50, 0.50},
      {0.25, 0.25, 0.25, 0.25, 0.25}
  };
    
  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int olivinedry_and87_ID = 0;
  std::vector<double> olivinedry_expt_and87_totTcon(density.size());
  constexpr int wadsleydry_and87_ID = 1;
  std::vector<double> wadsleydry_expt_and87_totTcon(density.size());
  constexpr int ringwoodry_and87_ID = 2;
  std::vector<double> ringwoodry_expt_and87_totTcon(density.size());
  constexpr int brigma90Mg_and87_ID = 3;
  std::vector<double> brigma90Mg_expt_and87_totTcon(density.size());
  constexpr int opxenstati_and87_ID = 4;
  std::vector<double> opxenstati_expt_and87_totTcon(density.size());
  constexpr int cpxdiopsid_and87_ID = 5;
  std::vector<double> cpxdiopsid_expt_and87_totTcon(density.size());
  constexpr int grtpyropes_and87_ID = 6;
  std::vector<double> grtpyropes_expt_and87_totTcon(density.size());
  constexpr int grtmajorit_and87_ID = 7;
  std::vector<double> grtmajorit_expt_and87_totTcon(density.size());
  constexpr int ferroper10_and87_ID = 8;
  std::vector<double> ferroper10_expt_and87_totTcon(density.size());
  constexpr int davemaoite_and87_ID = 9;
  std::vector<double> davemaoite_expt_and87_totTcon(density.size());

  unsigned int and87_index = davemaoite_and87_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_and87_latTcond(and87_index, std::vector<double>(density.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_and87_totTcond(and87_index, std::vector<double>(density.size(), 0.0)); // Total thermal conductivity

  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_and87_latTcon = {3.24422, 3.56864, 3.89306, 4.21748, 4.54190};
  expt_and87_latTcond[olivinedry_and87_ID] = olivinedry_expt_and87_latTcon;
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> wadsleydry_expt_and87_latTcon = {4.88141, 5.36955, 5.85769, 6.34584, 6.83398};
  expt_and87_latTcond[wadsleydry_and87_ID] = wadsleydry_expt_and87_latTcon;
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> ringwoodry_expt_and87_latTcon = {3.82613,	4.20875, 4.59136,	4.97397, 5.35659};
  expt_and87_latTcond[ringwoodry_and87_ID] = ringwoodry_expt_and87_latTcon;
  // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigma90Mg_expt_and87_latTcon = {2.67480,	2.94228, 3.20976,	3.47724, 3.74472};
  expt_and87_latTcond[brigma90Mg_and87_ID] = brigma90Mg_expt_and87_latTcon;
  // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> opxenstati_expt_and87_latTcon = {5.26634,	5.79298, 6.31961,	6.84625, 7.37288};
  expt_and87_latTcond[opxenstati_and87_ID] = opxenstati_expt_and87_latTcon;
  // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> cpxdiopsid_expt_and87_latTcon = {5.47445,	6.02190, 6.56934,	7.11679, 7.66423};
  expt_and87_latTcond[cpxdiopsid_and87_ID] = cpxdiopsid_expt_and87_latTcon;
  // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> grtpyropes_expt_and87_latTcon = {3.69955,	4.06951, 4.43946,	4.80942, 5.17937};
  expt_and87_latTcond[grtpyropes_and87_ID] = grtpyropes_expt_and87_latTcon;
  // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
  std::vector<double> grtmajorit_expt_and87_latTcon = {8.36928,	9.20621, 10.04314, 10.88007, 11.71700};
  expt_and87_latTcond[grtmajorit_and87_ID] = grtmajorit_expt_and87_latTcon;
  // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper10_expt_and87_latTcon = {7.81549,	8.59704, 9.37859,	10.16013,	10.94168};
  expt_and87_latTcond[ferroper10_and87_ID] = ferroper10_expt_and87_latTcon;
  // Davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> davemaoite_expt_and87_latTcon = {2.62107,	2.88318, 3.14528,	3.40739, 3.66950}; 
  expt_and87_latTcond[davemaoite_and87_ID] = davemaoite_expt_and87_latTcon;

 // Perform element-wise sum
  for (size_t row = 0; row < density.size(); ++row)
  {
    olivinedry_expt_and87_totTcon[row] = olivinedry_expt_and87_latTcon[row];
    expt_and87_totTcond[olivinedry_and87_ID] = olivinedry_expt_and87_totTcon;
    wadsleydry_expt_and87_totTcon[row] = wadsleydry_expt_and87_latTcon[row];
    expt_and87_totTcond[wadsleydry_and87_ID] = wadsleydry_expt_and87_totTcon;
    ringwoodry_expt_and87_totTcon[row] = ringwoodry_expt_and87_latTcon[row];
    expt_and87_totTcond[ringwoodry_and87_ID] = ringwoodry_expt_and87_totTcon;
    brigma90Mg_expt_and87_totTcon[row] = brigma90Mg_expt_and87_latTcon[row];
    expt_and87_totTcond[brigma90Mg_and87_ID] = brigma90Mg_expt_and87_totTcon;
    opxenstati_expt_and87_totTcon[row] = opxenstati_expt_and87_latTcon[row];
    expt_and87_totTcond[opxenstati_and87_ID] = opxenstati_expt_and87_totTcon;
    cpxdiopsid_expt_and87_totTcon[row] = cpxdiopsid_expt_and87_latTcon[row];
    expt_and87_totTcond[cpxdiopsid_and87_ID] = cpxdiopsid_expt_and87_totTcon;
    grtpyropes_expt_and87_totTcon[row] = grtpyropes_expt_and87_latTcon[row];
    expt_and87_totTcond[grtpyropes_and87_ID] = grtpyropes_expt_and87_totTcon;
    grtmajorit_expt_and87_totTcon[row] = grtmajorit_expt_and87_latTcon[row];
    expt_and87_totTcond[grtmajorit_and87_ID] = grtmajorit_expt_and87_totTcon;
    ferroper10_expt_and87_totTcon[row] = ferroper10_expt_and87_latTcon[row];
    expt_and87_totTcond[ferroper10_and87_ID] = ferroper10_expt_and87_totTcon;
    davemaoite_expt_and87_totTcon[row] = davemaoite_expt_and87_latTcon[row];
    expt_and87_totTcond[davemaoite_and87_ID] = davemaoite_expt_and87_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < and87_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_and87_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_and87_Tcond[row][col] = std::pow(expt_and87_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_and87_Tcond;

   INFO("Checking and87 thermal conductivity (k) for different densities (rho) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case olivinedry_and87_ID: // olivinedry
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case wadsleydry_and87_ID: // wadsleydry
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("wadsleydry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("wadsleydry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case ringwoodry_and87_ID: // ringwoodry
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ringwoodry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ringwoodry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case brigma90Mg_and87_ID: // brigmanite (10% Fe)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigma90Mg expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigma90Mg computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case opxenstati_and87_ID: // orthopyroxene (Enstatite)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("opxenstati expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("opxenstati computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case cpxdiopsid_and87_ID: // clinopyroxene (Diopside)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("cpxdiopsid expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("cpxdiopsid computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case grtpyropes_and87_ID: // garnet (Pyrope)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtpyropes expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtpyropes computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case grtmajorit_and87_ID: // garnet (Majorite)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtmajorit expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtmajorit computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
         case ferroper10_and87_ID: // ferropericlase (Mg90Fe10O)
         {
           INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ferroper10 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ferroper10 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
          }
          case davemaoite_and87_ID: // davemaoite
          {
            INFO("Conditions: rho= " << in.density[i] << "[kg m^-3] ; X= " << (in.composition[0][i])*100 << "[%]");
            INFO("davemaoite expected k= " << expected_conductivities[row][i] << "[W/m/K]");
            INFO("davemaoite computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
            REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
            break;
          }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Gerya 2021")
{
  aspect::MaterialModel::ThermalConductivity::gerya_2021<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 600, 900, 1200, 1500};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 1e6, 5e6, 10e6, 100e6};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
      {1.00, 1.00, 1.00, 1.00, 1.00},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.50, 0.50, 0.50, 0.50, 0.50},
      {0.25, 0.25, 0.25, 0.25, 0.25}
  };
    
  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int oceanicrus_ger21_ID = 0;
  std::vector<double> oceanicrus_expt_ger21_totTcon(temperatures.size());
  constexpr int lithomantl_ger21_ID = 1;
  std::vector<double> lithomantl_expt_ger21_totTcon(temperatures.size());
  constexpr int asthemantl_ger21_ID = 2;
  std::vector<double> asthemantl_expt_ger21_totTcon(temperatures.size());

  unsigned int ger21_index = asthemantl_ger21_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_ger21_latTcond(ger21_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_ger21_totTcond(ger21_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Oceanic Crust: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> oceanicrus_expt_ger21_latTcon = {2.43729, 1.88018, 1.66526, 1.55133, 1.48178};
  expt_ger21_latTcond[oceanicrus_ger21_ID] = oceanicrus_expt_ger21_latTcon;
  // Lithospheric Mantle: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> lithomantl_expt_ger21_latTcon = {4.15972, 2.63997, 2.05370, 1.74293, 1.55319};
  expt_ger21_latTcond[lithomantl_ger21_ID] = lithomantl_expt_ger21_latTcon;
  // Asthenospheric Mantle: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> asthemantl_expt_ger21_latTcon = {4.15972, 2.63997, 2.05370, 1.74293, 1.55319};  
  expt_ger21_latTcond[asthemantl_ger21_ID] = asthemantl_expt_ger21_latTcon;

 // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    oceanicrus_expt_ger21_totTcon[row] = oceanicrus_expt_ger21_latTcon[row];
    expt_ger21_totTcond[oceanicrus_ger21_ID] = oceanicrus_expt_ger21_totTcon;
    lithomantl_expt_ger21_totTcon[row] = lithomantl_expt_ger21_latTcon[row];
    expt_ger21_totTcond[lithomantl_ger21_ID] = lithomantl_expt_ger21_totTcon;
    asthemantl_expt_ger21_totTcon[row] = asthemantl_expt_ger21_latTcon[row];
    expt_ger21_totTcond[asthemantl_ger21_ID] = asthemantl_expt_ger21_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < ger21_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_ger21_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_ger21_Tcond[row][col] = std::pow(expt_ger21_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_ger21_Tcond;

   INFO("Checking ger21 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case oceanicrus_ger21_ID: // Oceanic Crust
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("oceanicrus expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("oceanicrus computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case lithomantl_ger21_ID: // Lithospheric Mantle
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("lithomantl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("lithomantl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case asthemantl_ger21_ID: // Asthenospheric Mantle
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("asthemantl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("asthemantl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Grose & Afonso 2019")
{
  aspect::MaterialModel::ThermalConductivity::grose_afonso_2019<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 600, 900, 1200, 1500};
  in.temperature = temperatures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int olivinedry_gra19_ID = 0;
  constexpr int opxenstati_gra19_ID = 1;
  constexpr int cpxdiopsid_gra19_ID = 2;
  constexpr int grtpyropes_gra19_ID = 3;

  unsigned int gra19_index = grtpyropes_gra19_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_gra19_radTcond(gra19_index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  
  // Dry Olivine: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_gra19_radTcon = {0.133625, 0.46500, 2.18212, 2.97168, 2.63066};
  expt_gra19_radTcond[olivinedry_gra19_ID] = olivinedry_expt_gra19_radTcon;
  // Orthopyroxene: expected lattice thermal conductivities (k) in [W/m/K]
  std::vector<double> opxenstati_expt_gra19_radTcon = {0.00653284, 0.173187, 0.84408, 2.11347, 4.08815};
  expt_gra19_radTcond[opxenstati_gra19_ID] = opxenstati_expt_gra19_radTcon;
  // Clinopyroxene: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> cpxdiopsid_expt_gra19_radTcon = {0.00270542, 0.146555, 0.823236, 2.28996, 4.61226};
  expt_gra19_radTcond[cpxdiopsid_gra19_ID] = cpxdiopsid_expt_gra19_radTcon;
  // Garnet: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> grtpyropes_expt_gra19_radTcon = {0.0247772, 0.0746832, 0.31400, 0.74358, 1.93803};
  expt_gra19_radTcond[grtpyropes_gra19_ID] = grtpyropes_expt_gra19_radTcon;

  // Loop over all mID values
  for (unsigned int mID = 0; mID < gra19_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_gra19_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_gra19_Tcond[row][col] = std::pow(expt_gra19_radTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_gra19_Tcond;

   INFO("Checking grose_afonso_2019 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case olivinedry_gra19_ID: // Dry Olivine
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case opxenstati_gra19_ID: // Orthopyroxene
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("opxenstati expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("opxenstati computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case cpxdiopsid_gra19_ID: // Clinopyroxene
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ;  X= " << (in.composition[0][i])*100 << "[%]");
           INFO("cpxdiopsid expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("cpxdiopsid computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case grtpyropes_gra19_ID: // Garnet
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtpyropes expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtpyropes computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Hofmeister 2005")
{
  aspect::MaterialModel::ThermalConductivity::hofmeister_2005<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 600, 900, 1200, 1500};
  in.temperature = temperatures;

  // Assigning an array of values to in.grainsize (d) in [cm]
  std::vector<double> grainsize = {0.1, 0.5, 1.0, 2.0, 50};
  in.grainsize = grainsize;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int olivinedry_hof05_ID = 0;

  unsigned int hof05_index = olivinedry_hof05_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_hof05_radTcond(hof05_index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  
  // Dry Olivine: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_hof05_radTcon = {0.0500328, 0.131068809, 0.625395051, 1.504829428, 1.534699e-4};
  expt_hof05_radTcond[olivinedry_hof05_ID] = olivinedry_expt_hof05_radTcon;

  // Loop over all mID values
  for (unsigned int mID = 0; mID < hof05_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_hof05_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_hof05_Tcond[row][col] = std::pow(expt_hof05_radTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_hof05_Tcond;

   INFO("Checking hofmeister_2005 thermal conductivity (k) for different temperatures (T), grain size (d) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case olivinedry_hof05_ID: // Dry Olivine
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; d= " << in.grainsize[i] << "[cm] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        }        
      }
    }
  }
}

TEST_CASE("Utilities:: Thermal Conductivity Hofmeister & Branlund 2015")
{
  aspect::MaterialModel::ThermalConductivity::hofmeister_branlund_2015<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 600, 900, 1200, 1500};
  in.temperature = temperatures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
      {1.00, 1.00, 1.00, 1.00, 1.00},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.50, 0.50, 0.50, 0.50, 0.50},
      {0.25, 0.25, 0.25, 0.25, 0.25}
  };
    
  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int oceanicrus_hbr15_ID = 0;
  std::vector<double> oceanicrus_expt_hbr15_totTcon(temperatures.size());

  unsigned int hbr15_index = oceanicrus_hbr15_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_hbr15_latTcond(hbr15_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_hbr15_totTcond(hbr15_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Oceanic Crust: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> oceanicrus_expt_hbr15_latTcon = {2.57469, 1.89585, 1.66956, 1.55642, 1.48854}; 
  expt_hbr15_latTcond[oceanicrus_hbr15_ID] = oceanicrus_expt_hbr15_latTcon;

 // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    oceanicrus_expt_hbr15_totTcon[row] = oceanicrus_expt_hbr15_latTcon[row];
    expt_hbr15_totTcond[oceanicrus_hbr15_ID] = oceanicrus_expt_hbr15_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < hbr15_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_hbr15_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_hbr15_Tcond[row][col] = std::pow(expt_hbr15_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_hbr15_Tcond;

   INFO("Checking hbr15 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case oceanicrus_hbr15_ID: // Oceanic Crust
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("oceanicrus expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("oceanicrus computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Stackhouse 2015")
{
  aspect::MaterialModel::ThermalConductivity::stackhouse_2015<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1800, 2000, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 30e9, 60e9, 90e9, 120e9};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
      {1.00, 1.00, 1.00, 1.00, 1.00},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.50, 0.50, 0.50, 0.50, 0.50},
      {0.25, 0.25, 0.25, 0.25, 0.25}
  };
    
  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int lowermantl_sta15_ID = 0;
  std::vector<double> lowermantl_expt_sta15_totTcon(temperatures.size());

  unsigned int sta15_index = lowermantl_sta15_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_sta15_latTcond(sta15_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_sta15_totTcond(sta15_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Lower Mantle Assemblage: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> lowermantl_expt_sta15_latTcon = {1.08579, 3.38751, 4.95177, 6.63374, 9.63497};
  expt_sta15_latTcond[lowermantl_sta15_ID] = lowermantl_expt_sta15_latTcon;

 // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    lowermantl_expt_sta15_totTcon[row] = lowermantl_expt_sta15_latTcon[row];
    expt_sta15_totTcond[lowermantl_sta15_ID] = lowermantl_expt_sta15_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < sta15_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_sta15_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_sta15_Tcond[row][col] = std::pow(expt_sta15_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_sta15_Tcond;

   INFO("Checking sta15 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case lowermantl_sta15_ID: // Lower Mantle Assemblage
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("lowermantl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("lowermantl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}

TEST_CASE("Utilities:: Thermal Conductivity Tosi 2016")
{
  aspect::MaterialModel::ThermalConductivity::tosi_2016<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1700, 1800, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 1e9, 5e9, 10e9, 100e9};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int uppermantl_tos16_ID = 0;
  std::vector<double> uppermantl_expt_tos16_totTcon(temperatures.size());
  constexpr int umantrazon_tos16_ID = 1;
  std::vector<double> umantrazon_expt_tos16_totTcon(temperatures.size());
  constexpr int lmantrazon_tos16_ID = 2;
  std::vector<double> lmantrazon_expt_tos16_totTcon(temperatures.size());
  constexpr int lowermantl_tos16_ID = 3;
  std::vector<double> lowermantl_expt_tos16_totTcon(temperatures.size());

  unsigned int tos16_index = lowermantl_tos16_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_tos16_latTcond(tos16_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_tos16_totTcond(tos16_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Upper Mantle: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> uppermantl_expt_tos16_latTcon = {2.46271, 1.24999, 1.78653, 2.43429, 11.71039};
  expt_tos16_latTcond[uppermantl_tos16_ID] = uppermantl_expt_tos16_latTcon;
  // Upper Mantle Transition Zone: expected lattice thermal conductivities (k) in [W/m/K]
  std::vector<double> umantrazon_expt_tos16_latTcon = {3.79686, 1.61966, 2.07866, 2.63431, 10.37773};
  expt_tos16_latTcond[umantrazon_tos16_ID] = umantrazon_expt_tos16_latTcon;
  // Lower Mantle Transition Zone: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> lmantrazon_expt_tos16_latTcon = {3.50678, 1.39227, 1.83968, 2.37777, 9.66447};
  expt_tos16_latTcond[lmantrazon_tos16_ID] = lmantrazon_expt_tos16_latTcon;
  // Lower Mantle: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> lowermantl_expt_tos16_latTcon = {3.47335, 2.13845, 2.37846, 2.68032, 7.56725};
  expt_tos16_latTcond[lowermantl_tos16_ID] = lowermantl_expt_tos16_latTcon;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    uppermantl_expt_tos16_totTcon[row] = uppermantl_expt_tos16_latTcon[row];
    expt_tos16_totTcond[uppermantl_tos16_ID] = uppermantl_expt_tos16_totTcon;
    umantrazon_expt_tos16_totTcon[row] = umantrazon_expt_tos16_latTcon[row];
    expt_tos16_totTcond[umantrazon_tos16_ID] = umantrazon_expt_tos16_totTcon;
    lmantrazon_expt_tos16_totTcon[row] = lmantrazon_expt_tos16_latTcon[row];
    expt_tos16_totTcond[lmantrazon_tos16_ID] = lmantrazon_expt_tos16_totTcon;
    lowermantl_expt_tos16_totTcon[row] = lowermantl_expt_tos16_latTcon[row];
    expt_tos16_totTcond[lowermantl_tos16_ID] = lowermantl_expt_tos16_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < tos16_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_tos16_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_tos16_Tcond[row][col] = std::pow(expt_tos16_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_tos16_Tcond;

   INFO("Checking tosi_2016 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case uppermantl_tos16_ID: // Upper Mantle (UM)
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("uppermantl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("uppermantl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case umantrazon_tos16_ID: // Upper Mantle Transition Zone (uMTZ)
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("umantrazon expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("umantrazon computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case lmantrazon_tos16_ID: // Lower Mantle Transition Zone (lMTZ)
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("lmantrazon expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("lmantrazon computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case lowermantl_tos16_ID: // Lower Mantle (LM)
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("lowermantl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("lowermantl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}


TEST_CASE("Utilities:: Thermal Conductivity Xu 2004")
{
  aspect::MaterialModel::ThermalConductivity::xu_2004<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1700, 1800, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [Pa]
  std::vector<double> pressures = {1e5, 1e9, 5e9, 10e9, 100e9};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int olivinedry_xu004_ID = 0;
  std::vector<double> olivinedry_expt_xu004_totTcon(temperatures.size());
  constexpr int wadsleydry_xu004_ID = 1;
  std::vector<double> wadsleydry_expt_xu004_totTcon(temperatures.size());
  constexpr int ringwoodry_xu004_ID = 2;
  std::vector<double> ringwoodry_expt_xu004_totTcon(temperatures.size());

  unsigned int xu004_index = ringwoodry_xu004_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_xu004_latTcond(xu004_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_xu004_totTcond(xu004_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Dry Olivine: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_xu004_latTcon = {4.11726, 1.83987, 2.00632, 2.21873, 5.46835};
  expt_xu004_latTcond[olivinedry_xu004_ID] = olivinedry_expt_xu004_latTcon;
  // Dry Wadsleyite: expected lattice thermal conductivities (k) in [W/m/K]
  std::vector<double> wadsleydry_expt_xu004_latTcon = {8.07501, 3.57699, 3.78227, 4.05482, 8.42667};
  expt_xu004_latTcond[wadsleydry_xu004_ID] = wadsleydry_expt_xu004_latTcon;
  // Dry Ringwoodite: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> ringwoodry_expt_xu004_latTcon = {9.51056, 4.20878, 4.43470, 4.73685, 9.62399};
  expt_xu004_latTcond[ringwoodry_xu004_ID] = ringwoodry_expt_xu004_latTcon;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    olivinedry_expt_xu004_totTcon[row] = olivinedry_expt_xu004_latTcon[row];
    expt_xu004_totTcond[olivinedry_xu004_ID] = olivinedry_expt_xu004_totTcon;
    wadsleydry_expt_xu004_totTcon[row] = wadsleydry_expt_xu004_latTcon[row];
    expt_xu004_totTcond[wadsleydry_xu004_ID] = wadsleydry_expt_xu004_totTcon;
    ringwoodry_expt_xu004_totTcon[row] = ringwoodry_expt_xu004_latTcon[row];
    expt_xu004_totTcond[ringwoodry_xu004_ID] = ringwoodry_expt_xu004_totTcon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < xu004_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_xu004_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_xu004_Tcond[row][col] = std::pow(expt_xu004_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_xu004_Tcond;

   INFO("Checking xu_2004 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case olivinedry_xu004_ID: // olivinedry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case wadsleydry_xu004_ID: // wadsleydry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("wadsleydry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("wadsleydry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ringwoodry_xu004_ID: // ringwoodry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ringwoodry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ringwoodry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}

TEST_CASE("Utilities:: nondimensional thermal conductivity")
{
  aspect::MaterialModel::ThermalConductivity::nondimensional_Tcond<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed

  // Assigning an array of values to in.temperature (T) in [K]
  std::vector<double> temperatures = {300, 1600, 1700, 1800, 3000};
  in.temperature = temperatures;

  // Assigning an array of values to in.pressure (P) in [GPa]
  std::vector<double> pressures = {1e5, 1e9, 5e9, 10e9, 100e9};
  in.pressure = pressures;

  // Assigning a matrix of volume fractions to in.composition (X) in [%]
  std::vector<std::vector<double>> compositions = 
  {
    {1.00, 1.00, 1.00, 1.00, 1.00},
    {0.75, 0.75, 0.75, 0.75, 0.75},
    {0.50, 0.50, 0.50, 0.50, 0.50},
    {0.25, 0.25, 0.25, 0.25, 0.25}
  };

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int nondimensional_olivine_ID = 0;
  std::vector<double> olivine_expt_nondimensional_totTcond(temperatures.size());

  unsigned int nondimensional_index = nondimensional_olivine_ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_nondimensional_latTcond(nondimensional_index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_nondimensional_radTcond(nondimensional_index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_nondimensional_totTcond(nondimensional_index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity
  
  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivine_expt_nondimensional_latTcond = {3.58888, 1.55432, 1.51642, 1.50279, 2.99729}; 
  std::vector<double> olivine_expt_nondimensional_radTcond = {0.00138, 2.23156, 2.34983, 2.45491, 3.12276};
  expt_nondimensional_latTcond[nondimensional_olivine_ID] = olivine_expt_nondimensional_latTcond;
  expt_nondimensional_radTcond[nondimensional_olivine_ID] = olivine_expt_nondimensional_radTcond;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    olivine_expt_nondimensional_totTcond[row] = olivine_expt_nondimensional_latTcond[row]+olivine_expt_nondimensional_radTcond[row];
    expt_nondimensional_totTcond[nondimensional_olivine_ID] = olivine_expt_nondimensional_totTcond;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < nondimensional_index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_nondimensional_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_nondimensional_Tcond[row][col] = std::pow(expt_nondimensional_totTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_nondimensional_Tcond;

   INFO("Checking xu_2004 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

   // Loop over the different compositions
   for (size_t row = 0; row < expected_conductivities.size(); ++row)
   {
     in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
     model.evaluate(in, out);                // Call the function to compute the thermal conductivities

     // Loop over the different combinations of pressures (P) and temperatures (T)
     for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
     {
       switch (mID) // Compare the computed thermal conductivity with the expected value
       {
         case nondimensional_olivine_ID: // olivinedry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[/] ; P= " << in.pressure[i] << "[/] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry (Nondimensional) expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry (Nondimensional) computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}

TEST_CASE("Utilities::AsciiDataLookup")
{
  using namespace dealii;

  //TODO: add support for setting data directly instead of relying on a file to load:
  aspect::Utilities::StructuredDataLookup<1> lookup(2 /*n_components*/, 1.0 /*scaling*/);
  const std::string data_filename = aspect::Utilities::expand_ASPECT_SOURCE_DIR("$ASPECT_SOURCE_DIR/data/boundary-velocity/ascii-data/test/box_2d_left.0.txt");
  lookup.load_file(data_filename, MPI_COMM_WORLD);

  INFO(lookup.get_data(Point<1>(330000./2.0),0));
  INFO(lookup.get_data(Point<1>(330000./2.0),1));
  INFO(lookup.get_gradients(Point<1>(330000./2.0),0));
  INFO(lookup.get_gradients(Point<1>(330000./2.0),1));

  REQUIRE(lookup.get_data(Point<1>(330000./2.0),0) == Approx(0.5));
  REQUIRE(lookup.get_data(Point<1>(330000./2.0),1) == Approx(0.0));
  REQUIRE(lookup.get_gradients(Point<1>(330000./2.0),0)[0] == Approx(-1.0/330000.));
  REQUIRE(lookup.get_gradients(Point<1>(330000./2.0),1)[0] == Approx(0.0));
}


TEST_CASE("Utilities::AsciiDataLookup manual dim=1")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<1> lookup(2 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"a", "b"};
  Table<1,double> table(2);
  std::vector<Table<1,double>> raw_data(2, table);

  std::vector<std::vector<double>> coordinate_values(1, std::vector<double>({1.0, 2.0}));
  // c1:
  raw_data[0](0) = 0.0;
  raw_data[0](1) = 1.0;
  // c2:
  raw_data[1](0) = 5.0;
  raw_data[1](1) = 3.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  INFO(lookup.get_data(Point<1>(1.5), 0));
  INFO(lookup.get_data(Point<1>(1.5), 1));

  REQUIRE(lookup.get_data(Point<1>(1.5),0) == Approx(0.5));
  REQUIRE(lookup.get_data(Point<1>(1.5),1) == Approx(4.0));

  REQUIRE(lookup.get_gradients(Point<1>(1.5),0)[0] == Approx(1.0));
  REQUIRE(lookup.get_gradients(Point<1>(1.5),1)[0] == Approx(-2.0));
}

TEST_CASE("Utilities::AsciiDataLookup manual dim=2")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<2> lookup(1 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"topography"};
  std::vector<Table<2,double>> raw_data(1, Table<2,double>(3,3));
  std::vector<std::vector<double>> coordinate_values(2, std::vector<double>(3, 0.));

  // x:
  coordinate_values[0] = {0., 1., 3.};
  // y:
  coordinate_values[1] = {5., 6., 7.};
  // c1:
  raw_data[0](0,0) = 1.0;
  raw_data[0](1,0) = 2.0;
  raw_data[0](2,0) = 3.0;
  raw_data[0](0,1) = 4.0;
  raw_data[0](1,1) = 5.0;
  raw_data[0](2,1) = 6.0;
  raw_data[0](0,2) = 0.0;
  raw_data[0](1,2) = 0.0;
  raw_data[0](2,2) = 0.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  REQUIRE(lookup.get_data(Point<2>(1.0,6.0),0) == Approx(5.0));
  REQUIRE(lookup.get_data(Point<2>(2.0,6.0),0) == Approx(5.5));
}

TEST_CASE("Utilities::AsciiDataLookup manual dim=2 equid")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<2> lookup(1 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"topography"};
  std::vector<Table<2,double>> raw_data(1, Table<2,double>(3,3));
  std::vector<std::vector<double>> coordinate_values(2, std::vector<double>(3, 0.));

  // x:
  coordinate_values[0] = {0., 1., 2.};
  // y:
  coordinate_values[1] = {5., 6., 7.};
  // c1:
  raw_data[0](0,0) = 1.0;
  raw_data[0](1,0) = 2.0;
  raw_data[0](2,0) = 3.0;
  raw_data[0](0,1) = 4.0;
  raw_data[0](1,1) = 5.0;
  raw_data[0](2,1) = 6.0;
  raw_data[0](0,2) = 0.0;
  raw_data[0](1,2) = 0.0;
  raw_data[0](2,2) = 0.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  REQUIRE(lookup.get_data(Point<2>(1.0,6.0),0) == Approx(5.0));
  REQUIRE(lookup.get_data(Point<2>(1.5,6.0),0) == Approx(5.5));
}

TEST_CASE("Random draw volume weighted average rotation matrix")
{
  std::vector<double> unsorted_volume_fractions = {2.,5.,1.,3.,6.,4.};
  std::vector<double> sorted_volume_fractions_ref = {1.,2.,3.,4.,5.,6.};
  const std::vector<std::size_t> permutation = aspect::Utilities::compute_sorting_permutation<double>(unsorted_volume_fractions);
  const std::vector<double> sorted_volume_fractions = aspect::Utilities::apply_permutation<double>(unsorted_volume_fractions,permutation);

  for (unsigned int i = 0; i < sorted_volume_fractions.size(); i++)
    {
      REQUIRE(sorted_volume_fractions[i] == Approx(sorted_volume_fractions_ref[i]));
    }

  const std::vector<dealii::Tensor<2,3>> unsorted_rotation_matrices =
  {
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}})
  };

  const std::vector<dealii::Tensor<2,3>> sorted_rotation_matrices_ref =
  {
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}})
  };
  const std::vector<dealii::Tensor<2,3>> sorted_rotation_matrices = aspect::Utilities::apply_permutation<dealii::Tensor<2,3>>(unsorted_rotation_matrices,permutation);
  for (unsigned int i = 0; i < sorted_rotation_matrices.size(); i++)
    {
      REQUIRE(sorted_rotation_matrices[i][0][0] == Approx(sorted_rotation_matrices_ref[i][0][0]));
    }

  std::mt19937 random_number_generator;
  random_number_generator.seed(5);
  const std::vector<dealii::Tensor<2,3>> result = aspect::Utilities::rotation_matrices_random_draw_volume_weighting(unsorted_volume_fractions,
                                                   unsorted_rotation_matrices,
                                                   25,
                                                   random_number_generator);

  const std::vector<dealii::Tensor<2,3>> result_ref =
  {
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
  };
  for (unsigned int i = 0; i < result.size(); i++)
    {
      REQUIRE(result[i][0][0] == Approx(result_ref[i][0][0]));
    }
}

TEST_CASE("wrap angle")
{
  REQUIRE(aspect::Utilities::wrap_angle(-720.) == 0.);
  REQUIRE(aspect::Utilities::wrap_angle(-540.) == 180.);
  REQUIRE(aspect::Utilities::wrap_angle(-361.) == 359.);
  REQUIRE(aspect::Utilities::wrap_angle(-360.) == 0.);
  REQUIRE(aspect::Utilities::wrap_angle(-359.) == 1.);
  REQUIRE(aspect::Utilities::wrap_angle(-270.) == 90.);
  REQUIRE(aspect::Utilities::wrap_angle(-180.) == 180.);
  REQUIRE(aspect::Utilities::wrap_angle(-90.) == 270.);
  REQUIRE(aspect::Utilities::wrap_angle(0.) == 0.);
  REQUIRE(aspect::Utilities::wrap_angle(90.) == 90.);
  REQUIRE(aspect::Utilities::wrap_angle(180.) == 180.);
  REQUIRE(aspect::Utilities::wrap_angle(270.) == 270.);
  REQUIRE(aspect::Utilities::wrap_angle(359.) == 359.);
  REQUIRE(aspect::Utilities::wrap_angle(360.) == 0.);
  REQUIRE(aspect::Utilities::wrap_angle(361.) == 1.);
  REQUIRE(aspect::Utilities::wrap_angle(540.) == 180.);
  REQUIRE(aspect::Utilities::wrap_angle(720.) == 0.);
}

TEST_CASE("CPO elastic tensor transform functions")
{
  dealii::SymmetricTensor<2,6> reference_elastic_tensor({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21});

// first test whether the functions are invertable
  {
    dealii::SymmetricTensor<2,6> result_up_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_up_down[i][j]));
          }
      }
  }
  {
    dealii::SymmetricTensor<2,6> result_down_up = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(aspect::Utilities::Tensors::to_voigt_stiffness_vector(reference_elastic_tensor));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_down_up[i][j]));
          }
      }
  }
  {
    dealii::SymmetricTensor<2,6> result_up_2down_up = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(aspect::Utilities::Tensors::to_voigt_stiffness_vector(aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));

    for (size_t i = 0; i < 6; i++)
      {
        for (size_t j = 0; j < 6; j++)
          {
            REQUIRE(reference_elastic_tensor[i][j] == Approx(result_up_2down_up[i][j]));
          }
      }
  }

  // test rotations
  // rotation matrix
  dealii::Tensor<2,3> rotation_tensor;

  {

    // fill the rotation matrix with a rotations in all directions
    {
      double radians = 0;
      double alpha = radians;
      double beta = radians;
      double gamma = radians;
      rotation_tensor[0][0] = std::cos(alpha) * std::cos(beta);
      rotation_tensor[0][1] = std::sin(alpha) * std::cos(beta);
      rotation_tensor[0][2] = -std::sin(beta);
      rotation_tensor[1][0] = std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha)*std::cos(gamma);
      rotation_tensor[1][1] = std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha)*std::cos(gamma);
      rotation_tensor[1][2] = std::cos(beta) * std::sin(gamma);
      rotation_tensor[2][0] = std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha)*std::sin(gamma);
      rotation_tensor[2][1] = std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha)*std::sin(gamma);
      rotation_tensor[2][2] = std::cos(beta) * std::cos(gamma);
    }

    {
      dealii::SymmetricTensor<4,3> result_full_stiffness_tensor = aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor);
      dealii::SymmetricTensor<4,3> result_full_stiffness_tensor_rotate_zero = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,result_full_stiffness_tensor);

      // first check that one the tensors didn't change with zero rotation
      for (size_t i = 0; i < 3; i++)
        {
          for (size_t j = 0; j < 3; j++)
            {
              for (size_t k = 0; k < 3; k++)
                {
                  for (size_t l = 0; l < 3; l++)
                    {
                      REQUIRE(result_full_stiffness_tensor[i][j][k][l] == Approx(result_full_stiffness_tensor_rotate_zero[i][j][k][l]));
                    }
                }
            }
        }
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(
                                                               aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,
                                                                   aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor,reference_elastic_tensor);

      // first check that one the tensors didn't change with zero rotation
      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_up_1_rotate_down[i][j] == Approx(reference_elastic_tensor[i][j]));
              REQUIRE(result_1_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }
    }

    // fill the rotation matrix with a rotations in all directions
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = radians;
      double beta = radians;
      double gamma = radians;
      rotation_tensor[0][0] = std::cos(alpha) * std::cos(beta);
      rotation_tensor[0][1] = std::sin(alpha) * std::cos(beta);
      rotation_tensor[0][2] = -std::sin(beta);
      rotation_tensor[1][0] = std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha)*std::cos(gamma);
      rotation_tensor[1][1] = std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha)*std::cos(gamma);
      rotation_tensor[1][2] = std::cos(beta) * std::sin(gamma);
      rotation_tensor[2][0] = std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha)*std::sin(gamma);
      rotation_tensor[2][1] = std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha)*std::sin(gamma);
      rotation_tensor[2][2] = std::cos(beta) * std::cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(
                                                               aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,
                                                                   aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor,reference_elastic_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::SymmetricTensor<4,3> result_up_10_rotate = aspect::Utilities::Tensors::to_full_stiffness_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,result_up_10_rotate);
          result_5_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, result_5_rotate);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              // This test doesn't work when rotating in all rotations at the same time.
              //REQUIRE(result_1_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }
  {
    // fill the rotation matrix with a rotations in the alpha direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = radians;
      double beta = 0;
      double gamma = 0;
      rotation_tensor[0][0] = std::cos(alpha) * std::cos(beta);
      rotation_tensor[0][1] = std::sin(alpha) * std::cos(beta);
      rotation_tensor[0][2] = -std::sin(beta);
      rotation_tensor[1][0] = std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha)*std::cos(gamma);
      rotation_tensor[1][1] = std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha)*std::cos(gamma);
      rotation_tensor[1][2] = std::cos(beta) * std::sin(gamma);
      rotation_tensor[2][0] = std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha)*std::sin(gamma);
      rotation_tensor[2][1] = std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha)*std::sin(gamma);
      rotation_tensor[2][2] = std::cos(beta) * std::cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(
                                                               aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,
                                                                   aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor,reference_elastic_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::SymmetricTensor<4,3> result_up_10_rotate = aspect::Utilities::Tensors::to_full_stiffness_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor, result_up_10_rotate);
          result_5_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, result_5_rotate);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }
  {
    // fill the rotation matrix with a rotations in the beta direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = 0;
      double beta = radians;
      double gamma = 0;
      rotation_tensor[0][0] = std::cos(alpha) * std::cos(beta);
      rotation_tensor[0][1] = std::sin(alpha) * std::cos(beta);
      rotation_tensor[0][2] = -std::sin(beta);
      rotation_tensor[1][0] = std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha)*std::cos(gamma);
      rotation_tensor[1][1] = std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha)*std::cos(gamma);
      rotation_tensor[1][2] = std::cos(beta) * std::sin(gamma);
      rotation_tensor[2][0] = std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha)*std::sin(gamma);
      rotation_tensor[2][1] = std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha)*std::sin(gamma);
      rotation_tensor[2][2] = std::cos(beta) * std::cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(
                                                               aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,
                                                                   aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, reference_elastic_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::SymmetricTensor<4,3> result_up_10_rotate = aspect::Utilities::Tensors::to_full_stiffness_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor, result_up_10_rotate);
          result_5_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, result_5_rotate);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }

  {
    // fill the rotation matrix with a rotations in the gamma direction
    {
      double radians = (dealii::numbers::PI/180.0)*(360/5); //0.35*dealii::numbers::PI; //(dealii::numbers::PI/180.0)*36;
      double alpha = 0;
      double beta = 0;
      double gamma = radians;
      rotation_tensor[0][0] = std::cos(alpha) * std::cos(beta);
      rotation_tensor[0][1] = std::sin(alpha) * std::cos(beta);
      rotation_tensor[0][2] = -std::sin(beta);
      rotation_tensor[1][0] = std::cos(alpha) * std::sin(beta) * std::sin(gamma) - std::sin(alpha)*std::cos(gamma);
      rotation_tensor[1][1] = std::sin(alpha) * std::sin(beta) * std::sin(gamma) + std::cos(alpha)*std::cos(gamma);
      rotation_tensor[1][2] = std::cos(beta) * std::sin(gamma);
      rotation_tensor[2][0] = std::cos(alpha) * std::sin(beta) * std::cos(gamma) + std::sin(alpha)*std::sin(gamma);
      rotation_tensor[2][1] = std::sin(alpha) * std::sin(beta) * std::cos(gamma) - std::cos(alpha)*std::sin(gamma);
      rotation_tensor[2][2] = std::cos(beta) * std::cos(gamma);
    }

    {
      dealii::SymmetricTensor<2,6> result_up_1_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(
                                                               aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor,
                                                                   aspect::Utilities::Tensors::to_full_stiffness_tensor(reference_elastic_tensor)));
      dealii::SymmetricTensor<2,6> result_1_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, reference_elastic_tensor);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_1_rotate[i][j] == Approx(result_up_1_rotate_down[i][j]));
            }
        }

      dealii::SymmetricTensor<4,3> result_up_10_rotate = aspect::Utilities::Tensors::to_full_stiffness_tensor(result_up_1_rotate_down);

      dealii::SymmetricTensor<2,6> result_5_rotate = result_1_rotate;

      for (size_t i = 0; i < 4; i++)
        {
          result_up_10_rotate = aspect::Utilities::Tensors::rotate_full_stiffness_tensor(rotation_tensor, result_up_10_rotate);
          result_5_rotate = aspect::Utilities::Tensors::rotate_voigt_stiffness_matrix(rotation_tensor, result_5_rotate);
        }

      dealii::SymmetricTensor<2,6> result_up_10_rotate_down = aspect::Utilities::Tensors::to_voigt_stiffness_matrix(result_up_10_rotate);

      for (size_t i = 0; i < 6; i++)
        {
          for (size_t j = 0; j < 6; j++)
            {
              REQUIRE(result_5_rotate[i][j] == Approx(result_up_10_rotate_down[i][j]));
              REQUIRE(result_5_rotate[i][j] == Approx(reference_elastic_tensor[i][j]));
            }
        }
    }
  }

  /**
   * test Levi-Cevita tensor function
   */
  {

    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][0][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][0][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][0][2] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][1][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][1][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][1][2] == Approx(1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][2][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][2][1] == Approx(-1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[0][2][2] == Approx(0.0));

    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][0][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][0][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][0][2] == Approx(-1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][1][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][1][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][1][2] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][2][0] == Approx(1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][2][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[1][2][2] == Approx(0.0));

    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][0][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][0][1] == Approx(1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][0][2] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][1][0] == Approx(-1.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][1][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][1][2] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][2][0] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][2][1] == Approx(0.0));
    REQUIRE(aspect::Utilities::Tensors::levi_civita<3>()[2][2][2] == Approx(0.0));
  }
}

TEST_CASE("Utilities::string_to_unsigned_int")
{
  CHECK(aspect::Utilities::string_to_unsigned_int("1234") == 1234);

  CHECK(aspect::Utilities::string_to_unsigned_int(std::vector<std::string>({"234","0","1"}))
        == std::vector<unsigned int>({234,0,1}));

  CHECK(aspect::Utilities::string_to_unsigned_int(std::vector<std::string>({"42"}))
        == std::vector<unsigned int>({42}));

  CHECK(aspect::Utilities::string_to_unsigned_int(std::vector<std::string>({}))
        == std::vector<unsigned int>());
}

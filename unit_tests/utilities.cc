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
  std::vector<double> pressures = {1e5, 1e9, 5e9, 10e9, 100e9};
  in.pressure = pressures;

  // Assigning a lithology index to in.composition 
  //(0: pyrolite, 1: harzburgite, 2: meta-MORB, 3: dunite, 99: test)
  std::vector<std::vector<double>> lithologies = 
  {
    {0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1},
    {2, 2, 2, 2, 2},
    {3, 3, 3, 3, 3}
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

  unsigned int mineralpar_index = akimotoite_exptID+1; // Number of minerals
  unsigned int rockspar_index = duniteOl_exptID+1; // Number of rocks

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> expt_minerals_latTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_minerals_radTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_minerals_totTcond(temperatures.size(), std::vector<double>(mineralpar_index, 0.0)); // Total thermal conductivity

  // Preallocate matrixes for storing thermal conductivities of rocks
  std::vector<std::vector<double>> expt_rocks_latTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> expt_rocks_radTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> expt_rocks_totTcond(temperatures.size(), std::vector<double>(rockspar_index, 0.0)); // Total thermal conductivity

  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> olivinedry_expt_latTcon = {3.58888, 1.58719, 2.36282, 3.67868, 4.25767};
  std::vector<double> olivinedry_expt_radTcon = {0.00138288, 2.23152, 2.34978, 2.45486, 3.12273};
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> wadsleydry_expt_latTcon = {5.88364, 3.08196, 3.21488, 3.22817, 2.77368};
  std::vector<double> wadsleydry_expt_radTcon = {1.1834e-9, 1.51515, 1.71067, 1.87995, 2.71457};
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> ringwoodry_expt_latTcon = {4.98456, 2.17063, 2.41846, 3.19795, 5.89102};
  std::vector<double> ringwoodry_expt_radTcon = {5.52667e-10, 0.74367, 0.85004, 0.94227, 1.38818};
  // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigm100Mg_expt_latTcon = {10.69504, 2.03810, 2.21218, 2.47430, 5.69003};
  std::vector<double> brigm100Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigma97Mg_expt_latTcon = {5.73751, 2.29118, 2.48780, 2.75655, 6.95242};
  std::vector<double> brigma97Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigma90Mg_expt_latTcon = {3.79122, 2.91572, 3.14193, 3.43191, 8.92076};
  std::vector<double> brigma90Mg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> brigmaAlMg_expt_latTcon = {6.30403, 2.31503, 2.56860, 2.91357, 7.88626};
  std::vector<double> brigmaAlMg_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> brigmaFeAl_expt_latTcon = {3.99962, 1.87753, 2.05369, 2.30265, 6.17198};
  std::vector<double> brigmaFeAl_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> opxenstati_expt_latTcon = {5.79950, 2.55142, 3.15043, 3.26008, 2.56641};
  std::vector<double> opxenstati_expt_radTcon = {5.45056e-5, 2.93307, 3.08161, 3.20925, 3.90772};
  // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> cpxdiopsid_expt_latTcon = {5.99273, 2.65322, 3.13819, 3.61301, 3.41681};
  std::vector<double> cpxdiopsid_expt_radTcon = {2.64987e-5, 2.99817, 3.14622, 3.27243, 3.94085};
  // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> grtpyropes_expt_latTcon = {4.38827, 2.15733, 2.64684, 3.54118, 4.22405};
  std::vector<double> grtpyropes_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> grtgrossul_expt_latTcon = {4.08838, 1.91337, 2.26557, 3.05964, 4.01413};
  std::vector<double> grtgrossul_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> grtalmandi_expt_latTcon = {3.39124, 1.70813, 2.17978, 3.13522, 4.07543};
  std::vector<double> grtalmandi_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
  std::vector<double> grtmajorit_expt_latTcon = {9.73983, 4.24079, 4.57076, 5.13035, 4.76249};
  std::vector<double> grtmajorit_expt_radTcon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> quartzpure_expt_latTcon = {9.53243, 1.84339, 2.49820, 2.47676, 1.49323};
  std::vector<double> quartzpure_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> coesitSiO2_expt_latTcon = {7.21196, 1.31776, 1.23920, 1.17013, 0.84983};
  std::vector<double> coesitSiO2_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> stishovite_expt_latTcon = {67.68617, 29.30897, 28.43361, 27.62648, 29.60473};
  std::vector<double> stishovite_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> stisho05Al_expt_latTcon = {24.18571, 10.48883, 10.35956, 10.44473, 15.11569};
  std::vector<double> stisho05Al_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> antigor010_expt_latTcon = {4.55589, 1.99618, 2.41198, 3.15852, 3.57446};
  std::vector<double> antigor010_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> antigor001_expt_latTcon = {1.06669, 0.49210, 1.01828, 1.51143, 1.48578};
  std::vector<double> antigor001_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> phaseDFeAl_expt_latTcon = {2.59325, 1.13915, 1.31955, 1.59985, 6.12649};
  std::vector<double> phaseDFeAl_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> phaseD02Al_expt_latTcon = {3.60666, 1.56888, 1.65207, 1.95489, 8.61776};
  std::vector<double> phaseD02Al_expt_radTcon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper08_expt_latTcon = {5.08425, 2.20657, 2.24939, 2.50821, 14.39941};
  std::vector<double> ferroper08_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper10_expt_latTcon = {4.48610, 1.94695, 1.98084, 2.19296, 12.65637};
  std::vector<double> ferroper10_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg80Fe20O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper20_expt_latTcon = {3.48647, 3.39116, 3.56449, 3.77742, 7.57324};
  std::vector<double> ferroper20_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> ferroper56_expt_latTcon = {2.69167, 1.23171, 1.55079, 2.02398, 7.04069};
  std::vector<double> ferroper56_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> davemaoite_expt_latTcon = {10.86634, 4.77342, 5.15289, 5.86997, 13.48443};
  std::vector<double> davemaoite_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> newhexAlph_expt_latTcon = {10.59581, 4.58812, 4.45113, 4.32578, 12.15184};
  std::vector<double> newhexAlph_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> akimotoite_expt_latTcon = {10.69504, 2.03810, 2.21218, 2.47430, 5.69003};
  std::vector<double> akimotoite_expt_radTcon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};

  // pyrolite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> pyrolite_expt_latTcon = {1, 1, 1, 1, 1};
  std::vector<double> pyrolite_expt_radTcon = {1, 1, 1, 1, 1};
  // harzburg: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> harzburg_expt_latTcon = {1, 1, 1, 1, 1};
  std::vector<double> harzburg_expt_radTcon = {1, 1, 1, 1, 1};
  // metaMORB: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> metaMORB_expt_latTcon = {1, 1, 1, 1, 1};
  std::vector<double> metaMORB_expt_radTcon = {1, 1, 1, 1, 1};
  // duniteOl: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> duniteOl_expt_latTcon = {1, 1, 1, 1, 1};
  std::vector<double> duniteOl_expt_radTcon = {1, 1, 1, 1, 1};

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    olivinedry_expt_totTcon[row] = olivinedry_expt_latTcon[row] + olivinedry_expt_radTcon[row];
    expt_minerals_latTcond[row][olivinedry_exptID] = olivinedry_expt_latTcon[row];
    expt_minerals_radTcond[row][olivinedry_exptID] = olivinedry_expt_radTcon[row];  
    expt_minerals_totTcond[row][olivinedry_exptID] = olivinedry_expt_totTcon[row];
    // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    wadsleydry_expt_totTcon[row] = wadsleydry_expt_latTcon[row] + wadsleydry_expt_radTcon[row];
    expt_minerals_latTcond[row][wadsleydry_exptID] = wadsleydry_expt_latTcon[row];
    expt_minerals_radTcond[row][wadsleydry_exptID] = wadsleydry_expt_radTcon[row];   
    expt_minerals_totTcond[row][wadsleydry_exptID] = wadsleydry_expt_totTcon[row];
    // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    ringwoodry_expt_totTcon[row] = ringwoodry_expt_latTcon[row] + ringwoodry_expt_radTcon[row];
    expt_minerals_latTcond[row][ringwoodry_exptID] = ringwoodry_expt_latTcon[row];
    expt_minerals_radTcond[row][ringwoodry_exptID] = ringwoodry_expt_radTcon[row]; 
    expt_minerals_totTcond[row][ringwoodry_exptID] = ringwoodry_expt_totTcon[row];
    // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigm100Mg_expt_totTcon[row] = brigm100Mg_expt_latTcon[row] + brigm100Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_radTcon[row]; 
    expt_minerals_totTcond[row][brigm100Mg_exptID] = brigm100Mg_expt_totTcon[row];
    // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigma97Mg_expt_totTcon[row] = brigma97Mg_expt_latTcon[row] + brigma97Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_radTcon[row];
    expt_minerals_totTcond[row][brigma97Mg_exptID] = brigma97Mg_expt_totTcon[row];
    // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigma90Mg_expt_totTcon[row] = brigma90Mg_expt_latTcon[row] + brigma90Mg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_radTcon[row];
    expt_minerals_totTcond[row][brigma90Mg_exptID] = brigma90Mg_expt_totTcon[row];
    // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    brigmaAlMg_expt_totTcon[row] = brigmaAlMg_expt_latTcon[row] + brigmaAlMg_expt_radTcon[row];
    expt_minerals_latTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_latTcon[row];
    expt_minerals_radTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_radTcon[row]; 
    expt_minerals_totTcond[row][brigmaAlMg_exptID] = brigmaAlMg_expt_totTcon[row];
    // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    brigmaFeAl_expt_totTcon[row] = brigmaFeAl_expt_latTcon[row] + brigmaFeAl_expt_radTcon[row];
    expt_minerals_latTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_latTcon[row];
    expt_minerals_radTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_radTcon[row];
    expt_minerals_totTcond[row][brigmaFeAl_exptID] = brigmaFeAl_expt_totTcon[row];
    // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    opxenstati_expt_totTcon[row] = opxenstati_expt_latTcon[row] + opxenstati_expt_radTcon[row];
    expt_minerals_latTcond[row][opxenstati_exptID] = opxenstati_expt_latTcon[row];
    expt_minerals_radTcond[row][opxenstati_exptID] = opxenstati_expt_radTcon[row];
    expt_minerals_totTcond[row][opxenstati_exptID] = opxenstati_expt_totTcon[row];
    // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    cpxdiopsid_expt_totTcon[row] = cpxdiopsid_expt_latTcon[row] + cpxdiopsid_expt_radTcon[row];
    expt_minerals_latTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_latTcon[row];
    expt_minerals_radTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_radTcon[row];
    expt_minerals_totTcond[row][cpxdiopsid_exptID] = cpxdiopsid_expt_totTcon[row];
    // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    grtpyropes_expt_totTcon[row] = grtpyropes_expt_latTcon[row] + grtpyropes_expt_radTcon[row];
    expt_minerals_latTcond[row][grtpyropes_exptID] = grtpyropes_expt_latTcon[row];
    expt_minerals_radTcond[row][grtpyropes_exptID] = grtpyropes_expt_radTcon[row];
    expt_minerals_totTcond[row][grtpyropes_exptID] = grtpyropes_expt_totTcon[row];
    // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    grtgrossul_expt_totTcon[row] = grtgrossul_expt_latTcon[row] + grtgrossul_expt_radTcon[row];
    expt_minerals_latTcond[row][grtgrossul_exptID] = grtgrossul_expt_latTcon[row];
    expt_minerals_radTcond[row][grtgrossul_exptID] = grtgrossul_expt_radTcon[row];
    expt_minerals_totTcond[row][grtgrossul_exptID] = grtgrossul_expt_totTcon[row];
    // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    grtalmandi_expt_totTcon[row] = grtalmandi_expt_latTcon[row] + grtalmandi_expt_radTcon[row];
    expt_minerals_latTcond[row][grtalmandi_exptID] = grtalmandi_expt_latTcon[row];
    expt_minerals_radTcond[row][grtalmandi_exptID] = grtalmandi_expt_radTcon[row];
    expt_minerals_totTcond[row][grtalmandi_exptID] = grtalmandi_expt_totTcon[row];
    // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
    grtmajorit_expt_totTcon[row] = grtmajorit_expt_latTcon[row] + grtmajorit_expt_radTcon[row];
    expt_minerals_latTcond[row][grtmajorit_exptID] = grtmajorit_expt_latTcon[row];
    expt_minerals_radTcond[row][grtmajorit_exptID] = grtmajorit_expt_radTcon[row];   
    expt_minerals_totTcond[row][grtmajorit_exptID] = grtmajorit_expt_totTcon[row];
    // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    quartzpure_expt_totTcon[row] = quartzpure_expt_latTcon[row] + quartzpure_expt_radTcon[row];
    expt_minerals_latTcond[row][quartzpure_exptID] = quartzpure_expt_latTcon[row];
    expt_minerals_radTcond[row][quartzpure_exptID] = quartzpure_expt_radTcon[row];
    expt_minerals_totTcond[row][quartzpure_exptID] = quartzpure_expt_totTcon[row];
    // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    coesitSiO2_expt_totTcon[row] = coesitSiO2_expt_latTcon[row] + coesitSiO2_expt_radTcon[row];
    expt_minerals_latTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_latTcon[row];
    expt_minerals_radTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_radTcon[row];
    expt_minerals_totTcond[row][coesitSiO2_exptID] = coesitSiO2_expt_totTcon[row];
    // stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    stishovite_expt_totTcon[row] = stishovite_expt_latTcon[row] + stishovite_expt_radTcon[row];
    expt_minerals_latTcond[row][stishovite_exptID] = stishovite_expt_latTcon[row];
    expt_minerals_radTcond[row][stishovite_exptID] = stishovite_expt_radTcon[row];    
    expt_minerals_totTcond[row][stishovite_exptID] = stishovite_expt_totTcon[row];
    // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    stisho05Al_expt_totTcon[row] = stisho05Al_expt_latTcon[row] + stisho05Al_expt_radTcon[row];
    expt_minerals_latTcond[row][stisho05Al_exptID] = stisho05Al_expt_latTcon[row];
    expt_minerals_radTcond[row][stisho05Al_exptID] = stisho05Al_expt_radTcon[row];
    expt_minerals_totTcond[row][stisho05Al_exptID] = stisho05Al_expt_totTcon[row];
    // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    antigor010_expt_totTcon[row] = antigor010_expt_latTcon[row] + antigor010_expt_radTcon[row];
    expt_minerals_latTcond[row][antigor010_exptID] = antigor010_expt_latTcon[row];
    expt_minerals_radTcond[row][antigor010_exptID] = antigor010_expt_radTcon[row];
    expt_minerals_totTcond[row][antigor010_exptID] = antigor010_expt_totTcon[row];
    // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    antigor001_expt_totTcon[row] = antigor001_expt_latTcon[row] + antigor001_expt_radTcon[row]; 
    expt_minerals_latTcond[row][antigor001_exptID] = antigor001_expt_latTcon[row];
    expt_minerals_radTcond[row][antigor001_exptID] = antigor001_expt_radTcon[row];
    expt_minerals_totTcond[row][antigor001_exptID] = antigor001_expt_totTcon[row];   
    // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    phaseDFeAl_expt_totTcon[row] = phaseDFeAl_expt_latTcon[row] + phaseDFeAl_expt_radTcon[row];
    expt_minerals_latTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_latTcon[row];
    expt_minerals_radTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_radTcon[row];
    expt_minerals_totTcond[row][phaseDFeAl_exptID] = phaseDFeAl_expt_totTcon[row]; 
    // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    phaseD02Al_expt_totTcon[row] = phaseD02Al_expt_latTcon[row] + phaseD02Al_expt_radTcon[row];
    expt_minerals_latTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_latTcon[row];
    expt_minerals_radTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_radTcon[row];
    expt_minerals_totTcond[row][phaseD02Al_exptID] = phaseD02Al_expt_totTcon[row];
    // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper08_expt_totTcon[row] = ferroper08_expt_latTcon[row] + ferroper08_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper08_exptID] = ferroper08_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper08_exptID] = ferroper08_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper08_exptID] = ferroper08_expt_totTcon[row];
    // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper10_expt_totTcon[row] = ferroper10_expt_latTcon[row] + ferroper10_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper10_exptID] = ferroper10_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper10_exptID] = ferroper10_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper10_exptID] = ferroper10_expt_totTcon[row];
    // Ferropericlase (Mg80Fe20O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper20_expt_totTcon[row] = ferroper20_expt_latTcon[row] + ferroper20_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper20_exptID] = ferroper20_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper20_exptID] = ferroper20_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper20_exptID] = ferroper20_expt_totTcon[row];
    // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    ferroper56_expt_totTcon[row] = ferroper56_expt_latTcon[row] + ferroper56_expt_radTcon[row];
    expt_minerals_latTcond[row][ferroper56_exptID] = ferroper56_expt_latTcon[row];
    expt_minerals_radTcond[row][ferroper56_exptID] = ferroper56_expt_radTcon[row];
    expt_minerals_totTcond[row][ferroper56_exptID] = ferroper56_expt_totTcon[row];
    // davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    davemaoite_expt_totTcon[row] = davemaoite_expt_latTcon[row] + davemaoite_expt_radTcon[row]; 
    expt_minerals_latTcond[row][davemaoite_exptID] = davemaoite_expt_latTcon[row];
    expt_minerals_radTcond[row][davemaoite_exptID] = davemaoite_expt_radTcon[row];
    expt_minerals_totTcond[row][davemaoite_exptID] = davemaoite_expt_totTcon[row];
    // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
    newhexAlph_expt_totTcon[row] = newhexAlph_expt_latTcon[row] + newhexAlph_expt_radTcon[row];
    expt_minerals_latTcond[row][newhexAlph_exptID] = newhexAlph_expt_latTcon[row];
    expt_minerals_radTcond[row][newhexAlph_exptID] = newhexAlph_expt_radTcon[row];
    expt_minerals_totTcond[row][newhexAlph_exptID] = newhexAlph_expt_totTcon[row];
    // akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    akimotoite_expt_totTcon[row] = akimotoite_expt_latTcon[row] + akimotoite_expt_radTcon[row];
    expt_minerals_latTcond[row][akimotoite_exptID] = akimotoite_expt_latTcon[row];
    expt_minerals_radTcond[row][akimotoite_exptID] = akimotoite_expt_radTcon[row]; 
    expt_minerals_totTcond[row][akimotoite_exptID] = akimotoite_expt_totTcon[row];

    // pyrolite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    pyrolite_expt_totTcon[row] = pyrolite_expt_latTcon[row] + pyrolite_expt_radTcon[row];
    expt_rocks_latTcond[row][pyrolite_exptID] = pyrolite_expt_latTcon[row];
    expt_rocks_radTcond[row][pyrolite_exptID] = pyrolite_expt_radTcon[row];
    expt_rocks_totTcond[row][pyrolite_exptID] = pyrolite_expt_totTcon[row];
    // harzburg: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    harzburg_expt_totTcon[row] = harzburg_expt_latTcon[row] + harzburg_expt_radTcon[row];
    expt_rocks_latTcond[row][harzburg_exptID] = harzburg_expt_latTcon[row];
    expt_rocks_radTcond[row][harzburg_exptID] = harzburg_expt_radTcon[row];
    expt_rocks_totTcond[row][harzburg_exptID] = harzburg_expt_totTcon[row];
    // metaMORB: expected lattice and radiative thermal conductivities (k) in [W/m/K]
    metaMORB_expt_totTcon[row] = metaMORB_expt_latTcon[row] + metaMORB_expt_radTcon[row];
    expt_rocks_latTcond[row][metaMORB_exptID] = metaMORB_expt_latTcon[row];
    expt_rocks_radTcond[row][metaMORB_exptID] = metaMORB_expt_radTcon[row];
    expt_rocks_totTcond[row][metaMORB_exptID] = metaMORB_expt_totTcon[row];
    // duniteOl: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
    duniteOl_expt_totTcon[row] = duniteOl_expt_latTcon[row] + duniteOl_expt_radTcon[row];
    expt_rocks_latTcond[row][duniteOl_exptID] = duniteOl_expt_latTcon[row];
    expt_rocks_radTcond[row][duniteOl_exptID] = duniteOl_expt_radTcon[row];
    expt_rocks_totTcond[row][duniteOl_exptID] = duniteOl_expt_totTcon[row];
  }

  // Loop over all lithologies
  for (unsigned int lID = 0; lID < lithologies.size(); ++lID)
  {
    
    INFO("Checking thermal conductivity (k) for different lithologies as a function of temperature (T) and pressure (P)");

    // Loop over the different combinations of pressures (P) and temperatures (T)
    for (size_t row = 0; row < expt_rocks_totTcond[lID].size(); ++row)
    {

      in.composition[0] = lithologies[row];  // Assign the current lithology as model in.composition input
      model.evaluate(in, out);  // Call the function to compute the thermal conductivities

      switch (lID) // Compare the computed thermal conductivity with the expected value
      {
       case pyrolite_exptID: // pyrolite
       {
         INFO("Lithology: " << in.composition[row][lID] << " (Pyrolite) ; Conditions: T = " << in.temperature[row] << "[K] ; P = " << in.pressure[row] << "[Pa]");
         INFO("pyrolite expected k = " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
         INFO("pyrolite computed k = " << out.thermal_conductivities[row] << "[W/m/K]");
         REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
         break;
       }
       case harzburg_exptID: // harzburgite
       {
         INFO("Conditions: T= " << in.temperature[row] << "[K] ; P= " << in.pressure[row] << "[Pa]");
         INFO("harzburg expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
         INFO("harzburg computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
         REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
         break;
        }
       case metaMORB_exptID: // metabasalt (MORB)
       {
         INFO("Conditions: T= " << in.temperature[row] << "[K] ; P= " << in.pressure[row] << "[Pa]");
         INFO("metaMORB expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
         INFO("metaMORB computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
         REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
         break;
        }
       case duniteOl_exptID: // dunite (100% olivine)
       {
         INFO("Conditions: T= " << in.temperature[row] << "[K] ; P= " << in.pressure[row] << "[Pa]");
         INFO("duniteOl expected k= " << expt_rocks_totTcond[row][lID] << "[W/m/K]");
         INFO("duniteOl computed k= " << out.thermal_conductivities[row] << "[W/m/K]");
         REQUIRE(out.thermal_conductivities[row] == Approx(expt_rocks_totTcond[row][lID]));
         break;
       }
      }
    }
  }

 
 /*
 // Loop over all mID values
 for (unsigned int mID = 0; mID < mineralpar_index; ++mID)
 {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_total_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
   {
     for (size_t col = 0; col < compositions[row].size(); ++col)
     {
       expected_total_Tcond[row][col] = std::pow(expt_minerals_totTcond[mID][col], compositions[row][col]);
     }
   }

   std::vector<std::vector<double>> expected_conductivities = expected_total_Tcond;

   INFO("Checking thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

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
         case olivinedry_exptID: // olivinedry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("olivinedry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("olivinedry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case wadsleydry_exptID: // wadsleydry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("wadsleydry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("wadsleydry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ringwoodry_exptID: // ringwoodry
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ringwoodry expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ringwoodry computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case brigm100Mg_exptID: // brigm100Mg
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigm100Mg expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigm100Mg computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case brigma97Mg_exptID: // brigma97Mg
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigma97Mg expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigma97Mg computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case brigma90Mg_exptID: // brigma90Mg
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigma90Mg expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigma90Mg computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case brigmaAlMg_exptID: // brigmaAlMg
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigmaAlMg expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigmaAlMg computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case brigmaFeAl_exptID: // brigmaFeAl
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("brigmaFeAl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("brigmaFeAl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case opxenstati_exptID: // opxenstati
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("opxenstati expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("opxenstati computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case cpxdiopsid_exptID: // cpxdiopsid
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("cpxdiopsid expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("cpxdiopsid computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case grtpyropes_exptID: // grtpyropes
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtpyropes expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtpyropes computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case grtgrossul_exptID: // grtgrossul
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtgrossul expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtgrossul computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case grtalmandi_exptID: // grtalmandi
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtalmandi expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtalmandi computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case grtmajorit_exptID: // grtmajorit
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("grtmajorit expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("grtmajorit computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case quartzpure_exptID: // quartzpure
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("quartzpure expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("quartzpure computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case coesitSiO2_exptID: // coesitSiO2
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("coesitSiO2 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("coesitSiO2 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case stishovite_exptID: // stishovite
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("stishovite expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("stishovite computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case stisho05Al_exptID: // stisho05Al
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("stisho05Al expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("stisho05Al computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case antigor010_exptID: // antigor010
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("antigor010 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("antigor010 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case antigor001_exptID: // antigor001
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("antigor001 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("antigor001 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case phaseDFeAl_exptID: // phaseDFeAl
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("phaseDFeAl expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("phaseDFeAl computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case phaseD02Al_exptID: // phaseD02Al
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("phaseD02Al expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("phaseD02Al computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ferroper08_exptID: // ferroper08
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ferroper08 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ferroper08 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ferroper10_exptID: // ferroper10
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ferroper10 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ferroper10 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ferroper20_exptID: // ferroper20
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ferroper20 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ferroper20 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case ferroper56_exptID: // ferroper56
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("ferroper56 expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("ferroper56 computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case davemaoite_exptID: // davemaoite
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("davemaoite expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("davemaoite computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case newhexAlph_exptID: // newhexAlph
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("newhexAlph expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("newhexAlph computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case akimotoite_exptID: // akimotoite
         {
           INFO("Conditions: T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("akimotoite expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("akimotoite computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
       } 
     }
   }
 } 
   */
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

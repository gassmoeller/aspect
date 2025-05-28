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
#include <aspect/material_model/thermal_conductivity/PT_dep_R_bounded.h>
#include <aspect/material_model/thermal_conductivity/Hofmeister1999.h>
// #include <aspect/material_model/thermal_conductivity/Anderson1987.h>
// #include <aspect/material_model/thermal_conductivity/Gerya2021.h>
// #include <aspect/material_model/thermal_conductivity/GroseAfonso2019.h>
// #include <aspect/material_model/thermal_conductivity/Hofmeister2005.h>
// #include <aspect/material_model/thermal_conductivity/HofmeisterBranlund2015.h>
// #include <aspect/material_model/thermal_conductivity/Stackhouse2015.h>
#include <aspect/material_model/thermal_conductivity/Tosi2016.h>
#include <aspect/material_model/thermal_conductivity/Xu2004.h>

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

TEST_CASE("Utilities::PT dependent thermal conductivity Enrico")
{
  aspect::MaterialModel::ThermalConductivity::PTdepRbounded<3> model;
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

  // unsigned int MineralPar_Index = 0; // Initialize the counter

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  constexpr int OlivineDry_ExptID = 0;
  std::vector<double> OlivineDry_Expt_TotTCon(temperatures.size());
  constexpr int WadsleyDry_ExptID = 1;
  std::vector<double> WadsleyDry_Expt_TotTCon(temperatures.size());
  constexpr int RingwooDry_ExptID = 2;
  std::vector<double> RingwooDry_Expt_TotTCon(temperatures.size());
  constexpr int En100Brigm_ExptID = 3;
  std::vector<double> En100Brigm_Expt_TotTCon(temperatures.size());
  constexpr int En97Brigma_ExptID = 4;
  std::vector<double> En97Brigma_Expt_TotTCon(temperatures.size());
  constexpr int En90Brigma_ExptID = 5;
  std::vector<double> En90Brigma_Expt_TotTCon(temperatures.size());
  constexpr int AlMgBrigma_ExptID = 6;
  std::vector<double> AlMgBrigma_Expt_TotTCon(temperatures.size());
  constexpr int FeAlBrigma_ExptID = 7;
  std::vector<double> FeAlBrigma_Expt_TotTCon(temperatures.size());
  constexpr int OpxEnstati_ExptID = 8;
  std::vector<double> OpxEnstati_Expt_TotTCon(temperatures.size());
  constexpr int CpxDiopsid_ExptID = 9;
  std::vector<double> CpxDiopsid_Expt_TotTCon(temperatures.size());
  constexpr int GrtPyropes_ExptID = 10;
  std::vector<double> GrtPyropes_Expt_TotTCon(temperatures.size());
  constexpr int GrtGrossul_ExptID = 11;
  std::vector<double> GrtGrossul_Expt_TotTCon(temperatures.size());
  constexpr int GrtAlmandi_ExptID = 12;
  std::vector<double> GrtAlmandi_Expt_TotTCon(temperatures.size());
  constexpr int GrtMajorit_ExptID = 13;
  std::vector<double> GrtMajorit_Expt_TotTCon(temperatures.size());
  constexpr int QuartzPure_ExptID = 14;
  std::vector<double> QuartzPure_Expt_TotTCon(temperatures.size());
  constexpr int CoesitSiO2_ExptID = 15;
  std::vector<double> CoesitSiO2_Expt_TotTCon(temperatures.size());
  constexpr int Stishovite_ExptID = 16;
  std::vector<double> Stishovite_Expt_TotTCon(temperatures.size());
  constexpr int Al05Stisho_ExptID = 17;
  std::vector<double> Al05Stisho_Expt_TotTCon(temperatures.size());
  constexpr int Antigor010_ExptID = 18;
  std::vector<double> Antigor010_Expt_TotTCon(temperatures.size());
  constexpr int Antigor001_ExptID = 19;
  std::vector<double> Antigor001_Expt_TotTCon(temperatures.size());
  constexpr int FeAlPhaseD_ExptID = 20;
  std::vector<double> FeAlPhaseD_Expt_TotTCon(temperatures.size());
  constexpr int Al02PhaseD_ExptID = 21;
  std::vector<double> Al02PhaseD_Expt_TotTCon(temperatures.size());
  constexpr int Ferroper08_ExptID = 22;
  std::vector<double> Ferroper08_Expt_TotTCon(temperatures.size());
  constexpr int Ferroper10_ExptID = 23;
  std::vector<double> Ferroper10_Expt_TotTCon(temperatures.size());
  constexpr int Ferroper20_ExptID = 24;
  std::vector<double> Ferroper20_Expt_TotTCon(temperatures.size());
  constexpr int Ferroper56_ExptID = 25;
  std::vector<double> Ferroper56_Expt_TotTCon(temperatures.size());
  constexpr int Davemaoite_ExptID = 26;
  std::vector<double> Davemaoite_Expt_TotTCon(temperatures.size());
  constexpr int NewHexAlPh_ExptID = 27;
  std::vector<double> NewHexAlPh_Expt_TotTCon(temperatures.size());
  constexpr int Akimotoite_ExptID = 28;
  std::vector<double> Akimotoite_Expt_TotTCon(temperatures.size());

  unsigned int MineralPar_Index = Akimotoite_ExptID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> Expt_Minerals_LatTcond(MineralPar_Index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> Expt_Minerals_RadTcond(MineralPar_Index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> Expt_Minerals_TotTcond(MineralPar_Index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity


  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> OlivineDry_Expt_LatTCon = {3.58888, 1.58719, 2.36282, 3.67868, 4.25767};
  std::vector<double> OlivineDry_Expt_RadTCon = {0.00138288, 2.23152, 2.34978, 2.45486, 3.12273};
  Expt_Minerals_LatTcond[OlivineDry_ExptID] = OlivineDry_Expt_LatTCon;
  Expt_Minerals_RadTcond[OlivineDry_ExptID] = OlivineDry_Expt_RadTCon;
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> WadsleyDry_Expt_LatTCon = {5.88364, 3.08196, 3.21488, 3.22817, 2.77368};
  std::vector<double> WadsleyDry_Expt_RadTCon = {1.1834e-9, 1.51515, 1.71067, 1.87995, 2.71457};
  Expt_Minerals_LatTcond[WadsleyDry_ExptID] = WadsleyDry_Expt_LatTCon;
  Expt_Minerals_RadTcond[WadsleyDry_ExptID] = WadsleyDry_Expt_RadTCon;
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> RingwooDry_Expt_LatTCon = {4.98456, 2.17063, 2.41846, 3.19795, 5.89102};
  std::vector<double> RingwooDry_Expt_RadTCon = {5.52667e-10, 0.74367, 0.85004, 0.94227, 1.38818};
  Expt_Minerals_LatTcond[RingwooDry_ExptID] = RingwooDry_Expt_LatTCon;
  Expt_Minerals_RadTcond[RingwooDry_ExptID] = RingwooDry_Expt_RadTCon;
  // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> En100Brigm_Expt_LatTCon = {10.69504, 2.03810, 2.21218, 2.47430, 5.69003};
  std::vector<double> En100Brigm_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[En100Brigm_ExptID] = En100Brigm_Expt_LatTCon;
  Expt_Minerals_RadTcond[En100Brigm_ExptID] = En100Brigm_Expt_RadTCon;
  // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> En97Brigma_Expt_LatTCon = {5.73751, 2.29118, 2.48780, 2.75655, 6.95242};
  std::vector<double> En97Brigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[En97Brigma_ExptID] = En97Brigma_Expt_LatTCon;
  Expt_Minerals_RadTcond[En97Brigma_ExptID] = En97Brigma_Expt_RadTCon;
  // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> En90Brigma_Expt_LatTCon = {3.79122, 2.91572, 3.14193, 3.43191, 8.92076};
  std::vector<double> En90Brigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[En90Brigma_ExptID] = En90Brigma_Expt_LatTCon;
  Expt_Minerals_RadTcond[En90Brigma_ExptID] = En90Brigma_Expt_RadTCon;
  // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> AlMgBrigma_Expt_LatTCon = {6.30403, 2.31503, 2.56860, 2.91357, 7.88626};
  std::vector<double> AlMgBrigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[AlMgBrigma_ExptID] = AlMgBrigma_Expt_LatTCon;
  Expt_Minerals_RadTcond[AlMgBrigma_ExptID] = AlMgBrigma_Expt_RadTCon;
  // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> FeAlBrigma_Expt_LatTCon = {3.99962, 1.87753, 2.05369, 2.30265, 6.17198};
  std::vector<double> FeAlBrigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[FeAlBrigma_ExptID] = FeAlBrigma_Expt_LatTCon;
  Expt_Minerals_RadTcond[FeAlBrigma_ExptID] = FeAlBrigma_Expt_RadTCon;
  // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> OpxEnstati_Expt_LatTCon = {5.79950, 2.55142, 3.15043, 3.26008, 2.56641};
  std::vector<double> OpxEnstati_Expt_RadTCon = {5.45056e-5, 2.93307, 3.08161, 3.20925, 3.90772};
  Expt_Minerals_LatTcond[OpxEnstati_ExptID] = OpxEnstati_Expt_LatTCon;
  Expt_Minerals_LatTcond[OpxEnstati_ExptID] = OpxEnstati_Expt_LatTCon;
  // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> CpxDiopsid_Expt_LatTCon = {5.99273, 2.65322, 3.13819, 3.61301, 3.41681};
  std::vector<double> CpxDiopsid_Expt_RadTCon = {2.64987e-5, 2.99817, 3.14622, 3.27243, 3.94085};
  Expt_Minerals_LatTcond[CpxDiopsid_ExptID] = CpxDiopsid_Expt_LatTCon;
  Expt_Minerals_LatTcond[CpxDiopsid_ExptID] = CpxDiopsid_Expt_LatTCon;
  // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> GrtPyropes_Expt_LatTCon = {4.38827, 2.15733, 2.64684, 3.54118, 4.22405};
  std::vector<double> GrtPyropes_Expt_RadTCon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  Expt_Minerals_LatTcond[GrtPyropes_ExptID] = GrtPyropes_Expt_LatTCon;
  Expt_Minerals_LatTcond[GrtPyropes_ExptID] = GrtPyropes_Expt_LatTCon;
  // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> GrtGrossul_Expt_LatTCon = {4.08838, 1.91337, 2.26557, 3.05964, 4.01413};
  std::vector<double> GrtGrossul_Expt_RadTCon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  Expt_Minerals_LatTcond[GrtGrossul_ExptID] = GrtGrossul_Expt_LatTCon;
  Expt_Minerals_LatTcond[GrtGrossul_ExptID] = GrtGrossul_Expt_LatTCon;
  // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> GrtAlmandi_Expt_LatTCon = {3.39124, 1.70813, 2.17978, 3.13522, 4.07543};
  std::vector<double> GrtAlmandi_Expt_RadTCon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  Expt_Minerals_LatTcond[GrtAlmandi_ExptID] = GrtAlmandi_Expt_LatTCon;
  Expt_Minerals_LatTcond[GrtAlmandi_ExptID] = GrtAlmandi_Expt_LatTCon;
  // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
  std::vector<double> GrtMajorit_Expt_LatTCon = {9.73983, 4.24079, 4.57076, 5.13035, 4.76249};
  std::vector<double> GrtMajorit_Expt_RadTCon = {2.55667e-5, 2.08101, 2.25936, 2.42089, 3.48797};
  Expt_Minerals_LatTcond[GrtMajorit_ExptID] = GrtMajorit_Expt_LatTCon;
  Expt_Minerals_LatTcond[GrtMajorit_ExptID] = GrtMajorit_Expt_LatTCon;
  // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> QuartzPure_Expt_LatTCon = {9.53243, 1.84339, 2.49820, 2.47676, 1.49323};
  std::vector<double> QuartzPure_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[QuartzPure_ExptID] = QuartzPure_Expt_LatTCon;
  Expt_Minerals_LatTcond[QuartzPure_ExptID] = QuartzPure_Expt_LatTCon;
  // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> CoesitSiO2_Expt_LatTCon = {7.21196, 1.31776, 1.23920, 1.17013, 0.84983};
  std::vector<double> CoesitSiO2_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[CoesitSiO2_ExptID] = CoesitSiO2_Expt_LatTCon;
  Expt_Minerals_LatTcond[CoesitSiO2_ExptID] = CoesitSiO2_Expt_LatTCon;
  // Stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Stishovite_Expt_LatTCon = {67.68617, 29.30897, 28.43361, 27.62648, 29.60473};
  std::vector<double> Stishovite_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[Stishovite_ExptID] = Stishovite_Expt_LatTCon;
  Expt_Minerals_LatTcond[Stishovite_ExptID] = Stishovite_Expt_LatTCon;
  // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Al05Stisho_Expt_LatTCon = {24.18571, 10.48883, 10.35956, 10.44473, 15.11569};
  std::vector<double> Al05Stisho_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[Al05Stisho_ExptID] = Al05Stisho_Expt_LatTCon;
  Expt_Minerals_LatTcond[Al05Stisho_ExptID] = Al05Stisho_Expt_LatTCon;
  // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Antigor010_Expt_LatTCon = {4.55589, 1.99618, 2.41198, 3.15852, 3.57446};
  std::vector<double> Antigor010_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[Antigor010_ExptID] = Antigor010_Expt_LatTCon;
  Expt_Minerals_LatTcond[Antigor010_ExptID] = Antigor010_Expt_LatTCon;
  // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Antigor001_Expt_LatTCon = {1.06669, 0.49210, 1.01828, 1.51143, 1.48578};
  std::vector<double> Antigor001_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[Antigor001_ExptID] = Antigor001_Expt_LatTCon;
  Expt_Minerals_LatTcond[Antigor001_ExptID] = Antigor001_Expt_LatTCon;
  // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> FeAlPhaseD_Expt_LatTCon = {2.59325, 1.13915, 1.31955, 1.59985, 6.12649};
  std::vector<double> FeAlPhaseD_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[FeAlPhaseD_ExptID] = FeAlPhaseD_Expt_LatTCon;
  Expt_Minerals_LatTcond[FeAlPhaseD_ExptID] = FeAlPhaseD_Expt_LatTCon;
  // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Al02PhaseD_Expt_LatTCon = {3.60666, 1.56888, 1.65207, 1.95489, 8.61776};
  std::vector<double> Al02PhaseD_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  Expt_Minerals_LatTcond[Al02PhaseD_ExptID] = Al02PhaseD_Expt_LatTCon;
  Expt_Minerals_LatTcond[Al02PhaseD_ExptID] = Al02PhaseD_Expt_LatTCon;
  // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Ferroper08_Expt_LatTCon = {5.08425, 2.20657, 2.24939, 2.50821, 14.39941};
  std::vector<double> Ferroper08_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Ferroper08_ExptID] = Ferroper08_Expt_LatTCon;
  Expt_Minerals_LatTcond[Ferroper08_ExptID] = Ferroper08_Expt_LatTCon;
  // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Ferroper10_Expt_LatTCon = {4.48610, 1.94695, 1.98084, 2.19296, 12.65637};
  std::vector<double> Ferroper10_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Ferroper10_ExptID] = Ferroper10_Expt_LatTCon;
  Expt_Minerals_LatTcond[Ferroper10_ExptID] = Ferroper10_Expt_LatTCon;
  // Ferropericlase (Mg80Fe20O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Ferroper20_Expt_LatTCon = {3.48647, 3.39116, 3.56449, 3.77742, 7.57324};
  std::vector<double> Ferroper20_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Ferroper20_ExptID] = Ferroper20_Expt_LatTCon;
  Expt_Minerals_LatTcond[Ferroper20_ExptID] = Ferroper20_Expt_LatTCon;
  // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> Ferroper56_Expt_LatTCon = {2.69167, 1.23171, 1.55079, 2.02398, 7.04069};
  std::vector<double> Ferroper56_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Ferroper56_ExptID] = Ferroper56_Expt_LatTCon;
  Expt_Minerals_LatTcond[Ferroper56_ExptID] = Ferroper56_Expt_LatTCon;
  // Davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> Davemaoite_Expt_LatTCon = {10.86634, 4.77342, 5.15289, 5.86997, 13.48443};
  std::vector<double> Davemaoite_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Davemaoite_ExptID] = Davemaoite_Expt_LatTCon;
  Expt_Minerals_LatTcond[Davemaoite_ExptID] = Davemaoite_Expt_LatTCon;
  // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> NewHexAlPh_Expt_LatTCon = {10.59581, 4.58812, 4.45113, 4.32578, 12.15184};
  std::vector<double> NewHexAlPh_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[NewHexAlPh_ExptID] = NewHexAlPh_Expt_LatTCon;
  Expt_Minerals_LatTcond[NewHexAlPh_ExptID] = NewHexAlPh_Expt_LatTCon;
  // Akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> Akimotoite_Expt_LatTCon = {10.69504, 2.03810, 2.21218, 2.47430, 5.69003};
  std::vector<double> Akimotoite_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  Expt_Minerals_LatTcond[Akimotoite_ExptID] = Akimotoite_Expt_LatTCon;
  Expt_Minerals_LatTcond[Akimotoite_ExptID] = Akimotoite_Expt_LatTCon;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    OlivineDry_Expt_TotTCon[row] = OlivineDry_Expt_LatTCon[row] + OlivineDry_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[OlivineDry_ExptID] = OlivineDry_Expt_TotTCon;
    WadsleyDry_Expt_TotTCon[row] = WadsleyDry_Expt_LatTCon[row] + WadsleyDry_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[WadsleyDry_ExptID] = WadsleyDry_Expt_TotTCon;
    RingwooDry_Expt_TotTCon[row] = RingwooDry_Expt_LatTCon[row] + RingwooDry_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[RingwooDry_ExptID] = RingwooDry_Expt_TotTCon;
    En100Brigm_Expt_TotTCon[row] = En100Brigm_Expt_LatTCon[row] + En100Brigm_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[En100Brigm_ExptID] = En100Brigm_Expt_TotTCon;
    En97Brigma_Expt_TotTCon[row] = En97Brigma_Expt_LatTCon[row] + En97Brigma_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[En97Brigma_ExptID] = En97Brigma_Expt_TotTCon;
    En90Brigma_Expt_TotTCon[row] = En90Brigma_Expt_LatTCon[row] + En90Brigma_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[En90Brigma_ExptID] = En90Brigma_Expt_TotTCon;
    AlMgBrigma_Expt_TotTCon[row] = AlMgBrigma_Expt_LatTCon[row] + AlMgBrigma_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[AlMgBrigma_ExptID] = AlMgBrigma_Expt_TotTCon;
    FeAlBrigma_Expt_TotTCon[row] = FeAlBrigma_Expt_LatTCon[row] + FeAlBrigma_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[FeAlBrigma_ExptID] = FeAlBrigma_Expt_TotTCon;
    OpxEnstati_Expt_TotTCon[row] = OpxEnstati_Expt_LatTCon[row] + OpxEnstati_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[OpxEnstati_ExptID] = OpxEnstati_Expt_TotTCon;
    CpxDiopsid_Expt_TotTCon[row] = CpxDiopsid_Expt_LatTCon[row] + CpxDiopsid_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[CpxDiopsid_ExptID] = CpxDiopsid_Expt_TotTCon;
    GrtPyropes_Expt_TotTCon[row] = GrtPyropes_Expt_LatTCon[row] + GrtPyropes_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[GrtPyropes_ExptID] = GrtPyropes_Expt_TotTCon;
    GrtGrossul_Expt_TotTCon[row] = GrtGrossul_Expt_LatTCon[row] + GrtGrossul_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[GrtGrossul_ExptID] = GrtGrossul_Expt_TotTCon;
    GrtAlmandi_Expt_TotTCon[row] = GrtAlmandi_Expt_LatTCon[row] + GrtAlmandi_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[GrtAlmandi_ExptID] = GrtAlmandi_Expt_TotTCon;
    GrtMajorit_Expt_TotTCon[row] = GrtMajorit_Expt_LatTCon[row] + GrtMajorit_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[GrtMajorit_ExptID] = GrtMajorit_Expt_TotTCon;
    QuartzPure_Expt_TotTCon[row] = QuartzPure_Expt_LatTCon[row] + QuartzPure_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[QuartzPure_ExptID] = QuartzPure_Expt_TotTCon;
    CoesitSiO2_Expt_TotTCon[row] = CoesitSiO2_Expt_LatTCon[row] + CoesitSiO2_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[CoesitSiO2_ExptID] = CoesitSiO2_Expt_TotTCon;
    Stishovite_Expt_TotTCon[row] = Stishovite_Expt_LatTCon[row] + Stishovite_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Stishovite_ExptID] = Stishovite_Expt_TotTCon;
    Al05Stisho_Expt_TotTCon[row] = Al05Stisho_Expt_LatTCon[row] + Al05Stisho_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Al05Stisho_ExptID] = Al05Stisho_Expt_TotTCon;
    Antigor010_Expt_TotTCon[row] = Antigor010_Expt_LatTCon[row] + Antigor010_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Antigor010_ExptID] = Antigor010_Expt_TotTCon;
    Antigor001_Expt_TotTCon[row] = Antigor001_Expt_LatTCon[row] + Antigor001_Expt_RadTCon[row]; 
    Expt_Minerals_TotTcond[Antigor001_ExptID] = Antigor001_Expt_TotTCon;   
    FeAlPhaseD_Expt_TotTCon[row] = FeAlPhaseD_Expt_LatTCon[row] + FeAlPhaseD_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[FeAlPhaseD_ExptID] = FeAlPhaseD_Expt_TotTCon; 
    Al02PhaseD_Expt_TotTCon[row] = Al02PhaseD_Expt_LatTCon[row] + Al02PhaseD_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Al02PhaseD_ExptID] = Al02PhaseD_Expt_TotTCon;
    Ferroper08_Expt_TotTCon[row] = Ferroper08_Expt_LatTCon[row] + Ferroper08_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Ferroper08_ExptID] = Ferroper08_Expt_TotTCon;
    Ferroper10_Expt_TotTCon[row] = Ferroper10_Expt_LatTCon[row] + Ferroper10_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Ferroper10_ExptID] = Ferroper10_Expt_TotTCon;
    Ferroper20_Expt_TotTCon[row] = Ferroper20_Expt_LatTCon[row] + Ferroper20_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Ferroper20_ExptID] = Ferroper20_Expt_TotTCon;
    Ferroper56_Expt_TotTCon[row] = Ferroper56_Expt_LatTCon[row] + Ferroper56_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Ferroper56_ExptID] = Ferroper56_Expt_TotTCon;
    Davemaoite_Expt_TotTCon[row] = Davemaoite_Expt_LatTCon[row] + Davemaoite_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Davemaoite_ExptID] = Davemaoite_Expt_TotTCon;
    NewHexAlPh_Expt_TotTCon[row] = NewHexAlPh_Expt_LatTCon[row] + NewHexAlPh_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[NewHexAlPh_ExptID] = NewHexAlPh_Expt_TotTCon;
    Akimotoite_Expt_TotTCon[row] = Akimotoite_Expt_LatTCon[row] + Akimotoite_Expt_RadTCon[row];
    Expt_Minerals_TotTcond[Akimotoite_ExptID] = Akimotoite_Expt_TotTCon;
  }

 // Loop over all mID values
 for (unsigned int mID = 0; mID < MineralPar_Index; ++mID)
 {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_total_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
   {
     for (size_t col = 0; col < compositions[row].size(); ++col)
     {
       expected_total_Tcond[row][col] = std::pow(Expt_Minerals_TotTcond[mID][col], compositions[row][col]);
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
         case OlivineDry_ExptID: // OlivineDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("OlivineDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("OlivineDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case WadsleyDry_ExptID: // WadsleyDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("WadsleyDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("WadsleyDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case RingwooDry_ExptID: // RingwooDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("RingwooDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("RingwooDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case En100Brigm_ExptID: // En100Brigm
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("En100Brigm Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("En100Brigm Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case En97Brigma_ExptID: // En97Brigma
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("En97Brigma Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("En97Brigma Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case En90Brigma_ExptID: // En90Brigma
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("En90Brigma Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("En90Brigma Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case AlMgBrigma_ExptID: // AlMgBrigma
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("AlMgBrigma Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("AlMgBrigma Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case FeAlBrigma_ExptID: // FeAlBrigma
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("FeAlBrigma Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("FeAlBrigma Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case OpxEnstati_ExptID: // OpxEnstati
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("OpxEnstati Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("OpxEnstati Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case CpxDiopsid_ExptID: // CpxDiopsid
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("CpxDiopsid Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("CpxDiopsid Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case GrtPyropes_ExptID: // GrtPyropes
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("GrtPyropes Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("GrtPyropes Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case GrtGrossul_ExptID: // GrtGrossul
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("GrtGrossul Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("GrtGrossul Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case GrtAlmandi_ExptID: // GrtAlmandi
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("GrtAlmandi Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("GrtAlmandi Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case GrtMajorit_ExptID: // GrtMajorit
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("GrtMajorit Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("GrtMajorit Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case QuartzPure_ExptID: // QuartzPure
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("QuartzPure Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("QuartzPure Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case CoesitSiO2_ExptID: // CoesitSiO2
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("CoesitSiO2 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("CoesitSiO2 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Stishovite_ExptID: // Stishovite
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Stishovite Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Stishovite Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Al05Stisho_ExptID: // Al05Stisho
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Al05Stisho Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Al05Stisho Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Antigor010_ExptID: // Antigor010
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Antigor010 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Antigor010 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Antigor001_ExptID: // Antigor001
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Antigor001 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Antigor001 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case FeAlPhaseD_ExptID: // FeAlPhaseD
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("FeAlPhaseD Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("FeAlPhaseD Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Al02PhaseD_ExptID: // Al02PhaseD
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Al02PhaseD Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Al02PhaseD Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Ferroper08_ExptID: // Ferroper08
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Ferroper08 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Ferroper08 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Ferroper10_ExptID: // Ferroper10
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Ferroper10 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Ferroper10 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Ferroper20_ExptID: // Ferroper20
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Ferroper20 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Ferroper20 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Ferroper56_ExptID: // Ferroper56
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Ferroper56 Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Ferroper56 Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Davemaoite_ExptID: // Davemaoite
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Davemaoite Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Davemaoite Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case NewHexAlPh_ExptID: // NewHexAlPh
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("NewHexAlPh Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("NewHexAlPh Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case Akimotoite_ExptID: // Akimotoite
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("Akimotoite Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("Akimotoite Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
       } 
     }
   }
 } 
}

TEST_CASE("Utilities:: Thermal Conductivity Hofmeister 1999")
{
  aspect::MaterialModel::ThermalConductivity::Hofmeister1999<3> model;
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
  constexpr int OlivineDry_Hof99ID = 0;
  std::vector<double> OlivineDry_Expt_Hof99_TotTCon(temperatures.size());
  constexpr int WadsleyDry_Hof99ID = 1;
  std::vector<double> WadsleyDry_Expt_Hof99_TotTCon(temperatures.size());
  constexpr int RingwooDry_Hof99ID = 2;
  std::vector<double> RingwooDry_Expt_Hof99_TotTCon(temperatures.size());

  unsigned int Hofm1999Par_Index = RingwooDry_Hof99ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> Expt_Hofm1999_LatTcond(Hofm1999Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> Expt_Hofm1999_RadTcond(Hofm1999Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Radiative thermal conductivity
  std::vector<std::vector<double>> Expt_Hofm1999_TotTcond(Hofm1999Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> OlivineDry_Expt_Hof99_LatTCon = {4.69000, 2.41949, 2.66506, 2.97359, 7.19819};
  std::vector<double> OlivineDry_Expt_Hof99_RadTCon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  Expt_Hofm1999_LatTcond[OlivineDry_Hof99ID] = OlivineDry_Expt_Hof99_LatTCon;
  Expt_Hofm1999_RadTcond[OlivineDry_Hof99ID] = OlivineDry_Expt_Hof99_RadTCon;
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> WadsleyDry_Expt_Hof99_LatTCon = {7.68434, 4.20366, 4.53114, 4.95031, 11.17803};
  std::vector<double> WadsleyDry_Expt_Hof99_RadTCon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  Expt_Hofm1999_LatTcond[WadsleyDry_Hof99ID] = WadsleyDry_Expt_Hof99_LatTCon;
  Expt_Hofm1999_RadTcond[WadsleyDry_Hof99ID] = WadsleyDry_Expt_Hof99_RadTCon;
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  std::vector<double> RingwooDry_Expt_Hof99_LatTCon = {7.68419, 4.14794, 4.47408, 4.89128, 11.00724};
  std::vector<double> RingwooDry_Expt_Hof99_RadTCon = {0.00572, 0.28688, 0.32277, 0.35967, 0.80725};
  Expt_Hofm1999_LatTcond[RingwooDry_Hof99ID] = RingwooDry_Expt_Hof99_LatTCon;
  Expt_Hofm1999_RadTcond[RingwooDry_Hof99ID] = RingwooDry_Expt_Hof99_RadTCon;

 // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    OlivineDry_Expt_Hof99_TotTCon[row] = OlivineDry_Expt_Hof99_LatTCon[row] + OlivineDry_Expt_Hof99_RadTCon[row];
    Expt_Hofm1999_TotTcond[OlivineDry_Hof99ID] = OlivineDry_Expt_Hof99_TotTCon;
    WadsleyDry_Expt_Hof99_TotTCon[row] = WadsleyDry_Expt_Hof99_LatTCon[row] + WadsleyDry_Expt_Hof99_RadTCon[row];
    Expt_Hofm1999_TotTcond[WadsleyDry_Hof99ID] = WadsleyDry_Expt_Hof99_TotTCon;
    RingwooDry_Expt_Hof99_TotTCon[row] = RingwooDry_Expt_Hof99_LatTCon[row] + RingwooDry_Expt_Hof99_RadTCon[row];
    Expt_Hofm1999_TotTcond[RingwooDry_Hof99ID] = RingwooDry_Expt_Hof99_TotTCon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < Hofm1999Par_Index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_Hof1999_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_Hof1999_Tcond[row][col] = std::pow(Expt_Hofm1999_TotTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_Hof1999_Tcond;

   INFO("Checking Hof1999 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

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
         case OlivineDry_Hof99ID: // OlivineDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("OlivineDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("OlivineDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case WadsleyDry_Hof99ID: // WadsleyDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("WadsleyDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("WadsleyDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case RingwooDry_Hof99ID: // RingwooDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("RingwooDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("RingwooDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
        } 
      }
    }
  }
}

/*
TEST_CASE("Utilities:: Thermal Conductivity Anderson 1987")
{
  aspect::MaterialModel::ThermalConductivity::Anderson1987<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/

/*
TEST_CASE("Utilities:: Thermal Conductivity Gerya 2021")
{
  aspect::MaterialModel::ThermalConductivity::Gerya2021<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/

/*
TEST_CASE("Utilities:: Thermal Conductivity Grose & Afonso 2019")
{
  aspect::MaterialModel::ThermalConductivity::GroseAfonso2019<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/

/*
TEST_CASE("Utilities:: Thermal Conductivity Hofmeister 2005")
{
  aspect::MaterialModel::ThermalConductivity::Hofmeister2005<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/

/*
TEST_CASE("Utilities:: Thermal Conductivity Hofmeister & Branlund 2015")
{
  aspect::MaterialModel::ThermalConductivity::HofmeisterBranlund2015<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/

/*
TEST_CASE("Utilities:: Thermal Conductivity Stackhouse 2015")
{
  aspect::MaterialModel::ThermalConductivity::Stackhouse2015<3> model;
  aspect::MaterialModel::MaterialModelInputs<3> in(5,1);    // Adjust the size of inputs as needed
  aspect::MaterialModel::MaterialModelOutputs<3> out(5,1);  // Adjust the size of outputs as needed
}
*/


TEST_CASE("Utilities:: Thermal Conductivity Tosi 2016")
{
  aspect::MaterialModel::ThermalConductivity::Tosi2016<3> model;
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
  constexpr int UpperMantl_Tos16ID = 0;
  std::vector<double> UpperMantl_Expt_Tos16_TotTCon(temperatures.size());
  constexpr int UManTraZon_Tos16ID = 1;
  std::vector<double> UManTraZon_Expt_Tos16_TotTCon(temperatures.size());
  constexpr int LManTraZon_Tos16ID = 2;
  std::vector<double> LManTraZon_Expt_Tos16_TotTCon(temperatures.size());
  constexpr int LowerMantl_Tos16ID = 3;
  std::vector<double> LowerMantl_Expt_Tos16_TotTCon(temperatures.size());

  unsigned int Tosi2016Par_Index = LowerMantl_Tos16ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> Expt_Tos16_LatTcond(Tosi2016Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> Expt_Tos16_TotTcond(Tosi2016Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Upper Mantle: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> UpperMantl_Expt_Tos16_LatTCon = {2.46271, 1.24999, 1.78653, 2.43429, 11.71039};
  Expt_Tos16_LatTcond[UpperMantl_Tos16ID] = UpperMantl_Expt_Tos16_LatTCon;
  // Upper Mantle Transition Zone: expected lattice thermal conductivities (k) in [W/m/K]
  std::vector<double> UManTraZon_Expt_Tos16_LatTCon = {3.79686, 1.61966, 2.07866, 2.63431, 10.37773};
  Expt_Tos16_LatTcond[UManTraZon_Tos16ID] = UManTraZon_Expt_Tos16_LatTCon;
  // Lower Mantle Transition Zone: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> LManTraZon_Expt_Tos16_LatTCon = {3.50678, 1.39227, 1.83968, 2.37777, 9.66447};
  Expt_Tos16_LatTcond[LManTraZon_Tos16ID] = LManTraZon_Expt_Tos16_LatTCon;
  // Lower Mantle: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> LowerMantl_Expt_Tos16_LatTCon = {3.47335, 2.13845, 2.37846, 2.68032, 7.56725};
  Expt_Tos16_LatTcond[LowerMantl_Tos16ID] = LowerMantl_Expt_Tos16_LatTCon;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    UpperMantl_Expt_Tos16_TotTCon[row] = UpperMantl_Expt_Tos16_LatTCon[row];
    Expt_Tos16_TotTcond[UpperMantl_Tos16ID] = UpperMantl_Expt_Tos16_TotTCon;
    UManTraZon_Expt_Tos16_TotTCon[row] = UManTraZon_Expt_Tos16_LatTCon[row];
    Expt_Tos16_TotTcond[UManTraZon_Tos16ID] = UManTraZon_Expt_Tos16_TotTCon;
    LManTraZon_Expt_Tos16_TotTCon[row] = LManTraZon_Expt_Tos16_LatTCon[row];
    Expt_Tos16_TotTcond[LManTraZon_Tos16ID] = LManTraZon_Expt_Tos16_TotTCon;
    LowerMantl_Expt_Tos16_TotTCon[row] = LowerMantl_Expt_Tos16_LatTCon[row];
    Expt_Tos16_TotTcond[LowerMantl_Tos16ID] = LowerMantl_Expt_Tos16_TotTCon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < Tosi2016Par_Index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_Tos16_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_Tos16_Tcond[row][col] = std::pow(Expt_Tos16_TotTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_Tos16_Tcond;

   INFO("Checking Tosi2016 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

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
         case UpperMantl_Tos16ID: // Upper Mantle (UM)
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("UpperMantl Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("UpperMantl Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case UManTraZon_Tos16ID: // Upper Mantle Transition Zone (uMTZ)
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("UManTraZon Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("UManTraZon Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case LManTraZon_Tos16ID: // Lower Mantle Transition Zone (lMTZ)
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("LManTraZon Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("LManTraZon Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case LowerMantl_Tos16ID: // Lower Mantle (LM)
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("LowerMantl Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("LowerMantl Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
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
  aspect::MaterialModel::ThermalConductivity::Xu2004<3> model;
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
  constexpr int OlivineDry_Xu004ID = 0;
  std::vector<double> OlivineDry_Expt_Xu004_TotTCon(temperatures.size());
  constexpr int WadsleyDry_Xu004ID = 1;
  std::vector<double> WadsleyDry_Expt_Xu004_TotTCon(temperatures.size());
  constexpr int RingwooDry_Xu004ID = 2;
  std::vector<double> RingwooDry_Expt_Xu004_TotTCon(temperatures.size());

  unsigned int Xu2004Par_Index = RingwooDry_Xu004ID+1; // Number of minerals

  // Preallocate matrixes for storing thermal conductivities of minerals
  std::vector<std::vector<double>> Expt_Xu2004_LatTcond(Xu2004Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Lattice thermal conductivity
  std::vector<std::vector<double>> Expt_Xu2004_TotTcond(Xu2004Par_Index, std::vector<double>(temperatures.size(), 0.0)); // Total thermal conductivity

  // Dry Olivine: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> OlivineDry_Expt_Xu004_LatTCon = {4.11726, 1.83987, 2.00632, 2.21873, 5.46835};
  Expt_Xu2004_LatTcond[OlivineDry_Xu004ID] = OlivineDry_Expt_Xu004_LatTCon;
  // Dry Wadsleyite: expected lattice thermal conductivities (k) in [W/m/K]
  std::vector<double> WadsleyDry_Expt_Xu004_LatTCon = {8.07501, 3.57699, 3.78227, 4.05482, 8.42667};
  Expt_Xu2004_LatTcond[WadsleyDry_Xu004ID] = WadsleyDry_Expt_Xu004_LatTCon;
  // Dry Ringwoodite: expected lattice thermal conductivities (k) in [W/m/K] 
  std::vector<double> RingwooDry_Expt_Xu004_LatTCon = {9.51056, 4.20878, 4.43470, 4.73685, 9.62399};
  Expt_Xu2004_LatTcond[RingwooDry_Xu004ID] = RingwooDry_Expt_Xu004_LatTCon;

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
    OlivineDry_Expt_Xu004_TotTCon[row] = OlivineDry_Expt_Xu004_LatTCon[row];
    Expt_Xu2004_TotTcond[OlivineDry_Xu004ID] = OlivineDry_Expt_Xu004_TotTCon;
    WadsleyDry_Expt_Xu004_TotTCon[row] = WadsleyDry_Expt_Xu004_LatTCon[row];
    Expt_Xu2004_TotTcond[WadsleyDry_Xu004ID] = WadsleyDry_Expt_Xu004_TotTCon;
    RingwooDry_Expt_Xu004_TotTCon[row] = RingwooDry_Expt_Xu004_LatTCon[row];
    Expt_Xu2004_TotTcond[RingwooDry_Xu004ID] = RingwooDry_Expt_Xu004_TotTCon;
  }

  // Loop over all mID values
  for (unsigned int mID = 0; mID < Xu2004Par_Index; ++mID)
  {
   in.Mineral_ID = mID; // Set the current mID

   // Initialize the expected value matrix with the same dimensions of the composition matrix
   std::vector<std::vector<double>> expected_Xu2004_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

   // Perform element-wise calculation
   for (size_t row = 0; row < compositions.size(); ++row)
    {
      for (size_t col = 0; col < compositions[row].size(); ++col)
      {
        expected_Xu2004_Tcond[row][col] = std::pow(Expt_Xu2004_TotTcond[mID][col], compositions[row][col]);
      }
    }

   std::vector<std::vector<double>> expected_conductivities = expected_Xu2004_Tcond;

   INFO("Checking Xu2004 thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X)");

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
         case OlivineDry_Xu004ID: // OlivineDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("OlivineDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("OlivineDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case WadsleyDry_Xu004ID: // WadsleyDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("WadsleyDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("WadsleyDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
           REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
           break;
         }
         case RingwooDry_Xu004ID: // RingwooDry
         {
           INFO("Conditions T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%]");
           INFO("RingwooDry Expected k= " << expected_conductivities[row][i] << "[W/m/K]");
           INFO("RingwooDry Computed k= " << out.thermal_conductivities[i] << "[W/m/K]");
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

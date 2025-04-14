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

  unsigned int MineralPar_Index = 0; // Initialize the counter

  // Preallocate the expected total thermal conductivities (k) in [W/m/K]
  unsigned int OlivineDry_ExptID = MineralPar_Index++;
  std::vector<double> OlivineDry_Expt_TotTCon(temperatures.size());
  // unsigned int WadsleyDry_ExptID = MineralPar_Index++;
  // std::vector<double> WadsleyDry_Expt_TotTCon(temperatures.size());
  // unsigned int RingwooDry_ExptID = MineralPar_Index++;
  // std::vector<double> RingwooDry_Expt_TotTCon(temperatures.size());
  // unsigned int En100Brigm_ExptID = MineralPar_Index++;
  // std::vector<double> En100Brigm_Expt_TotTCon(temperatures.size());
  // unsigned int En97Brigma_ExptID = MineralPar_Index++;
  // std::vector<double> En97Brigma_Expt_TotTCon(temperatures.size());
  // unsigned int En90Brigma_ExptID = MineralPar_Index++;
  // std::vector<double> En90Brigma_Expt_TotTCon(temperatures.size());
  // unsigned int AlMgBrigma_ExptID = MineralPar_Index++;
  // std::vector<double> AlMgBrigma_Expt_TotTCon(temperatures.size());
  // unsigned int FeAlBrigma_ExptID = MineralPar_Index++;
  // std::vector<double> FeAlBrigma_Expt_TotTCon(temperatures.size());
  // unsigned int OpxEnstati_ExptID = MineralPar_Index++;
  // std::vector<double> OpxEnstati_Expt_TotTCon(temperatures.size());
  // unsigned int CpxDiopsid_ExptID = MineralPar_Index++;
  // std::vector<double> CpxDiopsid_Expt_TotTCon(temperatures.size());
  // unsigned int GrtPyropes_ExptID = MineralPar_Index++;
  // std::vector<double> GrtPyropes_Expt_TotTCon(temperatures.size());
  // unsigned int GrtGrossul_ExptID = MineralPar_Index++;
  // std::vector<double> GrtGrossul_Expt_TotTCon(temperatures.size());
  // unsigned int GrtAlmandi_ExptID = MineralPar_Index++;
  // std::vector<double> GrtAlmandi_Expt_TotTCon(temperatures.size());
  // unsigned int GrtMajorit_ExptID = MineralPar_Index++;
  // std::vector<double> GrtMajorit_Expt_TotTCon(temperatures.size());
  // unsigned int QuartzPure_ExptID = MineralPar_Index++;
  // std::vector<double> QuartzPure_Expt_TotTCon(temperatures.size());
  // unsigned int CoesitSiO2_ExptID = MineralPar_Index++;
  // std::vector<double> CoesitSiO2_Expt_TotTCon(temperatures.size());
  // unsigned int Stishovite_ExptID = MineralPar_Index++;
  // std::vector<double> Stishovite_Expt_TotTCon(temperatures.size());
  // unsigned int Al05Stisho_ExptID = MineralPar_Index++;
  // std::vector<double> Al05Stisho_Expt_TotTCon(temperatures.size());
  // unsigned int Antigor010_ExptID = MineralPar_Index++;
  // std::vector<double> Antigor010_Expt_TotTCon(temperatures.size());
  // unsigned int Antigor001_ExptID = MineralPar_Index++;
  // std::vector<double> Antigor001_Expt_TotTCon(temperatures.size());
  // unsigned int FeAlPhaseD_ExptID = MineralPar_Index++;
  // std::vector<double> FeAlPhaseD_Expt_TotTCon(temperatures.size());
  // unsigned int Al02PhaseD_ExptID = MineralPar_Index++;
  // std::vector<double> Al02PhaseD_Expt_TotTCon(temperatures.size());
  // unsigned int Ferroper08_ExptID = MineralPar_Index++;
  // std::vector<double> Ferroper08_Expt_TotTCon(temperatures.size());
  // unsigned int Ferroper10_ExptID = MineralPar_Index++;
  // std::vector<double> Ferroper10_Expt_TotTCon(temperatures.size());
  // unsigned int Ferroper56_ExptID = MineralPar_Index++;
  // std::vector<double> Ferroper56_Expt_TotTCon(temperatures.size());
  // unsigned int Davemaoite_ExptID = MineralPar_Index++;
  // std::vector<double> Davemaoite_Expt_TotTCon(temperatures.size());
  // unsigned int NewHexAlPh_ExptID = MineralPar_Index++;
  // std::vector<double> NewHexAlPh_Expt_TotTCon(temperatures.size());
  // unsigned int Akimotoite_ExptID = MineralPar_Index++;
  // std::vector<double> Akimotoite_Expt_TotTCon(temperatures.size());


  // Olivine: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  std::vector<double> OlivineDry_Expt_LatTCon = {3.58888, 1.58719, 2.36282, 3.67868, 4.25767};
  std::vector<double> OlivineDry_Expt_RadTCon = {0.00138288, 2.23152, 2.34978, 2.45486, 3.12273};
  // Dry Wadsleyite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> WadsleyDry_Expt_LatTCon = {5.88364, 4.17228, 3.21488, 3.22817, 2.77368};
  // std::vector<double> WadsleyDry_Expt_RadTCon = {1.1834e-9, 1.51515, 1.71067, 1.87995, 2.71457};
  // Dry Ringwoodite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  // std::vector<double> RingwooDry_Expt_LatTCon = {4.98456, 2.17063, 2.41846, 3.19795, 5.89102};
  // std::vector<double> RingwooDry_Expt_RadTCon = {5.52667e-10, 0.74367, 0.85004, 0.94227, 1.38818};
  // Mg-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> En100Brigm_Expt_LatTCon = {8.59248, 1.78976, 2.38264, 2.76606, 2.37440};
  // std::vector<double> En100Brigm_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (3%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> En97Brigma_Expt_LatTCon = {4.87331, 2.06731, 2.62291, 3.01271, 3.08418};
  // std::vector<double> En97Brigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe-Bridgmanite (10%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> En90Brigma_Expt_LatTCon = {3.47292, 2.76473, 3.27554, 3.66762, 4.33882};
  // std::vector<double> En90Brigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  // std::vector<double> AlMgBrigma_Expt_LatTCon = {3.46328, 1.22708, 1.18182, 1.14069, 0.83111};
  // std::vector<double> AlMgBrigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Fe,Al-Bridgmanite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> FeAlBrigma_Expt_LatTCon = {3.38945, 1.69568, 2.17986, 2.53249, 2.75512};
  // std::vector<double> FeAlBrigma_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Orthopyroxene (Enstatite): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> OpxEnstati_Expt_LatTCon = {5.79950, 2.55142, 3.15043, 3.26008, 2.56641};
  // std::vector<double> OpxEnstati_Expt_RadTCon = {5.45056e-5, 2.93307, 3.08161, 3.20925, 3.90772};
  // Clinopyroxene (Diopside): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> CpxDiopsid_Expt_LatTCon = {6.01129, 6.14634, 7.49354, 8.87744, 10.83836};
  // std::vector<double> CpxDiopsid_Expt_RadTCon = {2.64987e-5, 2.99817, 3.14622, 3.27243, 3.94085};
  // Garnet (Pyrope): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> GrtPyropes_Expt_LatTCon = {4.38827, 2.15733, 2.64684, 3.54118, 4.22405};
  // std::vector<double> GrtPyropes_Expt_RadTCon = {7.97451e-5, 2.20910, 2.37779, 2.52981, 3.52532};
  // Garnet (Grossular): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> GrtGrossul_Expt_LatTCon = {4.08838, 1.91337, 2.26557, 3.05964, 4.01413};
  // std::vector<double> GrtGrossul_Expt_RadTCon = {7.97451e-5, 2.20910, 2.37779, 2.52981, 3.52532};
  // Garnet (Almandine): expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  // std::vector<double> GrtAlmandi_Expt_LatTCon = {3.39124, 1.70813, 2.17978, 3.13522, 4.07543};
  // std::vector<double> GrtAlmandi_Expt_RadTCon = {7.97451e-5, 2.20910, 2.37779, 2.52981, 3.52532};
  // Garnet (Majorite): expected lattice and radiative thermal conductivities (k) in [W/m/K]  
  // std::vector<double> GrtMajorit_Expt_LatTCon = {9.73983, 4.24079, 4.57076, 5.13035, 4.76249};
  // std::vector<double> GrtMajorit_Expt_RadTCon = {7.97451e-5, 2.20910, 2.37779, 2.52981, 3.52532};
  // Quartz: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> QuartzPure_Expt_LatTCon = {9.53243, 1.84339, 2.49820, 2.47676, 1.49323};
  // std::vector<double> QuartzPure_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Coesite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> CoesitSiO2_Expt_LatTCon = {7.21196, 1.31776, 1.23920, 1.17013, 0.84983};
  // std::vector<double> CoesitSiO2_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Stishovite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Stishovite_Expt_LatTCon = {67.68617, 29.30897, 28.43361, 27.62648, 18.97792};
  // std::vector<double> Stishovite_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-stishovite (5 vol%): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Al05Stisho_Expt_LatTCon = {24.18571, 10.48883, 10.35956, 10.44473, 15.11569};
  // std::vector<double> Al05Stisho_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (010): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Antigor010_Expt_LatTCon = {4.55589, 1.99618, 2.41198, 3.15852, 3.57446};
  // std::vector<double> Antigor010_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Antigorite (001): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Antigor001_Expt_LatTCon = {1.06669, 0.49210, 1.01828, 1.51143, 1.48578};
  // std::vector<double> Antigor001_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Fe,Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> FeAlPhaseD_Expt_LatTCon = {2.59325, 1.13915, 1.31955, 1.59985, 1.78759};
  // std::vector<double> FeAlPhaseD_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Al-phase D (Dense Hydrous Magnesium Silicate): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Al02PhaseD_Expt_LatTCon = {3.60666, 1.56888, 1.65207, 1.95489, 8.61776};
  // std::vector<double> Al02PhaseD_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Ferropericlase (Mg92Fe8O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Ferroper08_Expt_LatTCon = {5.08425, 2.20657, 2.24939, 2.50821, 14.39941};
  // std::vector<double> Ferroper08_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Ferropericlase (Mg90Fe10O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Ferroper10_Expt_LatTCon = {4.49930, 4.32469, 4.52851, 5.15144, 37.89530};
  // std::vector<double> Ferroper10_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Ferropericlase (Mg56Fe44O): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Ferroper56_Expt_LatTCon = {2.69167, 1.23171, 1.55079, 2.02398, 7.04069};
  // std::vector<double> Ferroper56_Expt_RadTCon = {9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11, 9.87998e-11};
  // Davemaoite: expected lattice and radiative thermal conductivities (k) in [W/m/K] 
  // std::vector<double> Davemaoite_Expt_LatTCon = {10.86634, 4.77342, 5.15289, 5.86997, 13.48443};
  // std::vector<double> Davemaoite_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // New-hexagonal-alluminium-phase (FeNAL): expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> NewHexAlPh_Expt_LatTCon = {10.59581, 4.58812, 4.45113, 4.32578, 12.15184};
  // std::vector<double> NewHexAlPh_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Akimotoite: expected lattice and radiative thermal conductivities (k) in [W/m/K]
  // std::vector<double> Akimotoite_Expt_LatTCon = {8.59248, 1.78976, 2.38264, 2.76606, 2.37440};
  // std::vector<double> Akimotoite_Expt_RadTCon = {0.71149, 0.68558, 0.66935, 0.64541, 0.02321};
  // Preallocate a vector for storing thermal conductivities of minerals
  // std::vector<double> Expt_Minerals_LatTcond(MineralPar_Index, 0.0); // Lattice thermal conductivity
  // std::vector<double> Expt_Minerals_RadTcond(MineralPar_Index, 0.0); // Radiative thermal conductivity
  // std::vector<double> Expt_Minerals_TotTcond(MineralPar_Index, 0.0); // Total thermal conductivity

  // Perform element-wise sum
  for (size_t row = 0; row < temperatures.size(); ++row)
  {
      OlivineDry_Expt_TotTCon[row] = OlivineDry_Expt_LatTCon[row] + OlivineDry_Expt_RadTCon[row];
      // WadsleyDry_Expt_TotTCon[row] = WadsleyDry_Expt_LatTCon[row] + WadsleyDry_Expt_RadTCon[row];
      // RingwooDry_Expt_TotTCon[row] = RingwooDry_Expt_LatTCon[row] + RingwooDry_Expt_RadTCon[row];
      // En100Brigm_Expt_TotTCon[row] = En100Brigm_Expt_LatTCon[row] + En100Brigm_Expt_RadTCon[row];
      // En97Brigma_Expt_TotTCon[row] = En97Brigma_Expt_LatTCon[row] + En97Brigma_Expt_RadTCon[row];
      // En90Brigma_Expt_TotTCon[row] = En90Brigma_Expt_LatTCon[row] + En90Brigma_Expt_RadTCon[row];
      // AlMgBrigma_Expt_TotTCon[row] = AlMgBrigma_Expt_LatTCon[row] + AlMgBrigma_Expt_RadTCon[row];
      // FeAlBrigma_Expt_TotTCon[row] = FeAlBrigma_Expt_LatTCon[row] + FeAlBrigma_Expt_RadTCon[row];
      // OpxEnstati_Expt_TotTCon[row] = OpxEnstati_Expt_LatTCon[row] + OpxEnstati_Expt_RadTCon[row];
      // CpxDiopsid_Expt_TotTCon[row] = CpxDiopsid_Expt_LatTCon[row] + CpxDiopsid_Expt_RadTCon[row];
      // GrtPyropes_Expt_TotTCon[row] = GrtPyropes_Expt_LatTCon[row] + GrtPyropes_Expt_RadTCon[row];
      // GrtGrossul_Expt_TotTCon[row] = GrtGrossul_Expt_LatTCon[row] + GrtGrossul_Expt_RadTCon[row];
      // GrtAlmandi_Expt_TotTCon[row] = GrtAlmandi_Expt_LatTCon[row] + GrtAlmandi_Expt_RadTCon[row];
      // GrtMajorit_Expt_TotTCon[row] = GrtMajorit_Expt_LatTCon[row] + GrtMajorit_Expt_RadTCon[row];
      // QuartzPure_Expt_TotTCon[row] = QuartzPure_Expt_LatTCon[row] + QuartzPure_Expt_RadTCon[row];
      // CoesitSiO2_Expt_TotTCon[row] = CoesitSiO2_Expt_LatTCon[row] + CoesitSiO2_Expt_RadTCon[row];
      // Stishovite_Expt_TotTCon[row] = Stishovite_Expt_LatTCon[row] + Stishovite_Expt_RadTCon[row];
      // Al05Stisho_Expt_TotTCon[row] = Al05Stisho_Expt_LatTCon[row] + Al05Stisho_Expt_RadTCon[row];
      // Antigor010_Expt_TotTCon[row] = Antigor010_Expt_LatTCon[row] + Antigor010_Expt_RadTCon[row];
      // Antigor001_Expt_TotTCon[row] = Antigor001_Expt_LatTCon[row] + Antigor001_Expt_RadTCon[row];    
      // FeAlPhaseD_Expt_TotTCon[row] = FeAlPhaseD_Expt_LatTCon[row] + FeAlPhaseD_Expt_RadTCon[row];
      // Al02PhaseD_Expt_TotTCon[row] = Al02PhaseD_Expt_LatTCon[row] + Al02PhaseD_Expt_RadTCon[row];
      // Ferroper08_Expt_TotTCon[row] = Ferroper08_Expt_LatTCon[row] + Ferroper08_Expt_RadTCon[row];
      // Ferroper10_Expt_TotTCon[row] = Ferroper10_Expt_LatTCon[row] + Ferroper10_Expt_RadTCon[row];
      // Ferroper56_Expt_TotTCon[row] = Ferroper56_Expt_LatTCon[row] + Ferroper56_Expt_RadTCon[row];
      // Davemaoite_Expt_TotTCon[row] = Davemaoite_Expt_LatTCon[row] + Davemaoite_Expt_RadTCon[row];
      // NewHexAlPh_Expt_TotTCon[row] = NewHexAlPh_Expt_LatTCon[row] + NewHexAlPh_Expt_RadTCon[row];
      // Akimotoite_Expt_TotTCon[row] = Akimotoite_Expt_LatTCon[row] + Akimotoite_Expt_RadTCon[row];
  }

  // Assigning an array of values to in.composition (X) in [fraction]
   std::vector<std::vector<double>> compositions = {
      {1.00, 1.00, 1.00, 1.00, 1.00},
      {0.75, 0.75, 0.75, 0.75, 0.75},
      {0.50, 0.50, 0.50, 0.50, 0.50},
      {0.25, 0.25, 0.25, 0.25, 0.25}
  };
  
  // Initialize the expected value matrix with the same dimensions of the composition matrix
  std::vector<std::vector<double>> expected_total_Tcond(compositions.size(), std::vector<double>(compositions[0].size()));

  // Perform element-wise calculation
  for (size_t row = 0; row < compositions.size(); ++row)
  {
    for (size_t col = 0; col < compositions[row].size(); ++col)
    {
      expected_total_Tcond[row][col] = std::pow(OlivineDry_Expt_TotTCon[col], compositions[row][col]);
    }
  }

  std::vector<std::vector<double>> expected_conductivities = expected_total_Tcond;

  INFO("Checking thermal conductivity (k) for different temperatures (T), pressures (P) and compositions (X) values");

  // Loop over the different compositions
  for (size_t row = 0; row < expected_conductivities.size(); ++row)
  {
    in.composition[0] = compositions[row];  // Assign the current row of composition as model inputs
    model.evaluate(in, out);                // Call the function to compute the thermal conductivities

    // Loop over the different combinations of pressures (P) and temperatures (T)
    for (size_t i = 0; i < expected_conductivities[row].size(); ++i)
    {
      INFO("Computed OlivineDry k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed WadsleyDry k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed RingwooDry k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed En100Brigm k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed En97Brigma k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed En90Brigma k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed AlMgBrigma k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]"); 
      // INFO("Computed FeAlBrigma k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed OpxEnstati k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed CpxDiopsid k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed GrtPyropes k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed GrtGrossul k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed GrtAlmandi k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed GrtMajorit k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed QuartzPure k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed CoesitSiO2 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Stishovite k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Al05Stisho k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Antigor010 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Antigor001 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]"); 
      // INFO("Computed FeAlPhaseD k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Al02PhaseD k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Ferroper08 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Ferroper10 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Ferroper56 k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Davemaoite k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed NewHexAlPh k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      // INFO("Computed Akimotoite k at T= " << in.temperature[i] << "[K] ; P= " << in.pressure[i] << "[Pa] ; X= " << (in.composition[0][i])*100 << "[%] is " << out.thermal_conductivities[i] << "[W/m/K]");
      
      // Compare the computed thermal conductivity with the expected value
      REQUIRE(out.thermal_conductivities[i] == Approx(expected_conductivities[row][i]));
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

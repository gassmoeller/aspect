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

#include <aspect/material_model/thermal_conductivity/PT_dep_R_bounded.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      // Helper function: Compute lattice thermal conductivity
      double compute_lattice_thermal_conductivity(double a0, double b1, double ymin, double ymax, double P_log, double T_mod, double T_room, double n_exp)
      {
        // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
        double zSimpl = a0 + b1 * P_log;
        double ySimpl = std::exp(zSimpl);
        double yPrime = ySimpl / (1 + ySimpl);
        double yRealS = ymin + (ymax - ymin) * yPrime;
        double PDep_LatTCon = std::exp(yRealS);
        return PDep_LatTCon * std::pow((T_room / T_mod), n_exp);
      }
     
      // Helper function: Compute radiative thermal conductivity
      double compute_radiative_thermal_conductivity(double c0, double d1, double jmin, double jmax, double T_log)
      {
        // Compute the radiative thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
        double zSimpl = c0 + d1 * T_log;
        double ySimpl = std::exp(zSimpl);
        double yPrime = ySimpl / (1 + ySimpl);
        double yRealS = jmin + (jmax - jmin) * yPrime;
        return std::exp(yRealS);
      }
     
      // Helper function: Compute total thermal conductivity
      double compute_total_thermal_conductivity(double lattice_conductivity, double radiative_conductivity)
      {
        double thermal_conductivity = lattice_conductivity + radiative_conductivity;
        return thermal_conductivity;
       }
     
      // Helper function: Compute aggregate thermal conductivity
      // double compute_aggregate_thermal_conductivity(const std::vector<std::vector<double>> &thermal_conductivities, double min_frac, int col)
      // {
        // return std::pow(thermal_conductivities[3][col], min_frac);
      // }
     
      // Main function: 
      template <int dim>
      void
      PTdepRbounded<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        unsigned int MineralPar_Index = 0; // Initialize the counter

        // Define coefficients for lattice thermal conductivity for different minerals 

        // Coefficients for dry olivine 
        // retreived from fitting TDTR dataset of
        // [Chang et al., 2017, PNAS, vol 114, p. 4078-4081]
        // https://doi.org/10.1073/pnas.1616216114
        // mineral composition [Mg1.8 Fe0.2 SiO4]    
        unsigned int OlivineDry_Index = MineralPar_Index++;
        const double OlivineDry_LatTC_a0 =   -4.1241;
        const double OlivineDry_LatTC_b1 =    2.1469;
        const double OlivineDry_LatTC_ymin =  1.28093384543429;
        const double OlivineDry_LatTC_ymax =  2.60726381956037;
        const double OlivineDry_TDep_n_Exp =  0.5;     
        
        // Coefficients for dry wadsleyite
        // retreived from fitting dataset of
        // [Xu et al., 2004, PEPI, vol 143, pp. 321-336]
        // https://doi.org/10.1016/j.pepi.2004.03.005
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        unsigned int WadsleyDry_Index = MineralPar_Index++;
        const double WadsleyDry_LatTC_a0 =  -0.6656;
        const double WadsleyDry_LatTC_b1 =   0.3082;
        const double WadsleyDry_LatTC_ymin = 1.75735804249439;
        const double WadsleyDry_LatTC_ymax = 2.37090451537473; 
        const double WadsleyDry_TDep_n_Exp = 0.5; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Marzotto et al., 2020, GRL, vol 47, issue 13]
        // https://doi.org/10.1029/2020GL087607
        // mineral composition [(Mg1.79Fe0.17)Si1.02O4]
        unsigned int RingwooDry_Index = MineralPar_Index++;
        const double RingwooDry_LatTC_a0 =  -5.4624;
        const double RingwooDry_LatTC_b1 =   2.0791;
        const double RingwooDry_LatTC_ymin = 1.60943791241410;
        const double RingwooDry_LatTC_ymax = 2.94939766245070; 
        const double RingwooDry_TDep_n_Exp = 0.5; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [MgSiO3]
        unsigned int En100Brigm_Index = MineralPar_Index++;
        const double En100Brigm_LatTC_a0 =  -4.3687;
        const double En100Brigm_LatTC_b1 =   1.0766; 
        const double En100Brigm_LatTC_ymin = 2.376025820; 
        const double En100Brigm_LatTC_ymax = 5.010635294;  
        const double En100Brigm_TDep_n_Exp = 1.01000;

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.03Mg0.97SiO3]
        unsigned int En97Brigma_Index = MineralPar_Index++;
        const double En97Brigma_LatTC_a0 =  -4.520600000;
        const double En97Brigma_LatTC_b1 =   1.019900000; 
        const double En97Brigma_LatTC_ymin = 1.750524121; 
        const double En97Brigma_LatTC_ymax = 4.499809670;  
        const double En97Brigma_TDep_n_Exp = 0.56605;

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.1Mg0.9SiO3]
        unsigned int En90Brigma_Index = MineralPar_Index++;
        const double En90Brigma_LatTC_a0 =  -4.883100000;
        const double En90Brigma_LatTC_b1 =   0.980900000; 
        const double En90Brigma_LatTC_ymin = 1.333739493; 
        const double En90Brigma_LatTC_ymax = 4.382026635;  
        const double En90Brigma_TDep_n_Exp = 0.17054;

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Al,Mg)SiO3]
        unsigned int AlMgBrigma_Index = MineralPar_Index++;
        const double AlMgBrigma_LatTC_a0 =  -4.331500000;
        const double AlMgBrigma_LatTC_b1 =   1.027000000; 
        const double AlMgBrigma_LatTC_ymin = 1.845020046; 
        const double AlMgBrigma_LatTC_ymax = 4.605170186;  
        const double AlMgBrigma_TDep_n_Exp = 0.61983;

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Fe,Al,Mg)SiO3]
        unsigned int FeAlBrigma_Index = MineralPar_Index++;
        const double FeAlBrigma_LatTC_a0 =  -4.510600000;
        const double FeAlBrigma_LatTC_b1 =   1.066800000; 
        const double FeAlBrigma_LatTC_ymin = 1.389093953; 
        const double FeAlBrigma_LatTC_ymax = 3.912023005;  
        const double FeAlBrigma_TDep_n_Exp = 0.46815;
         
        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of 
        // [Schloessin & Dvorak, 1972, GJI, 27(5), 499-516]
        // https://doi.org/10.1111/j.1365-246X.1972.tb06105.x
        // mineral composition [Mg2Si2O6]
        unsigned int OpxEnstati_Index = MineralPar_Index++;
        const double OpxEnstati_LatTC_a0 =   -3.0047;
        const double OpxEnstati_LatTC_b1 =    2.6;
        const double OpxEnstati_LatTC_ymin =  1.760865151; 
        const double OpxEnstati_LatTC_ymax =  2.096937429;
        const double OpxEnstati_TDep_n_Exp =  0.5;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of 
        // [Wang et al., 2014, JGR: Solid Earth, 119(8), 6277-6287]
        // https://doi.org/10.1002/2014JB011208
        // mineral composition [CaMgSi2O6]
        unsigned int CpxDiopsid_Index = MineralPar_Index++;
        const double CpxDiopsid_LatTC_a0 =   -3.251100000;
        const double CpxDiopsid_LatTC_b1 =    1.689100000;
        const double CpxDiopsid_LatTC_ymin =  1.793640135; 
        const double CpxDiopsid_LatTC_ymax =  2.389462023;
        const double CpxDiopsid_TDep_n_Exp =  0.5;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [Mg3Al2Si3O12]
        unsigned int GrtPyropes_Index = MineralPar_Index++;
        const double GrtPyropes_LatTC_a0 =   -4.3637;
        const double GrtPyropes_LatTC_b1 =    2.0368;
        const double GrtPyropes_LatTC_ymin =  1.481604541; 
        const double GrtPyropes_LatTC_ymax =  2.443131606;
        const double GrtPyropes_TDep_n_Exp =  0.4314;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Ca0.986Fe0.014)3Al2(SiO4)3]
        unsigned int GrtGrossul_Index = MineralPar_Index++;
        const double GrtGrossul_LatTC_a0 =  -4.7584;
        const double GrtGrossul_LatTC_b1 =   2.0816;
        const double GrtGrossul_LatTC_ymin = 1.410986974; 
        const double GrtGrossul_LatTC_ymax = 2.457625992;
        const double GrtGrossul_TDep_n_Exp = 0.4589;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3]
        unsigned int GrtAlmandi_Index = MineralPar_Index++;
        const double GrtAlmandi_LatTC_a0 =  -4.5047;
        const double GrtAlmandi_LatTC_b1 =   2.0988;
        const double GrtAlmandi_LatTC_ymin = 1.223775432; 
        const double GrtAlmandi_LatTC_ymax = 2.374762159;
        const double GrtAlmandi_TDep_n_Exp = 0.4172;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Giesting et al.2004  EPSL, 218(1-2), 45-56]
        // https://doi.org/10.1016/S0012-821X(03)00630-7
        // mineral composition [Mg3(MgSi)(SiO4)3]
        unsigned int GrtMajorit_Index = MineralPar_Index++;
        const double GrtMajorit_LatTC_a0 =  -4.3637;
        const double GrtMajorit_LatTC_b1 =   2.0368;
        const double GrtMajorit_LatTC_ymin = 2.279316466; 
        const double GrtMajorit_LatTC_ymax = 2.718047842;
        const double GrtMajorit_TDep_n_Exp = 0.5;

        // Coefficients for quartz 
        // retreived from fitting dataset of
        // [Xiong et al., 2019 - Journal of Applied Physics, 126(21)]
        // https://doi.org/10.1063/1.5114992
        // mineral composition [SiO2]
        unsigned int QuartzPure_Index = MineralPar_Index++;
        const double QuartzPure_LatTC_a0 =   -2.0203;
        const double QuartzPure_LatTC_b1 =    2.4456;
        const double QuartzPure_LatTC_ymin =  2.260981081; 
        const double QuartzPure_LatTC_ymax =  2.745391462;
        const double QuartzPure_TDep_n_Exp =  1.015433333;

        // Coefficients for coesite
        // retreived from fitting dataset of
        // [Yukutake & Shimada, 1978, PEPI, 17(3), 193-200]
        // https://doi.org/10.1016/0031-9201(78)90036-5
        // mineral composition [SiO2]
        unsigned int CoesitSiO2_Index = MineralPar_Index++;
        const double CoesitSiO2_LatTC_a0 =   -12.728;
        const double CoesitSiO2_LatTC_b1 =    2.9998;
        const double CoesitSiO2_LatTC_ymin =  1.982022416; 
        const double CoesitSiO2_LatTC_ymax =  2.249036030;
        const double CoesitSiO2_TDep_n_Exp =  1.015433333;

        // Coefficients for stishovite
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584, 117477]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [SiO2]
        unsigned int Stishovite_Index = MineralPar_Index++;
        // < 52 [GPa]
        const double Stishovite_LatTC_a0 =  16.917;
        const double Stishovite_LatTC_b1 = -4.6187;
        const double Stishovite_LatTC_ymin = 4.096113064; 
        const double Stishovite_LatTC_ymax = 4.217974805;
        // 52-56 [GPa]
        // const double Stishovite_LatTC_c0 = -156.12;
        // const double Stishovite_LatTC_d1 =  39.182;
        // const double Stishovite_LatTC_ymin = 4.077505176;
        // const double Stishovite_LatTC_ymax = 4.264199335;
        // > 56 [GPa]
        // const double Stishovite_LatTC_e0 = -12.728;
        // const double Stishovite_LatTC_f1 =  2.9998;
        // const double Stishovite_LatTC_ymin = 3.960844211;
        // const double Stishovite_LatTC_ymax = 4.738489125;
        // Temperature
        const double Stishovite_TDep_n_Exp = 0.5;

        // Coefficients for Al-stishovite (5 vol%)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584, 117477]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [(Al,Si)O2]
        unsigned int Al05Stisho_Index = MineralPar_Index++;
        const double Al05Stisho_LatTC_a0 = -6.4411;
        const double Al05Stisho_LatTC_b1 =  1.5885;
        const double Al05Stisho_LatTC_ymin = 3.188855035;
        const double Al05Stisho_LatTC_ymax = 4.154336189;
        const double Al05Stisho_TDep_n_Exp = 0.5;

        // Coefficients for antigorite (serpentine)
        // retreived from fitting dataset of
        // [Chien et al., 2024, Nature Communications, 15(1), 5198.]
        // https://doi.org/10.1038/s41467-024-49418-3
        // mineral composition [(Mg2.80Fe0.05)Si2.08O5(OH)3.77]
        // 010 direction
        unsigned int Antigor010_Index = MineralPar_Index++;
        const double Antigor010_LatTC_a0 = -4.3374;
        const double Antigor010_LatTC_b1 =  2.0217;
        const double Antigor010_LatTC_ymin = 1.519513205;
        const double Antigor010_LatTC_ymax = 2.434491480;
        const double Antigor010_TDep_n_Exp = 0.5;
        // 001 direction
        unsigned int Antigor001_Index = MineralPar_Index++;
        const double Antigor001_LatTC_a0 = -3.1109;
        const double Antigor001_LatTC_b1 =  2.0644;
        const double Antigor001_LatTC_ymin = 0.067658648;
        const double Antigor001_LatTC_ymax = 1.552797578;
        const double Antigor001_TDep_n_Exp = 0.5;

        // Coefficients for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, JGR: Solid Earth, vol. 127(6), e2022JB024556]
        // https://doi.org/10.1029/2022JB024556
        // mineral composition [Mg1.19Fe0.12Al0.174Si1.71H2.02O6]
        unsigned int FeAlPhaseD_Index = MineralPar_Index++;
        // (Fe,Al)-Phase D - 0-24 [GPa]
        const double FeAlPhaseD_LatTC_a0 = -3.9909;
        const double FeAlPhaseD_LatTC_b1 =  1.7710;
        const double FeAlPhaseD_LatTC_ymin = 0.956005323; 
        const double FeAlPhaseD_LatTC_ymax = 1.747361025;
        // (Fe,Al)-Phase D - 24-38 [GPa]
        // const double FeAlPhaseD_LatTC_c0 = -32.890;
        // const double FeAlPhaseD_LatTC_d1 =  9.6282;
        // const double FeAlPhaseD_LatTC_ymin = 1.442700096; 
        // const double FeAlPhaseD_LatTC_ymax = 3.512683919;
        // (Fe,Al)-Phase D - 38-48 [GPa]
        // const double FeAlPhaseD_LatTC_e0 =  141.88;
        // const double FeAlPhaseD_LatTC_f1 = -37.409;
        // const double FeAlPhaseD_LatTC_ymin = 1.789742436; 
        // const double FeAlPhaseD_LatTC_ymax = 3.270312073;
        // (Fe,Al)-Phase D - > 48 [GPa]
        // const double FeAlPhaseD_LatTC_g0 = -23.986;
        // const double FeAlPhaseD_LatTC_h1 =  6.1139;
        // const double FeAlPhaseD_LatTC_ymin = 1.313988596; 
        // const double FeAlPhaseD_LatTC_ymax = 2.992561000;
        // Temperature-dependency
        const double FeAlPhaseD_TDep_n_Exp = 0.5;

        // Coefficients for Al-phase D (Dense Hydrous Magnesium Silicate)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, JGR: Solid Earth, vol. 127(6), e2022JB024556]
        // https://doi.org/10.1029/2022JB024556
        // mineral composition [Mg1.29Al0.17Si1.73H1.98O6]
        unsigned int Al02PhaseD_Index = MineralPar_Index++;
        const double Al02PhaseD_LatTC_a0 = -6.1829;
        const double Al02PhaseD_LatTC_b1 =  1.8514;
        const double Al02PhaseD_LatTC_ymin = 1.285874399; 
        const double Al02PhaseD_LatTC_ymax = 3.502412041;
        const double Al02PhaseD_TDep_n_Exp = 0.5;

        // Coefficients for ferropericlase (Mg1-xFexO)
        // retreived from fitting dataset of
        // [Hsieh et al., 2018, PNAS, vol. 115, no. 16, p. 4099-4104]
        // mineral composition [Mg0.92Fe0.08O] - (8% Iron)
        unsigned int Ferroper08_Index = MineralPar_Index++;
        const double Ferroper08_LatTC_a0 = -6.9942;
        const double Ferroper08_LatTC_b1 =  1.953;
        const double Ferroper08_LatTC_ymin = 1.629240539; 
        const double Ferroper08_LatTC_ymax = 4.118362306;
        const double Ferroper08_TDep_n_Exp = 0.5;
        // mineral composition [Mg0.90Fe0.10O] - (10% Iron)
        unsigned int Ferroper10_Index = MineralPar_Index++;
        const double Ferroper10_LatTC_a0 = -5.2408;
        const double Ferroper10_LatTC_b1 =  0.9649;
        const double Ferroper10_LatTC_ymin = 1.2490430868; 
        const double Ferroper10_LatTC_ymax = 3.9318256327;
        const double Ferroper10_TDep_n_Exp = 0.025;
        // mineral composition [Mg0.44Fe0.56O] (56% Iron)
        unsigned int Ferroper56_Index = MineralPar_Index++;
        const double Ferroper56_LatTC_a0 = -3.8298;
        const double Ferroper56_LatTC_b1 =  1.1507;
        const double Ferroper56_LatTC_ymin = 0.993251773; 
        const double Ferroper56_LatTC_ymax = 3.592193222;
        const double Ferroper56_TDep_n_Exp = 0.5;

        // Coefficients for davemaoite 
        // retreived from fitting dataset of
        // [Zhang et al., 2021, Physical Review B, vol. 104, 184101]
        // https://doi.org/10.1103/PhysRevB.104.184101
        // mineral composition [CaSiO3]
        unsigned int Davemaoite_Index = MineralPar_Index++;
        const double Davemaoite_LatTC_a0 = -4.7377;
        const double Davemaoite_LatTC_b1 =  1.3661;
        const double Davemaoite_LatTC_ymin = 2.388762789; 
        const double Davemaoite_LatTC_ymax = 4.045106030;
        const double Davemaoite_TDep_n_Exp = 0.5;

        // Coefficients for new-hexagonal-alluminium-phase (FeNAL) 
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [Na0.71Mg2.05Al4.62Si1.16Fe(2+)0.09Fe(3+)0.17O12]
        unsigned int NewHexAlPh_Index = MineralPar_Index++;
        const double NewHexAlPh_LatTC_a0 = -29.421;
        const double NewHexAlPh_LatTC_b1 =  7.7792;
        const double NewHexAlPh_LatTC_ymin = 2.363551955; 
        const double NewHexAlPh_LatTC_ymax = 3.653998874;
        const double NewHexAlPh_TDep_n_Exp = 0.5;

        // Coefficients for Akimotoite
        // assumed to be equal to En100-Bridgmanite
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [MgSiO3]
        unsigned int Akimotoite_Index = MineralPar_Index++;
        const double Akimotoite_LatTC_a0 =  -4.368700000;
        const double Akimotoite_LatTC_b1 =   1.076600000; 
        const double Akimotoite_LatTC_ymin = 2.376025820; 
        const double Akimotoite_LatTC_ymax = 5.010635294;  
        const double Akimotoite_TDep_n_Exp = 0.5;

        // Define coefficients for radiative thermal conductivity of different minerals

        // Coefficients for dry olivine
        // retreived from fitting dataset of
        // [Marzotto et al. 2025, Nature Communication, under review]
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double OlivineDry_RadTC_c0 =   -10.00900000;
        const double OlivineDry_RadTC_d1 =    1.883900000;
        const double OlivineDry_RadTC_jmin = -23.02585093;
        const double OlivineDry_RadTC_jmax =  1.289885976;

        // Coefficients for dry wadsleyite
        // retreived from fitting dataset of
        // [Thomas et al., 2012, EPSL, vol. 357, p. 130-136.]
        // https://doi.org/10.1016/j.epsl.2012.09.035
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double WadsleyDry_RadTC_c0 =   -21.717;
        const double WadsleyDry_RadTC_d1 =    3.4271;
        const double WadsleyDry_RadTC_jmin = -23.025850930;
        const double WadsleyDry_RadTC_jmax =  1.0773006810; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Thomas et al., 2012, EPSL, vol. 357, p. 130-136.]
        // https://doi.org/10.1016/j.epsl.2012.09.035
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double RingwooDry_RadTC_c0 =   -23.067000000;
        const double RingwooDry_RadTC_d1 =    3.5985000000;
        const double RingwooDry_RadTC_jmin = -23.025850930;
        const double RingwooDry_RadTC_jmax =  0.4027750000; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [MgSiO3]
        const double En100Brigm_RadTC_c0 =  66.278;
        const double En100Brigm_RadTC_d1 = -8.2756; 
        const double En100Brigm_RadTC_jmin = -7.2568958208;         
        const double En100Brigm_RadTC_jmax = -0.3403920329;         

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Fe0.03Mg0.97SiO3]
        const double En97Brigma_RadTC_c0 =  66.278;
        const double En97Brigma_RadTC_d1 = -8.2756; 
        const double En97Brigma_RadTC_jmin = -7.2568958208;         
        const double En97Brigma_RadTC_jmax = -0.3403920329;  

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Fe0.1Mg0.9SiO3]
        const double En90Brigma_RadTC_c0 =  66.278;
        const double En90Brigma_RadTC_d1 = -8.2756; 
        const double En90Brigma_RadTC_jmin = -7.2568958208;         
        const double En90Brigma_RadTC_jmax = -0.3403920329;  

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [(Al,Mg)SiO3]
        const double AlMgBrigma_RadTC_c0 =  66.278;
        const double AlMgBrigma_RadTC_d1 = -8.2756; 
        const double AlMgBrigma_RadTC_jmin = -7.2568958208;         
        const double AlMgBrigma_RadTC_jmax = -0.3403920329;  

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [(Fe,Al,Mg)SiO3]
        const double FeAlBrigma_RadTC_c0 =  66.278;
        const double FeAlBrigma_RadTC_d1 = -8.2756; 
        const double FeAlBrigma_RadTC_jmin = -7.2568958208;         
        const double FeAlBrigma_RadTC_jmax = -0.3403920329;  

        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg2Si2O6]
        const double OpxEnstati_RadTC_c0 =   -13.532000000;
        const double OpxEnstati_RadTC_d1 =    2.4004000000;
        const double OpxEnstati_RadTC_jmin = -23.025850930; 
        const double OpxEnstati_RadTC_jmax =  1.4456685920;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [CaMgSi2O6]
        const double CpxDiopsid_RadTC_c0 =   -14.286000000;       
        const double CpxDiopsid_RadTC_d1 =    2.5119000000; 
        const double CpxDiopsid_RadTC_jmin = -23.025850930; 
        const double CpxDiopsid_RadTC_jmax =  1.4434214740;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3Al2(SiO4)3]
        const double GrtPyropes_RadTC_c0 =   -11.782000000;
        const double GrtPyropes_RadTC_d1 =    2.0718000000;
        const double GrtPyropes_RadTC_jmin = -23.025850930; 
        const double GrtPyropes_RadTC_jmax =  1.4479836950;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Ca3Al2(SiO4)3]
        const double GrtGrossul_RadTC_c0 =   -11.782000000;
        const double GrtGrossul_RadTC_d1 =    2.0718000000;
        const double GrtGrossul_RadTC_jmin = -23.025850930; 
        const double GrtGrossul_RadTC_jmax =  1.4479836950;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3 ]
        const double GrtAlmandi_RadTC_c0  =  -11.782000000;
        const double GrtAlmandi_RadTC_d1 =    2.0718000000;
        const double GrtAlmandi_RadTC_jmin = -23.025850930; 
        const double GrtAlmandi_RadTC_jmax =  1.4479836950;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3(MgSi)(SiO4)3]
        const double GrtMajorit_RadTC_c0  =  -11.782000000;
        const double GrtMajorit_RadTC_d1 =    2.0718000000;
        const double GrtMajorit_RadTC_jmin = -23.025850930;
        const double GrtMajorit_RadTC_jmax =  1.4479836950;

        // Coefficients for quartz
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double QuartzPure_RadTC_c0  =   0;
        const double QuartzPure_RadTC_d1  =   0;
        const double QuartzPure_RadTC_jmin = -23.025850930; 
        const double QuartzPure_RadTC_jmax = -23.050000000;

        // Coefficients for coesite
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double CoesitSiO2_RadTC_c0  =   0;
        const double CoesitSiO2_RadTC_d1  =   0;
        const double CoesitSiO2_RadTC_jmin = -23.025850930; 
        const double CoesitSiO2_RadTC_jmax = -23.050000000;

        // Coefficients for stishovite
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double Stishovite_RadTC_c0  =   0;
        const double Stishovite_RadTC_d1  =   0;
        const double Stishovite_RadTC_jmin = -23.025850930; 
        const double Stishovite_RadTC_jmax = -23.050000000;

        // Coefficients for Al-stishovite (5 vol%)
        // assumed 0 for now - no data available
        // mineral composition [(Al,Si)O2]
        const double Al05Stisho_RadTC_c0  =   0;
        const double Al05Stisho_RadTC_d1  =   0;
        const double Al05Stisho_RadTC_jmin = -23.025850930; 
        const double Al05Stisho_RadTC_jmax = -23.050000000;

        // Coefficients for antigorite (serpentine)
        // assumed 0 for now - no data available
        // mineral composition [///]
        // 010 direction
        const double Antigor010_RadTC_c0  =   0;
        const double Antigor010_RadTC_d1  =   0;
        const double Antigor010_RadTC_jmin = -23.025850930; 
        const double Antigor010_RadTC_jmax = -23.050000000;
        // 001 direction
        const double Antigor001_RadTC_c0  =   0;
        const double Antigor001_RadTC_d1  =   0;
        const double Antigor001_RadTC_jmin = -23.025850930; 
        const double Antigor001_RadTC_jmax = -23.050000000;

        // Coefficients for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
        // assumed 0 for now - no data available
        // mineral composition [Mg1.19Fe0.12Al0.174Si1.71H2.02O6]
        const double FeAlPhaseD_RadTC_c0  =   0;
        const double FeAlPhaseD_RadTC_d1  =   0;
        const double FeAlPhaseD_RadTC_jmin = -23.025850930; 
        const double FeAlPhaseD_RadTC_jmax = -23.050000000;

        // Coefficients for Al-phase D (Dense Hydrous Magnesium Silicate)
        // assumed 0 for now - no data available
        // mineral composition [Mg1.29Al0.17Si1.73H1.98O6]
        const double Al02PhaseD_RadTC_c0  =   0;
        const double Al02PhaseD_RadTC_d1  =   0;
        const double Al02PhaseD_RadTC_jmin = -23.025850930; 
        const double Al02PhaseD_RadTC_jmax = -23.050000000;

        // Coefficients for ferropericlase (Mg1-xFexO)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Mg0.92Fe0.08O] - (8% Iron)
        const double Ferroper08_RadTC_c0 =  66.278;
        const double Ferroper08_RadTC_d1 = -8.2756; 
        const double Ferroper08_RadTC_jmin = -7.2568958208;         
        const double Ferroper08_RadTC_jmax = -0.3403920329;  
        // mineral composition [Mg0.90Fe0.10O] - (10% Iron)
        const double Ferroper10_RadTC_c0 =  66.278;
        const double Ferroper10_RadTC_d1 = -8.2756; 
        const double Ferroper10_RadTC_jmin = -7.2568958208;         
        const double Ferroper10_RadTC_jmax = -0.3403920329;  
        // mineral composition [Mg0.44Fe0.56O] (56% Iron)
        const double Ferroper56_RadTC_c0 =  66.278;
        const double Ferroper56_RadTC_d1 = -8.2756; 
        const double Ferroper56_RadTC_jmin = -7.2568958208;         
        const double Ferroper56_RadTC_jmax = -0.3403920329;  

        // Coefficients for davemaoite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [CaSiO3] 
        const double Davemaoite_RadTC_c0 =  66.278;
        const double Davemaoite_RadTC_d1 = -8.2756; 
        const double Davemaoite_RadTC_jmin = -7.2568958208;         
        const double Davemaoite_RadTC_jmax = -0.3403920329;  

        // Coefficients for new-hexagonal-alluminium-phase (FeNAL) 
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Na0.71Mg2.05Al4.62Si1.16Fe(2+)0.09Fe(3+)0.17O12] 
        const double NewHexAlPh_RadTC_c0 =  66.278;
        const double NewHexAlPh_RadTC_d1 = -8.2756; 
        const double NewHexAlPh_RadTC_jmin = -7.2568958208;         
        const double NewHexAlPh_RadTC_jmax = -0.3403920329;  

        // Coefficients for akimotoite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [MgSiO3] 
        const double Akimotoite_RadTC_c0 =  66.278;
        const double Akimotoite_RadTC_d1 = -8.2756; 
        const double Akimotoite_RadTC_jmin = -7.2568958208;         
        const double Akimotoite_RadTC_jmax = -0.3403920329;  

        // Preallocate a vector for mineral fractions of different rocks
        // Pyrolite (58% olivine, 13% pyrope, 18% ensatite, 11% diopside)
        std::vector<double> MinFract_Pyrolite = {0.58,0.13,0.18,0.11}; 
        // Harzburgite (80% olivine, 20% ensatite)
        std::vector<double> MinFract_Harzburg = {0.80,0.20}; 
        // Meta-basaltic crust MORB (80% diopside, 20% pyrope)
        std::vector<double> MinFract_MetaMORB = {0.80,0.20}; 
        // Dunite (100% olivine)
        std::vector<double> MinFract_DuniteOl = {1.00}; 

        #include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

        // Check if the sum of Rock Mineral Fraction is equal to 1
        double sum_min_fract_pyrolite = std::accumulate(MinFract_Pyrolite.begin(), MinFract_Pyrolite.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_pyrolite - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of MinFract_Pyrolite must be equal to 1."));
        double sum_min_fract_harzburg = std::accumulate(MinFract_Harzburg.begin(), MinFract_Harzburg.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_harzburg - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of MinFract_Harzburg must be equal to 1."));
        double sum_min_fract_metaMORB = std::accumulate(MinFract_MetaMORB.begin(), MinFract_MetaMORB.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_metaMORB - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of MinFract_MetaMORB must be equal to 1."));
        double sum_min_fract_DuniteOl = std::accumulate(MinFract_DuniteOl.begin(), MinFract_DuniteOl.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_DuniteOl - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of MinFract_DuniteOl must be equal to 1."));

        // Define room temperature [K] 
        const double T_room = 298.15; 

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {

          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> All_Minerals_LatTcond(MineralPar_Index, 0.0); // Lattice thermal conductivity
          std::vector<double> All_Minerals_RadTcond(MineralPar_Index, 0.0); // Radiative thermal conductivity
          std::vector<double> All_Minerals_TotTcond(MineralPar_Index, 0.0); // Total thermal conductivity
          // Preallocate a matrix for storing thermal conductivities of minerals
          std::vector<std::vector<double>> All_Minerals_TConds(MineralPar_Index, std::vector<double>(3, 0.0));

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;

          // Compute natural logarithm of pressure and temperature 
          double P_log = std::log(P_GPa);
          double T_log = std::log(in.temperature[i]);

          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];

          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];

          // Compute lattice thermal conductivities for DryOlivine
          double OlivineDry_LatTCon = compute_lattice_thermal_conductivity(
            OlivineDry_LatTC_a0, OlivineDry_LatTC_b1, OlivineDry_LatTC_ymin, OlivineDry_LatTC_ymax,
            P_log, T_mod, T_room, OlivineDry_TDep_n_Exp);
          All_Minerals_LatTcond[OlivineDry_Index] = OlivineDry_LatTCon;
          // Compute lattice thermal conductivities for Dry Wadsleyite 
          double WadsleyDry_LatTCon = compute_lattice_thermal_conductivity(
            WadsleyDry_LatTC_a0, WadsleyDry_LatTC_b1, WadsleyDry_LatTC_ymin, WadsleyDry_LatTC_ymax,
            P_log, T_mod, T_room, WadsleyDry_TDep_n_Exp);
          All_Minerals_LatTcond[WadsleyDry_Index] = WadsleyDry_LatTCon;
          // Compute lattice thermal conductivities for Dry Ringwoodite  
          double RingwooDry_LatTCon = compute_lattice_thermal_conductivity(
            RingwooDry_LatTC_a0, RingwooDry_LatTC_b1, RingwooDry_LatTC_ymin, RingwooDry_LatTC_ymax,
            P_log, T_mod, T_room, RingwooDry_TDep_n_Exp);
          All_Minerals_LatTcond[RingwooDry_Index] = RingwooDry_LatTCon;
          // Compute lattice thermal conductivities for Mg-Bridgmanite
          double En100Brigm_LatTCon = compute_lattice_thermal_conductivity(
            En100Brigm_LatTC_a0, En100Brigm_LatTC_b1, En100Brigm_LatTC_ymin, En100Brigm_LatTC_ymax,
            P_log, T_mod, T_room, En100Brigm_TDep_n_Exp);
          All_Minerals_LatTcond[En100Brigm_Index] = En100Brigm_LatTCon;
          // Compute lattice thermal conductivities for Fe-Bridgmanite (3%)
          double En97Brigma_LatTCon = compute_lattice_thermal_conductivity(
            En97Brigma_LatTC_a0, En97Brigma_LatTC_b1, En97Brigma_LatTC_ymin, En97Brigma_LatTC_ymax,
            P_log, T_mod, T_room, En97Brigma_TDep_n_Exp);
          All_Minerals_LatTcond[En97Brigma_Index] = En97Brigma_LatTCon;
          // Compute lattice thermal conductivities for Fe-Bridgmanite (10%)
          double En90Brigma_LatTCon = compute_lattice_thermal_conductivity(
            En90Brigma_LatTC_a0, En90Brigma_LatTC_b1, En90Brigma_LatTC_ymin, En90Brigma_LatTC_ymax,
            P_log, T_mod, T_room, En90Brigma_TDep_n_Exp);
          All_Minerals_LatTcond[En90Brigma_Index] = En90Brigma_LatTCon;
          // Compute lattice thermal conductivities for Al-Bridgmanite
          double AlMgBrigma_LatTCon = compute_lattice_thermal_conductivity(
            AlMgBrigma_LatTC_a0, AlMgBrigma_LatTC_b1, AlMgBrigma_LatTC_ymin, AlMgBrigma_LatTC_ymax,
            P_log, T_mod, T_room, AlMgBrigma_TDep_n_Exp);
          All_Minerals_LatTcond[AlMgBrigma_Index] = AlMgBrigma_LatTCon;
          // Compute lattice thermal conductivities for Fe,Al-Bridgmanite
          double FeAlBrigma_LatTCon = compute_lattice_thermal_conductivity(
            FeAlBrigma_LatTC_a0, FeAlBrigma_LatTC_b1, FeAlBrigma_LatTC_ymin, FeAlBrigma_LatTC_ymax,
            P_log, T_mod, T_room, FeAlBrigma_TDep_n_Exp);
          All_Minerals_LatTcond[FeAlBrigma_Index] = FeAlBrigma_LatTCon;
          // Compute lattice thermal conductivities for Orthopyroxene (Enstatite)
          double OpxEnstati_LatTCon = compute_lattice_thermal_conductivity(
            OpxEnstati_LatTC_a0, OpxEnstati_LatTC_b1, OpxEnstati_LatTC_ymin, OpxEnstati_LatTC_ymax,
            P_log, T_mod, T_room, OpxEnstati_TDep_n_Exp);
          All_Minerals_LatTcond[OpxEnstati_Index] = OpxEnstati_LatTCon;
          // Compute lattice thermal conductivities for Clinopyroxene (Diopside)
          double CpxDiopsid_LatTCon = compute_lattice_thermal_conductivity(
             CpxDiopsid_LatTC_a0, CpxDiopsid_LatTC_b1, CpxDiopsid_LatTC_ymin, CpxDiopsid_LatTC_ymax,
             P_log, T_mod, T_room, CpxDiopsid_TDep_n_Exp);
          All_Minerals_LatTcond[CpxDiopsid_Index] = CpxDiopsid_LatTCon;
          // Compute lattice thermal conductivities for Garnet (Pyrope)
          double GrtPyropes_LatTCon = compute_lattice_thermal_conductivity(
            GrtPyropes_LatTC_a0, GrtPyropes_LatTC_b1, GrtPyropes_LatTC_ymin, GrtPyropes_LatTC_ymax,
            P_log, T_mod, T_room, GrtPyropes_TDep_n_Exp);
          All_Minerals_LatTcond[GrtPyropes_Index] = GrtPyropes_LatTCon;
          // Compute lattice thermal conductivities for Garnet (Grossular)
          double GrtGrossul_LatTCon = compute_lattice_thermal_conductivity(
            GrtGrossul_LatTC_a0, GrtGrossul_LatTC_b1, GrtGrossul_LatTC_ymin, GrtGrossul_LatTC_ymax,
            P_log, T_mod, T_room, GrtGrossul_TDep_n_Exp);
          All_Minerals_LatTcond[GrtGrossul_Index] = GrtGrossul_LatTCon;
          // Compute lattice thermal conductivities for Garnet (Almandine)
          double GrtAlmandi_LatTCon = compute_lattice_thermal_conductivity(
            GrtAlmandi_LatTC_a0, GrtAlmandi_LatTC_b1, GrtAlmandi_LatTC_ymin, GrtAlmandi_LatTC_ymax,
            P_log, T_mod, T_room, GrtAlmandi_TDep_n_Exp);
          All_Minerals_LatTcond[GrtAlmandi_Index] = GrtAlmandi_LatTCon;
          // Compute lattice thermal conductivities for Garnet (Majorite)
          double GrtMajorit_LatTCon = compute_lattice_thermal_conductivity(
            GrtMajorit_LatTC_a0, GrtMajorit_LatTC_b1, GrtMajorit_LatTC_ymin, GrtMajorit_LatTC_ymax,
            P_log, T_mod, T_room, GrtMajorit_TDep_n_Exp);
          All_Minerals_LatTcond[GrtMajorit_Index] = GrtMajorit_LatTCon;
          // Compute lattice thermal conductivities for Quartz
          double QuartzPure_LatTCon = compute_lattice_thermal_conductivity(
            QuartzPure_LatTC_a0, QuartzPure_LatTC_b1, QuartzPure_LatTC_ymin, QuartzPure_LatTC_ymax,
            P_log, T_mod, T_room, QuartzPure_TDep_n_Exp);
          All_Minerals_LatTcond[QuartzPure_Index] = QuartzPure_LatTCon;
          // Compute lattice thermal conductivities for Coesite
          double CoesitSiO2_LatTCon = compute_lattice_thermal_conductivity(
            CoesitSiO2_LatTC_a0, CoesitSiO2_LatTC_b1, CoesitSiO2_LatTC_ymin, CoesitSiO2_LatTC_ymax,
            P_log, T_mod, T_room, CoesitSiO2_TDep_n_Exp);
          All_Minerals_LatTcond[CoesitSiO2_Index] = CoesitSiO2_LatTCon;
          // Compute lattice thermal conductivities for Stishovite
          double Stishovite_LatTCon = compute_lattice_thermal_conductivity(
            Stishovite_LatTC_a0, Stishovite_LatTC_b1, Stishovite_LatTC_ymin, Stishovite_LatTC_ymax,
            P_log, T_mod, T_room, Stishovite_TDep_n_Exp);
          All_Minerals_LatTcond[Stishovite_Index] = Stishovite_LatTCon;
          // Compute lattice thermal conductivities for Al-stishovite (5 vol%)
          double Al05Stisho_LatTCon = compute_lattice_thermal_conductivity(
            Al05Stisho_LatTC_a0, Al05Stisho_LatTC_b1, Al05Stisho_LatTC_ymin, Al05Stisho_LatTC_ymax,
            P_log, T_mod, T_room, Al05Stisho_TDep_n_Exp);
          All_Minerals_LatTcond[Al05Stisho_Index] = Al05Stisho_LatTCon;
          // Compute lattice thermal conductivities for Antigorite (010)
          double Antigor010_LatTCon = compute_lattice_thermal_conductivity(
            Antigor010_LatTC_a0, Antigor010_LatTC_b1, Antigor010_LatTC_ymin, Antigor010_LatTC_ymax,
            P_log, T_mod, T_room, Antigor010_TDep_n_Exp);
          All_Minerals_LatTcond[Antigor010_Index] = Antigor010_LatTCon;
          // Compute lattice thermal conductivities for Antigorite (001)
          double Antigor001_LatTCon = compute_lattice_thermal_conductivity(
            Antigor001_LatTC_a0, Antigor001_LatTC_b1, Antigor001_LatTC_ymin, Antigor001_LatTC_ymax,
            P_log, T_mod, T_room, Antigor001_TDep_n_Exp);
          All_Minerals_LatTcond[Antigor001_Index] = Antigor001_LatTCon;
          // Compute lattice thermal conductivities for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
          double FeAlPhaseD_LatTCon = compute_lattice_thermal_conductivity(
            FeAlPhaseD_LatTC_a0, FeAlPhaseD_LatTC_b1, FeAlPhaseD_LatTC_ymin, FeAlPhaseD_LatTC_ymax,
            P_log, T_mod, T_room, FeAlPhaseD_TDep_n_Exp);
          All_Minerals_LatTcond[FeAlPhaseD_Index] = FeAlPhaseD_LatTCon;
          // Compute lattice thermal conductivities for Al-phase D (Dense Hydrous Magnesium Silicate)
          double Al02PhaseD_LatTCon = compute_lattice_thermal_conductivity(
            Al02PhaseD_LatTC_a0, Al02PhaseD_LatTC_b1, Al02PhaseD_LatTC_ymin, Al02PhaseD_LatTC_ymax,
            P_log, T_mod, T_room, Al02PhaseD_TDep_n_Exp);
          All_Minerals_LatTcond[Al02PhaseD_Index] = Al02PhaseD_LatTCon;
          // Compute lattice thermal conductivities for Ferropericlase (Mg92Fe8O)
          double Ferroper08_LatTCon = compute_lattice_thermal_conductivity(
            Ferroper08_LatTC_a0, Ferroper08_LatTC_b1, Ferroper08_LatTC_ymin, Ferroper08_LatTC_ymax,
            P_log, T_mod, T_room, Ferroper08_TDep_n_Exp);
          All_Minerals_LatTcond[Ferroper08_Index] = Ferroper08_LatTCon;
          // Compute lattice thermal conductivities for Ferropericlase (Mg90Fe10O)
          double Ferroper10_LatTCon = compute_lattice_thermal_conductivity(
            Ferroper10_LatTC_a0, Ferroper10_LatTC_b1, Ferroper10_LatTC_ymin, Ferroper10_LatTC_ymax,
            P_log, T_mod, T_room, Ferroper10_TDep_n_Exp);
          All_Minerals_LatTcond[Ferroper10_Index] = Ferroper10_LatTCon;
          // Compute lattice thermal conductivities for Ferropericlase (Mg56Fe44O)
          double Ferroper56_LatTCon = compute_lattice_thermal_conductivity(
            Ferroper56_LatTC_a0, Ferroper56_LatTC_b1, Ferroper56_LatTC_ymin, Ferroper56_LatTC_ymax,
            P_log, T_mod, T_room, Ferroper56_TDep_n_Exp);
          All_Minerals_LatTcond[Ferroper56_Index] = Ferroper56_LatTCon;
          // Compute lattice thermal conductivities for Davemaoite
          double Davemaoite_LatTCon = compute_lattice_thermal_conductivity(
            Davemaoite_LatTC_a0, Davemaoite_LatTC_b1, Davemaoite_LatTC_ymin, Davemaoite_LatTC_ymax,
            P_log, T_mod, T_room, Davemaoite_TDep_n_Exp);
          All_Minerals_LatTcond[Davemaoite_Index] = Davemaoite_LatTCon;
          // Compute lattice thermal conductivities for New-hexagonal-alluminium-phase (FeNAL)
          double NewHexAlPh_LatTCon = compute_lattice_thermal_conductivity(
            NewHexAlPh_LatTC_a0, NewHexAlPh_LatTC_b1, NewHexAlPh_LatTC_ymin, NewHexAlPh_LatTC_ymax,
            P_log, T_mod, T_room, NewHexAlPh_TDep_n_Exp);
          All_Minerals_LatTcond[NewHexAlPh_Index] = NewHexAlPh_LatTCon;
          // Compute lattice thermal conductivities for Akimotoite
          double Akimotoite_LatTCon = compute_lattice_thermal_conductivity(
            Akimotoite_LatTC_a0, Akimotoite_LatTC_b1, Akimotoite_LatTC_ymin, Akimotoite_LatTC_ymax,
            P_log, T_mod, T_room, Akimotoite_TDep_n_Exp);
          All_Minerals_LatTcond[Akimotoite_Index] = Akimotoite_LatTCon;

          // Compute radiative thermal conductivities for DryOlivine
          double OlivineDry_RadTCon = compute_radiative_thermal_conductivity(
            OlivineDry_RadTC_c0, OlivineDry_RadTC_d1, OlivineDry_RadTC_jmin, OlivineDry_RadTC_jmax, T_log);
          All_Minerals_RadTcond[OlivineDry_Index] = OlivineDry_RadTCon;
          // Compute radiative thermal conductivities for Dry Wadsleyite 
          double WadsleyDry_RadTCon = compute_radiative_thermal_conductivity(
            WadsleyDry_RadTC_c0, WadsleyDry_RadTC_d1, WadsleyDry_RadTC_jmin, WadsleyDry_RadTC_jmax, T_log);
          All_Minerals_RadTcond[WadsleyDry_Index] = WadsleyDry_RadTCon;
          // Compute radiative thermal conductivities for Dry Ringwoodite 
          double RingwooDry_RadTCon = compute_radiative_thermal_conductivity(
            RingwooDry_RadTC_c0, RingwooDry_RadTC_d1, RingwooDry_RadTC_jmin, RingwooDry_RadTC_jmax, T_log);
          All_Minerals_RadTcond[RingwooDry_Index] = RingwooDry_RadTCon;
          // Compute radiative thermal conductivities for Mg-Bridgmanite
          double En100Brigm_RadTCon = compute_radiative_thermal_conductivity(
            En100Brigm_RadTC_c0, En100Brigm_RadTC_d1, En100Brigm_RadTC_jmin, En100Brigm_RadTC_jmax, T_log);
          All_Minerals_RadTcond[En100Brigm_Index] = En100Brigm_RadTCon;
          // Compute radiative thermal conductivities for Fe-Bridgmanite (3%)
          double En97Brigma_RadTCon = compute_radiative_thermal_conductivity(
            En97Brigma_RadTC_c0, En97Brigma_RadTC_d1, En97Brigma_RadTC_jmin, En97Brigma_RadTC_jmax, T_log);
          All_Minerals_RadTcond[En97Brigma_Index] = En97Brigma_RadTCon;
          // Compute radiative thermal conductivities for Fe-Bridgmanite (10%)
          double En90Brigma_RadTCon = compute_radiative_thermal_conductivity(
            En90Brigma_RadTC_c0, En90Brigma_RadTC_d1, En90Brigma_RadTC_jmin, En90Brigma_RadTC_jmax, T_log);
          All_Minerals_RadTcond[En90Brigma_Index] = En90Brigma_RadTCon;
          // Compute radiative thermal conductivities for Al-Bridgmanite
          double AlMgBrigma_RadTCon = compute_radiative_thermal_conductivity(
            AlMgBrigma_RadTC_c0, AlMgBrigma_RadTC_d1, AlMgBrigma_RadTC_jmin, AlMgBrigma_RadTC_jmax, T_log);
          All_Minerals_RadTcond[AlMgBrigma_Index] = AlMgBrigma_RadTCon;
          // Compute radiative thermal conductivities for Fe,Al-Bridgmanite
          double FeAlBrigma_RadTCon = compute_radiative_thermal_conductivity(
            FeAlBrigma_RadTC_c0, FeAlBrigma_RadTC_d1, FeAlBrigma_RadTC_jmin, FeAlBrigma_RadTC_jmax, T_log);
          All_Minerals_RadTcond[FeAlBrigma_Index] = FeAlBrigma_RadTCon;
          // Compute radiative thermal conductivities for Orthopyroxene (Enstatite)
          double OpxEnstati_RadTCon = compute_radiative_thermal_conductivity(
            OpxEnstati_RadTC_c0, OpxEnstati_RadTC_d1, OpxEnstati_RadTC_jmin, OpxEnstati_RadTC_jmax, T_log);
          All_Minerals_RadTcond[OpxEnstati_Index] = OpxEnstati_RadTCon;
          // Compute radiative thermal conductivities for Clinopyroxene (Diopside)
          double CpxDiopsid_RadTCon = compute_radiative_thermal_conductivity(
            CpxDiopsid_RadTC_c0, CpxDiopsid_RadTC_d1, CpxDiopsid_RadTC_jmin, CpxDiopsid_RadTC_jmax, T_log);
          All_Minerals_RadTcond[CpxDiopsid_Index] = CpxDiopsid_RadTCon;
          // Compute radiative thermal conductivities for Garnet (Pyrope)
          double GrtPyropes_RadTCon = compute_radiative_thermal_conductivity(
            GrtPyropes_RadTC_c0, GrtPyropes_RadTC_d1, GrtPyropes_RadTC_jmin, GrtPyropes_RadTC_jmax, T_log);
          All_Minerals_RadTcond[GrtPyropes_Index] = GrtPyropes_RadTCon;
          // Compute radiative thermal conductivities for Garnet (Grossular)
          double GrtGrossul_RadTCon = compute_radiative_thermal_conductivity(
            GrtGrossul_RadTC_c0, GrtGrossul_RadTC_d1, GrtGrossul_RadTC_jmin, GrtGrossul_RadTC_jmax, T_log);
          All_Minerals_RadTcond[GrtGrossul_Index] = GrtGrossul_RadTCon;
          // Compute radiative thermal conductivities for Garnet (Almandine)
          double GrtAlmandi_RadTCon = compute_radiative_thermal_conductivity(
            GrtAlmandi_RadTC_c0, GrtAlmandi_RadTC_d1, GrtAlmandi_RadTC_jmin, GrtAlmandi_RadTC_jmax, T_log);
          All_Minerals_RadTcond[GrtAlmandi_Index] = GrtAlmandi_RadTCon;
          // Compute radiative thermal conductivities for Garnet (Majorite)
          double GrtMajorit_RadTCon = compute_radiative_thermal_conductivity(
            GrtMajorit_RadTC_c0, GrtMajorit_RadTC_d1, GrtMajorit_RadTC_jmin, GrtMajorit_RadTC_jmax, T_log);
          All_Minerals_RadTcond[GrtMajorit_Index] = GrtMajorit_RadTCon;
          // Compute radiative thermal conductivities for Quartz
          double QuartzPure_RadTCon = compute_radiative_thermal_conductivity(
            QuartzPure_RadTC_c0, QuartzPure_RadTC_d1, QuartzPure_RadTC_jmin, QuartzPure_RadTC_jmax, T_log);
          All_Minerals_RadTcond[QuartzPure_Index] = QuartzPure_RadTCon;
          // Compute radiative thermal conductivities for Coesite
          double CoesitSiO2_RadTCon = compute_radiative_thermal_conductivity(
            CoesitSiO2_RadTC_c0, CoesitSiO2_RadTC_d1, CoesitSiO2_RadTC_jmin, CoesitSiO2_RadTC_jmax, T_log);
          All_Minerals_RadTcond[CoesitSiO2_Index] = CoesitSiO2_RadTCon;
          // Compute radiative thermal conductivities for Stishovite
          double Stishovite_RadTCon = compute_radiative_thermal_conductivity(
            Stishovite_RadTC_c0, Stishovite_RadTC_d1, Stishovite_RadTC_jmin, Stishovite_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Stishovite_Index] = Stishovite_RadTCon;
          // Compute radiative thermal conductivities for Al-stishovite (5 vol%)
          double Al05Stisho_RadTCon = compute_radiative_thermal_conductivity(
            Al05Stisho_RadTC_c0, Al05Stisho_RadTC_d1, Al05Stisho_RadTC_jmin, Al05Stisho_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Al05Stisho_Index] = Al05Stisho_RadTCon;
          // Compute radiative thermal conductivities for Antigorite (010)
          double Antigor010_RadTCon = compute_radiative_thermal_conductivity(
            Antigor010_RadTC_c0, Antigor010_RadTC_d1, Antigor010_RadTC_jmin, Antigor010_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Antigor010_Index] = Antigor010_RadTCon;
          // Compute radiative thermal conductivities for Antigorite (001)
          double Antigor001_RadTCon = compute_radiative_thermal_conductivity(
            Antigor001_RadTC_c0, Antigor001_RadTC_d1, Antigor001_RadTC_jmin, Antigor001_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Antigor001_Index] = Antigor001_RadTCon;
          // Compute radiative thermal conductivities for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
          double FeAlPhaseD_RadTCon = compute_radiative_thermal_conductivity(
            FeAlPhaseD_RadTC_c0, FeAlPhaseD_RadTC_d1, FeAlPhaseD_RadTC_jmin, FeAlPhaseD_RadTC_jmax, T_log);
          All_Minerals_RadTcond[FeAlPhaseD_Index] = FeAlPhaseD_RadTCon;
          // Compute radiative thermal conductivities for Al-phase D (Dense Hydrous Magnesium Silicate)
          double Al02PhaseD_RadTCon = compute_radiative_thermal_conductivity(
            Al02PhaseD_RadTC_c0, Al02PhaseD_RadTC_d1, Al02PhaseD_RadTC_jmin, Al02PhaseD_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Al02PhaseD_Index] = Al02PhaseD_RadTCon;
          // Compute radiative thermal conductivities for Ferropericlase (Mg92Fe8O)
          double Ferroper08_RadTCon = compute_radiative_thermal_conductivity(
            Ferroper08_RadTC_c0, Ferroper08_RadTC_d1, Ferroper08_RadTC_jmin, Ferroper08_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Ferroper08_Index] = Ferroper08_RadTCon;
          // Compute radiative thermal conductivities for Ferropericlase (Mg90Fe10O)
          double Ferroper10_RadTCon = compute_radiative_thermal_conductivity(
            Ferroper10_RadTC_c0, Ferroper10_RadTC_d1, Ferroper10_RadTC_jmin, Ferroper10_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Ferroper10_Index] = Ferroper10_RadTCon;
          // Compute radiative thermal conductivities for Ferropericlase (Mg56Fe44O)
          double Ferroper56_RadTCon = compute_radiative_thermal_conductivity(
            Ferroper56_RadTC_c0, Ferroper56_RadTC_d1, Ferroper56_RadTC_jmin, Ferroper56_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Ferroper56_Index] = Ferroper56_RadTCon;
          // Compute radiative thermal conductivities for Davemaoite
          double Davemaoite_RadTCon = compute_radiative_thermal_conductivity(
            Davemaoite_RadTC_c0, Davemaoite_RadTC_d1, Davemaoite_RadTC_jmin, Davemaoite_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Davemaoite_Index] = Davemaoite_RadTCon;
          // Compute radiative thermal conductivities for New-hexagonal-alluminium-phase (FeNAL)
          double NewHexAlPh_RadTCon = compute_radiative_thermal_conductivity(
            NewHexAlPh_RadTC_c0, NewHexAlPh_RadTC_d1, NewHexAlPh_RadTC_jmin, NewHexAlPh_RadTC_jmax, T_log);
          All_Minerals_RadTcond[NewHexAlPh_Index] = NewHexAlPh_RadTCon;
          // Compute radiative thermal conductivities for Akimotoite
          double Akimotoite_RadTCon = compute_radiative_thermal_conductivity(
            Akimotoite_RadTC_c0, Akimotoite_RadTC_d1, Akimotoite_RadTC_jmin, Akimotoite_RadTC_jmax, T_log);
          All_Minerals_RadTcond[Akimotoite_Index] = Akimotoite_RadTCon;

          // Compute total thermal conductivities for DryOlivine
          double OlivineDry_TotTCon = compute_total_thermal_conductivity(
            OlivineDry_LatTCon, OlivineDry_RadTCon);
          All_Minerals_TotTcond[OlivineDry_Index] = OlivineDry_TotTCon;
          // Compute total thermal conductivities for Dry Wadsleyite
          double WadsleyDry_TotTCon = compute_total_thermal_conductivity(
            WadsleyDry_LatTCon, WadsleyDry_RadTCon);
          All_Minerals_TotTcond[WadsleyDry_Index] = WadsleyDry_TotTCon;
          // Compute total thermal conductivities for Dry Ringwoodite 
          double RingwooDry_TotTCon = compute_total_thermal_conductivity(
            RingwooDry_LatTCon, RingwooDry_RadTCon);
          All_Minerals_TotTcond[RingwooDry_Index] = RingwooDry_TotTCon;
          // Compute total thermal conductivities for Mg-Bridgmanite
          double En100Brigm_TotTCon = compute_total_thermal_conductivity(
            En100Brigm_LatTCon, En100Brigm_RadTCon);
          All_Minerals_TotTcond[En100Brigm_Index] = En100Brigm_TotTCon;
          // Compute total thermal conductivities for Fe-Bridgmanite (3%)
          double En97Brigma_TotTCon = compute_total_thermal_conductivity(
            En97Brigma_LatTCon, En97Brigma_RadTCon);
          All_Minerals_TotTcond[En97Brigma_Index] = En97Brigma_TotTCon;
          // Compute total thermal conductivities for Fe-Bridgmanite (10%)
          double En90Brigma_TotTCon = compute_total_thermal_conductivity(
            En90Brigma_LatTCon, En90Brigma_RadTCon);
          All_Minerals_TotTcond[En90Brigma_Index] = En90Brigma_TotTCon;
          // Compute total thermal conductivities for Al-Bridgmanite
          double AlMgBrigma_TotTCon = compute_total_thermal_conductivity(
            AlMgBrigma_LatTCon, AlMgBrigma_RadTCon);
          All_Minerals_TotTcond[AlMgBrigma_Index] = AlMgBrigma_TotTCon;
          // Compute total thermal conductivities for Fe,Al-Bridgmanite
          double FeAlBrigma_TotTCon = compute_total_thermal_conductivity(
            FeAlBrigma_LatTCon, FeAlBrigma_RadTCon);
          All_Minerals_TotTcond[FeAlBrigma_Index] = FeAlBrigma_TotTCon;
          // Compute total thermal conductivities for Orthopyroxene (Enstatite)
          double OpxEnstati_TotTCon = compute_total_thermal_conductivity(
            OpxEnstati_LatTCon, OpxEnstati_RadTCon);
          All_Minerals_TotTcond[OpxEnstati_Index] = OpxEnstati_TotTCon;
          // Compute total thermal conductivities for Clinopyroxene (Diopside)
          double CpxDiopsid_TotTCon = compute_total_thermal_conductivity(
            CpxDiopsid_LatTCon, CpxDiopsid_RadTCon);
          All_Minerals_TotTcond[CpxDiopsid_Index] = CpxDiopsid_TotTCon;
          // Compute total thermal conductivities for Garnet (Pyrope)
          double GrtPyropes_TotTCon = compute_total_thermal_conductivity(
            GrtPyropes_LatTCon, GrtPyropes_RadTCon);
          All_Minerals_TotTcond[GrtPyropes_Index] = GrtPyropes_TotTCon;
          // Compute total thermal conductivities for Garnet (Grossular)
          double GrtGrossul_TotTCon = compute_total_thermal_conductivity(
            GrtGrossul_LatTCon, GrtGrossul_RadTCon);
          All_Minerals_TotTcond[GrtGrossul_Index] = GrtGrossul_TotTCon;
          // Compute total thermal conductivities for Garnet (Almandine)
          double GrtAlmandi_TotTCon = compute_total_thermal_conductivity(
            GrtAlmandi_LatTCon, GrtAlmandi_RadTCon);
          All_Minerals_TotTcond[GrtAlmandi_Index] = GrtAlmandi_TotTCon;
          // Compute total thermal conductivities for Garnet (Majorite)
          double GrtMajorit_TotTCon = compute_total_thermal_conductivity(
            GrtMajorit_LatTCon, GrtMajorit_RadTCon);
          All_Minerals_TotTcond[GrtMajorit_Index] = GrtMajorit_TotTCon;
          // Compute total thermal conductivities for Quartz
          double QuartzPure_TotTCon = compute_total_thermal_conductivity(
            QuartzPure_LatTCon, QuartzPure_RadTCon);
          All_Minerals_TotTcond[QuartzPure_Index] = QuartzPure_TotTCon;
          // Compute total thermal conductivities for Coesite
          double CoesitSiO2_TotTCon = compute_total_thermal_conductivity(
            CoesitSiO2_LatTCon, CoesitSiO2_RadTCon);
          All_Minerals_TotTcond[CoesitSiO2_Index] = CoesitSiO2_TotTCon;
          // Compute total thermal conductivities for Stishovite
          double Stishovite_TotTCon = compute_total_thermal_conductivity(
            Stishovite_LatTCon, Stishovite_RadTCon);
          All_Minerals_TotTcond[Stishovite_Index] = Stishovite_TotTCon;
          // Compute total thermal conductivities for Al-stishovite (5 vol%)
          double Al05Stisho_TotTCon = compute_total_thermal_conductivity(
            Al05Stisho_LatTCon, Al05Stisho_RadTCon);
          All_Minerals_TotTcond[Al05Stisho_Index] = Al05Stisho_TotTCon;
          // Compute total thermal conductivities for Antigorite (010)
          double Antigor010_TotTCon = compute_total_thermal_conductivity(
            Antigor010_LatTCon, Antigor010_RadTCon);
          All_Minerals_TotTcond[Antigor010_Index] = Antigor010_TotTCon;
          // Compute total thermal conductivities for Antigorite (001)
          double Antigor001_TotTCon = compute_total_thermal_conductivity(
            Antigor001_LatTCon, Antigor001_RadTCon);
          All_Minerals_TotTcond[Antigor001_Index] = Antigor001_TotTCon;
          // Compute total thermal conductivities for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
          double FeAlPhaseD_TotTCon = compute_total_thermal_conductivity(
            FeAlPhaseD_LatTCon, FeAlPhaseD_RadTCon);
          All_Minerals_TotTcond[FeAlPhaseD_Index] = FeAlPhaseD_TotTCon;
          // Compute total thermal conductivities for Al-phase D (Dense Hydrous Magnesium Silicate)
          double Al02PhaseD_TotTCon = compute_total_thermal_conductivity(
            Al02PhaseD_LatTCon, Al02PhaseD_RadTCon);
          All_Minerals_TotTcond[Al02PhaseD_Index] = Al02PhaseD_TotTCon;
          // Compute total thermal conductivities for Ferropericlase (Mg92Fe8O)
          double Ferroper08_TotTCon = compute_total_thermal_conductivity(
            Ferroper08_LatTCon, Ferroper08_RadTCon);
          All_Minerals_TotTcond[Ferroper08_Index] = Ferroper08_TotTCon;
          // Compute total thermal conductivities for Ferropericlase (Mg90Fe10O)
          double Ferroper10_TotTCon = compute_total_thermal_conductivity(
            Ferroper10_LatTCon, Ferroper10_RadTCon);
          All_Minerals_TotTcond[Ferroper10_Index] = Ferroper10_TotTCon;
          // Compute total thermal conductivities for Ferropericlase (Mg56Fe44O)
          double Ferroper56_TotTCon = compute_total_thermal_conductivity(
            Ferroper56_LatTCon, Ferroper56_RadTCon);
          All_Minerals_TotTcond[Ferroper56_Index] = Ferroper56_TotTCon;
          // Compute total thermal conductivities for Davemaoite
          double Davemaoite_TotTCon = compute_total_thermal_conductivity(
            Davemaoite_LatTCon, Davemaoite_RadTCon);
          All_Minerals_TotTcond[Davemaoite_Index] = Davemaoite_TotTCon;
          // Compute total thermal conductivities for New-hexagonal-alluminium-phase (FeNAL)
          double NewHexAlPh_TotTCon = compute_total_thermal_conductivity(
            NewHexAlPh_LatTCon, NewHexAlPh_RadTCon);
          All_Minerals_TotTcond[NewHexAlPh_Index] = NewHexAlPh_TotTCon;
          // Compute total thermal conductivities for Akimotoite
          double Akimotoite_TotTCon = compute_total_thermal_conductivity(
            Akimotoite_LatTCon, Akimotoite_RadTCon);
          All_Minerals_TotTcond[Akimotoite_Index] = Akimotoite_TotTCon;
    
          // Fill the matrix column by column
          for (unsigned int row = 0; row < MineralPar_Index; ++row)
          {
            All_Minerals_TConds[row][0] = All_Minerals_LatTcond[row]; // Column 0: Lattice conductivities
            All_Minerals_TConds[row][1] = All_Minerals_RadTcond[row]; // Column 1: Radiative conductivities
            All_Minerals_TConds[row][2] = All_Minerals_TotTcond[row]; // Column 2: Total conductivities
          }

          // Compute P,T-dependent thermal conductivities of aggregate rocks 

          // Pyrolite (58% olivine, 13% pyrope, 18% ensatite, 11% diopside)
          // double AggRock_Pyrolite_TCond = std::pow(OlivineDry_TotTCon,MinFract_Pyrolite[0])*std::pow(GrtPyropes_TotTCon,MinFract_Pyrolite[1])*std::pow(OpxEnstati_TotTCon,MinFract_Pyrolite[2])*std::pow(CpxDiopsid_TotTCon,MinFract_Pyrolite[3]);   
          // Harzburgite (80% olivine, 20% ensatite)
          // double AggRock_Harzburg_TCond = std::pow(OlivineDry_TotTCon,MinFract_Harzburg[0])*std::pow(OpxEnstati_TotTCon,MinFract_Harzburg[1]);
          // Meta-basaltic crust MORB (80% diopside, 20% pyrope)
          // double AggRock_MetaMORB_TCond = std::pow(CpxDiopsid_TotTCon,MinFract_MetaMORB[0])*std::pow(GrtPyropes_TotTCon,MinFract_MetaMORB[1]);
          // Dunite (100% olivine)
          // double AggRock_DuniteOl_TCond = std::pow(OlivineDry_TotTCon,MinFract_DuniteOl[0]);
          // Test Case
          double AggRock_TestCase_TCond = std::pow(All_Minerals_TConds[10][2], min_frac);

          out.thermal_conductivities[i] = AggRock_TestCase_TCond;
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
  template class PTdepRbounded<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
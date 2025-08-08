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

#include <aspect/material_model/thermal_conductivity/marzotto_2025.h>

// Helper functions in anonymous namespace to compute thermal conductivities using the Marzotto et al. (2025) formulations
namespace 
{
  // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_lattice_thermal_conductivity_mar2025(double a0, double b1, double ymin, double ymax, double P_log, double T_mod, double T_room, double n_exp)
  { 
    double zsimpl = a0 + b1 * P_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = ymin + (ymax - ymin) * yprime;
    double lattice_thermal_conductivity = std::exp(yreals);
    return lattice_thermal_conductivity * std::pow((T_room / T_mod), n_exp);
  }
  
  // Compute the radiative thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_radiative_thermal_conductivity_mar2025(double c0, double d1, double jmin, double jmax, double T_log)
  {
    double zsimpl = c0 + d1 * T_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = jmin + (jmax - jmin) * yprime;
    return std::exp(yreals);
  }
     
  double compute_total_thermal_conductivity_mar2025(double lattice_thermal_conductivity, double radiative_thermal_conductivity)
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
      marzotto_2025<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        // Define coefficients for lattice thermal conductivity for different minerals 

        // Coefficients for dry olivine 
        // retreived from fitting TDTR dataset of
        // [Chang et al., 2017, PNAS, vol 114, p. 4078-4081]
        // https://doi.org/10.1073/pnas.1616216114
        // mineral composition [Mg1.8 Fe0.2 SiO4]    
        constexpr int olivinedry_index = 0;
        const double olivinedry_latTC_a0 =   -4.1241;
        const double olivinedry_latTC_b1 =    2.1469;
        const double olivinedry_latTC_ymin =  1.28093384543429;
        const double olivinedry_latTC_ymax =  2.60726381956037;
        const double olivinedry_Tdep_n_exp =  0.5;    
        
        // Coefficients for dry wadsleyite
        // retreived from fitting dataset of
        // [Xu et al., 2004, PEPI, vol 143, pp. 321-336]
        // https://doi.org/10.1016/j.pepi.2004.03.005
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        constexpr int wadsleydry_index = 1;
        const double wadsleydry_latTC_a0 =  -0.6656;
        const double wadsleydry_latTC_b1 =   0.3082;
        const double wadsleydry_latTC_ymin = 1.75735804249439;
        const double wadsleydry_latTC_ymax = 2.37090451537473; 
        const double wadsleydry_Tdep_n_exp = 0.5; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Marzotto et al., 2020, GRL, vol 47, issue 13]
        // https://doi.org/10.1029/2020GL087607
        // mineral composition [(Mg1.79Fe0.17)Si1.02O4]
        constexpr int ringwoodry_index = 2;
        const double ringwoodry_latTC_a0 =  -5.4624;
        const double ringwoodry_latTC_b1 =   2.0791;
        const double ringwoodry_latTC_ymin = 1.60943791241410;
        const double ringwoodry_latTC_ymax = 2.94939766245070; 
        const double ringwoodry_Tdep_n_exp = 0.5; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [MgSiO3]
        constexpr int brigm100Mg_index = 3;
        const double brigm100Mg_latTC_a0 =  -4.3687;
        const double brigm100Mg_latTC_b1 =   1.0766; 
        const double brigm100Mg_latTC_ymin = 2.376025820; 
        const double brigm100Mg_latTC_ymax = 5.010635294;  
        const double brigm100Mg_Tdep_n_exp = 1.01000;

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.03Mg0.97SiO3]
        constexpr int brigma97Mg_index = 4;
        const double brigma97Mg_latTC_a0 =  -4.520600000;
        const double brigma97Mg_latTC_b1 =   1.019900000; 
        const double brigma97Mg_latTC_ymin = 1.750524121; 
        const double brigma97Mg_latTC_ymax = 4.499809670;  
        const double brigma97Mg_Tdep_n_exp = 0.56605;

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.1Mg0.9SiO3]
        constexpr int brigma90Mg_index = 5;
        const double brigma90Mg_latTC_a0 =  -4.883100000;
        const double brigma90Mg_latTC_b1 =   0.980900000; 
        const double brigma90Mg_latTC_ymin = 1.333739493; 
        const double brigma90Mg_latTC_ymax = 4.382026635;  
        const double brigma90Mg_Tdep_n_exp = 0.17054;

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Al,Mg)SiO3]
        constexpr int brigmaAlMg_index = 6;
        const double brigmaAlMg_latTC_a0 =  -4.331500000;
        const double brigmaAlMg_latTC_b1 =   1.027000000; 
        const double brigmaAlMg_latTC_ymin = 1.845020046; 
        const double brigmaAlMg_latTC_ymax = 4.605170186;  
        const double brigmaAlMg_Tdep_n_exp = 0.61983;

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Fe,Al,Mg)SiO3]
        constexpr int brigmaFeAl_index = 7;
        const double brigmaFeAl_latTC_a0 =  -4.510600000;
        const double brigmaFeAl_latTC_b1 =   1.066800000; 
        const double brigmaFeAl_latTC_ymin = 1.389093953; 
        const double brigmaFeAl_latTC_ymax = 3.912023005;  
        const double brigmaFeAl_Tdep_n_exp = 0.46815;
         
        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of 
        // [Schloessin & Dvorak, 1972, GJI, 27(5), 499-516]
        // https://doi.org/10.1111/j.1365-246X.1972.tb06105.x
        // mineral composition [Mg2Si2O6]
        constexpr int opxenstati_index = 8;
        const double opxenstati_latTC_a0 =   -3.0047;
        const double opxenstati_latTC_b1 =    2.6;
        const double opxenstati_latTC_ymin =  1.760865151; 
        const double opxenstati_latTC_ymax =  2.096937429;
        const double opxenstati_Tdep_n_exp =  0.5;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of 
        // [Wang et al., 2014, JGR: Solid Earth, 119(8), 6277-6287]
        // https://doi.org/10.1002/2014JB011208
        // mineral composition [CaMgSi2O6]
        constexpr int cpxdiopsid_index = 9;
        const double cpxdiopsid_latTC_a0 =   -3.251100000;
        const double cpxdiopsid_latTC_b1 =    1.689100000;
        const double cpxdiopsid_latTC_ymin =  1.793640135; 
        const double cpxdiopsid_latTC_ymax =  2.389462023;
        const double cpxdiopsid_Tdep_n_exp =  0.5;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [Mg3Al2Si3O12]
        constexpr int grtpyropes_index = 10;
        const double grtpyropes_latTC_a0 =   -4.3637;
        const double grtpyropes_latTC_b1 =    2.0368;
        const double grtpyropes_latTC_ymin =  1.481604541; 
        const double grtpyropes_latTC_ymax =  2.443131606;
        const double grtpyropes_Tdep_n_exp =  0.4314;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Ca0.986Fe0.014)3Al2(SiO4)3]
        constexpr int grtgrossul_index = 11;
        const double grtgrossul_latTC_a0 =  -4.7584;
        const double grtgrossul_latTC_b1 =   2.0816;
        const double grtgrossul_latTC_ymin = 1.410986974; 
        const double grtgrossul_latTC_ymax = 2.457625992;
        const double grtgrossul_Tdep_n_exp = 0.4589;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3]
        constexpr int grtalmandi_index = 12;
        const double grtalmandi_latTC_a0 =  -4.5047;
        const double grtalmandi_latTC_b1 =   2.0988;
        const double grtalmandi_latTC_ymin = 1.223775432; 
        const double grtalmandi_latTC_ymax = 2.374762159;
        const double grtalmandi_Tdep_n_exp = 0.4172;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Giesting et al.2004  EPSL, 218(1-2), 45-56]
        // https://doi.org/10.1016/S0012-821X(03)00630-7
        // mineral composition [Mg3(MgSi)(SiO4)3]
        constexpr int grtmajorit_index = 13;
        const double grtmajorit_latTC_a0 =  -4.3637;
        const double grtmajorit_latTC_b1 =   2.0368;
        const double grtmajorit_latTC_ymin = 2.279316466; 
        const double grtmajorit_latTC_ymax = 2.718047842;
        const double grtmajorit_Tdep_n_exp = 0.5;

        // Coefficients for quartz 
        // retreived from fitting dataset of
        // [Xiong et al., 2019 - Journal of Applied Physics, 126(21)]
        // https://doi.org/10.1063/1.5114992
        // mineral composition [SiO2]
        constexpr int quartzpure_index = 14;
        const double quartzpure_latTC_a0 =   -2.0203;
        const double quartzpure_latTC_b1 =    2.4456;
        const double quartzpure_latTC_ymin =  2.260981081; 
        const double quartzpure_latTC_ymax =  2.745391462;
        const double quartzpure_Tdep_n_exp =  1.015433333;

        // Coefficients for coesite
        // retreived from fitting dataset of
        // [Yukutake & Shimada, 1978, PEPI, 17(3), 193-200]
        // https://doi.org/10.1016/0031-9201(78)90036-5
        // mineral composition [SiO2]
        constexpr int coesitSiO2_index = 15;
        const double coesitSiO2_latTC_a0 =   -12.728;
        const double coesitSiO2_latTC_b1 =    2.9998;
        const double coesitSiO2_latTC_ymin =  1.982022416; 
        const double coesitSiO2_latTC_ymax =  2.249036030;
        const double coesitSiO2_Tdep_n_exp =  1.015433333;

        // Coefficients for stishovite
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584, 117477]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [SiO2]
        constexpr int stishovite_index = 16;
        // Assign coefficients based on pressure ranges
        // Pressure < 52 [GPa]
        const double stishovite_latTC_1_a0 =  16.917;
        const double stishovite_latTC_1_b1 = -4.6187;
        const double stishovite_latTC_1_ymin = 4.096113064; 
        const double stishovite_latTC_1_ymax = 4.217974805;
        // Pressure between 52 and 56 [GPa]
        const double stishovite_latTC_2_a0 = -156.12;
        const double stishovite_latTC_2_b1 =  39.182;
        const double stishovite_latTC_2_ymin = 4.077505176;
        const double stishovite_latTC_2_ymax = 4.264199335;
        // Pressure > 56 [GPa]
        const double stishovite_latTC_3_a0 = -12.728;
        const double stishovite_latTC_3_b1 =  2.9998;
        const double stishovite_latTC_3_ymin = 3.960844211;
        const double stishovite_latTC_3_ymax = 4.738489125;
        // Temperature
        const double stishovite_Tdep_n_exp = 0.5;

        // Coefficients for Al-stishovite (5 vol%)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584, 117477]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [(Al,Si)O2]
        constexpr int stisho05Al_index = 17;
        const double stisho05Al_latTC_a0 = -6.4411;
        const double stisho05Al_latTC_b1 =  1.5885;
        const double stisho05Al_latTC_ymin = 3.188855035;
        const double stisho05Al_latTC_ymax = 4.154336189;
        const double stisho05Al_Tdep_n_exp = 0.5;

        // Coefficients for antigorite (serpentine)
        // retreived from fitting dataset of
        // [Chien et al., 2024, Nature Communications, 15(1), 5198.]
        // https://doi.org/10.1038/s41467-024-49418-3
        // mineral composition [(Mg2.80Fe0.05)Si2.08O5(OH)3.77]
        // 010 direction
        constexpr int antigor010_index = 18;
        const double antigor010_latTC_a0 = -4.3374;
        const double antigor010_latTC_b1 =  2.0217;
        const double antigor010_latTC_ymin = 1.519513205;
        const double antigor010_latTC_ymax = 2.434491480;
        const double antigor010_Tdep_n_exp = 0.5;
        // 001 direction
        constexpr int antigor001_index = 19;
        const double antigor001_latTC_a0 = -3.1109;
        const double antigor001_latTC_b1 =  2.0644;
        const double antigor001_latTC_ymin = 0.067658648;
        const double antigor001_latTC_ymax = 1.552797578;
        const double antigor001_Tdep_n_exp = 0.5;

        // Coefficients for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, JGR: Solid Earth, vol. 127(6), e2022JB024556]
        // https://doi.org/10.1029/2022JB024556
        // mineral composition [Mg1.19Fe0.12Al0.174Si1.71H2.02O6]
        constexpr int phaseDFeAl_index = 20;
        // (Fe,Al)-Phase D - 0-24 [GPa]
        const double phaseDFeAl_latTC_1_a0 = -3.9909;
        const double phaseDFeAl_latTC_1_b1 =  1.7710;
        const double phaseDFeAl_latTC_1_ymin = 0.956005323; 
        const double phaseDFeAl_latTC_1_ymax = 1.747361025;
        // (Fe,Al)-Phase D - 24-38 [GPa]
        const double phaseDFeAl_latTC_2_a0 = -32.890;
        const double phaseDFeAl_latTC_2_b1 =  9.6282;
        const double phaseDFeAl_latTC_2_ymin = 1.442700096; 
        const double phaseDFeAl_latTC_2_ymax = 3.512683919;
        // (Fe,Al)-Phase D - 38-48 [GPa]
        const double phaseDFeAl_latTC_3_a0 =  141.88;
        const double phaseDFeAl_latTC_3_b1 = -37.409;
        const double phaseDFeAl_latTC_3_ymin = 1.789742436; 
        const double phaseDFeAl_latTC_3_ymax = 3.270312073;
        // (Fe,Al)-Phase D - > 48 [GPa]
        const double phaseDFeAl_latTC_4_a0 = -23.986;
        const double phaseDFeAl_latTC_4_b1 =  6.1139;
        const double phaseDFeAl_latTC_4_ymin = 1.313988596; 
        const double phaseDFeAl_latTC_4_ymax = 2.992561000;
        // Temperature-dependency
        const double phaseDFeAl_Tdep_n_exp = 0.5;

        // Coefficients for Al-phase D (Dense Hydrous Magnesium Silicate)
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, JGR: Solid Earth, vol. 127(6), e2022JB024556]
        // https://doi.org/10.1029/2022JB024556
        // mineral composition [Mg1.29Al0.17Si1.73H1.98O6]
        constexpr int phaseD02Al_index = 21;
        const double phaseD02Al_latTC_a0 = -6.1829;
        const double phaseD02Al_latTC_b1 =  1.8514;
        const double phaseD02Al_latTC_ymin = 1.285874399; 
        const double phaseD02Al_latTC_ymax = 3.502412041;
        const double phaseD02Al_Tdep_n_exp = 0.5;

        // Coefficients for ferropericlase (Mg1-xFexO)
        // retreived from fitting dataset of
        // [Hsieh et al., 2018, PNAS, vol. 115, no. 16, p. 4099-4104] - 8%,10%,56% Iron
        // https://doi.org/10.1073/pnas.1718557115
        // [Zhang et al., 2023., GRL, 50(7), e2022GL101769] - 20% Iron
        // https://doi.org/10.1029/2022GL101769
        // mineral composition [Mg0.92Fe0.08O] - (8% Iron)
        constexpr int ferroper08_index = 22;
        const double ferroper08_latTC_a0 = -6.9942;
        const double ferroper08_latTC_b1 =  1.953;
        const double ferroper08_latTC_ymin = 1.629240539; 
        const double ferroper08_latTC_ymax = 4.118362306;
        const double ferroper08_Tdep_n_exp = 0.5;
        // mineral composition [Mg0.90Fe0.10O] - (10% Iron)
        constexpr int ferroper10_index = 23;
        const double ferroper10_latTC_a0 = -7.0133;
        const double ferroper10_latTC_b1 =  1.9321;
        const double ferroper10_latTC_ymin = 1.5040773968; 
        const double ferroper10_latTC_ymax = 4.0250359042;
        const double ferroper10_Tdep_n_exp = 0.5;
        // mineral composition [Mg0.80Fe0.20O] (20% Iron)
        constexpr int ferroper20_index = 24;
        const double ferroper20_latTC_a0 = -5.2408;
        const double ferroper20_latTC_b1 =  0.9649;
        const double ferroper20_latTC_ymin = 1.2490430868; 
        const double ferroper20_latTC_ymax = 3.9318256327;
        const double ferroper20_Tdep_n_exp = 0.025;
        // mineral composition [Mg0.44Fe0.56O] (56% Iron)
        constexpr int ferroper56_index = 25;
        const double ferroper56_latTC_a0 = -3.8298;
        const double ferroper56_latTC_b1 =  1.1507;
        const double ferroper56_latTC_ymin = 0.993251773; 
        const double ferroper56_latTC_ymax = 3.592193222;
        const double ferroper56_Tdep_n_exp = 0.5;

        // Coefficients for davemaoite 
        // retreived from fitting dataset of
        // [Zhang et al., 2021, Physical Review B, vol. 104, 184101]
        // https://doi.org/10.1103/PhysRevB.104.184101
        // mineral composition [CaSiO3]
        constexpr int davemaoite_index = 26;
        const double davemaoite_latTC_a0 = -4.7377;
        const double davemaoite_latTC_b1 =  1.3661;
        const double davemaoite_latTC_ymin = 2.388762789; 
        const double davemaoite_latTC_ymax = 4.045106030;
        const double davemaoite_Tdep_n_exp = 0.5;

        // Coefficients for new-hexagonal-alluminium-phase (FeNAL) 
        // retreived from fitting dataset of
        // [Hsieh et al., 2022, EPSL, vol. 584]
        // https://doi.org/10.1016/j.epsl.2022.117477
        // mineral composition [Na0.71Mg2.05Al4.62Si1.16Fe(2+)0.09Fe(3+)0.17O12]
        constexpr int newhexAlph_index = 27;
        const double newhexAlph_latTC_a0 = -29.421;
        const double newhexAlph_latTC_b1 =  7.7792;
        const double newhexAlph_latTC_ymin = 2.363551955; 
        const double newhexAlph_latTC_ymax = 3.653998874;
        const double newhexAlph_Tdep_n_exp = 0.5;

        // Coefficients for akimotoite
        // assumed to be equal to En100-Bridgmanite
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [MgSiO3]
        constexpr int akimotoite_index = 28;
        const double akimotoite_latTC_a0 =  -4.368700000;
        const double akimotoite_latTC_b1 =   1.076600000; 
        const double akimotoite_latTC_ymin = 2.376025820; 
        const double akimotoite_latTC_ymax = 5.010635294;  
        const double akimotoite_Tdep_n_exp = 1.01000;

        unsigned int mineralpar_index = akimotoite_index+1; // Number of minerals

        // Define coefficients for radiative thermal conductivity of different minerals

        // Coefficients for dry olivine
        // retreived from fitting dataset of
        // [Marzotto et al. 2025, Nature Communication, 16, 6058]
        // https://doi.org/10.1038/s41467-025-61148-8
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double olivinedry_radTC_c0 =   -10.00900000;
        const double olivinedry_radTC_d1 =    1.883900000;
        const double olivinedry_radTC_jmin = -23.02585093;
        const double olivinedry_radTC_jmax =  1.289885976;

        // Coefficients for dry wadsleyite
        // retreived from fitting dataset of
        // [Thomas et al., 2012, EPSL, vol. 357, p. 130-136.]
        // https://doi.org/10.1016/j.epsl.2012.09.035
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double wadsleydry_radTC_c0 =   -21.717;
        const double wadsleydry_radTC_d1 =    3.4271;
        const double wadsleydry_radTC_jmin = -23.025850930;
        const double wadsleydry_radTC_jmax =  1.0773006810; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Thomas et al., 2012, EPSL, vol. 357, p. 130-136.]
        // https://doi.org/10.1016/j.epsl.2012.09.035
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double ringwoodry_radTC_c0 =   -23.067000000;
        const double ringwoodry_radTC_d1 =    3.5985000000;
        const double ringwoodry_radTC_jmin = -23.025850930;
        const double ringwoodry_radTC_jmax =  0.4027750000; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [MgSiO3]
        const double brigm100Mg_radTC_c0 =  66.278;
        const double brigm100Mg_radTC_d1 = -8.2756; 
        const double brigm100Mg_radTC_jmin = -7.2568958208;         
        const double brigm100Mg_radTC_jmax = -0.3403920329;         

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Fe0.03Mg0.97SiO3]
        const double brigma97Mg_radTC_c0 =  66.278;
        const double brigma97Mg_radTC_d1 = -8.2756; 
        const double brigma97Mg_radTC_jmin = -7.2568958208;         
        const double brigma97Mg_radTC_jmax = -0.3403920329;  

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Fe0.1Mg0.9SiO3]
        const double brigma90Mg_radTC_c0 =  66.278;
        const double brigma90Mg_radTC_d1 = -8.2756; 
        const double brigma90Mg_radTC_jmin = -7.2568958208;         
        const double brigma90Mg_radTC_jmax = -0.3403920329;  

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [(Al,Mg)SiO3]
        const double brigmaAlMg_radTC_c0 =  66.278;
        const double brigmaAlMg_radTC_d1 = -8.2756; 
        const double brigmaAlMg_radTC_jmin = -7.2568958208;         
        const double brigmaAlMg_radTC_jmax = -0.3403920329;  

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [(Fe,Al,Mg)SiO3]
        const double brigmaFeAl_radTC_c0 =  66.278;
        const double brigmaFeAl_radTC_d1 = -8.2756; 
        const double brigmaFeAl_radTC_jmin = -7.2568958208;         
        const double brigmaFeAl_radTC_jmax = -0.3403920329;  

        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg2Si2O6]
        const double opxenstati_radTC_c0 =   -13.532000000;
        const double opxenstati_radTC_d1 =    2.4004000000;
        const double opxenstati_radTC_jmin = -23.025850930; 
        const double opxenstati_radTC_jmax =  1.4456685920;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [CaMgSi2O6]
        const double cpxdiopsid_radTC_c0 =   -14.286000000;       
        const double cpxdiopsid_radTC_d1 =    2.5119000000; 
        const double cpxdiopsid_radTC_jmin = -23.025850930; 
        const double cpxdiopsid_radTC_jmax =  1.4434214740;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3Al2(SiO4)3]
        const double grtpyropes_radTC_c0 =   -11.782000000;
        const double grtpyropes_radTC_d1 =    2.0718000000;
        const double grtpyropes_radTC_jmin = -23.025850930; 
        const double grtpyropes_radTC_jmax =  1.4479836950;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Ca3Al2(SiO4)3]
        const double grtgrossul_radTC_c0 =   -11.782000000;
        const double grtgrossul_radTC_d1 =    2.0718000000;
        const double grtgrossul_radTC_jmin = -23.025850930; 
        const double grtgrossul_radTC_jmax =  1.4479836950;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3 ]
        const double grtalmandi_radTC_c0  =  -11.782000000;
        const double grtalmandi_radTC_d1 =    2.0718000000;
        const double grtalmandi_radTC_jmin = -23.025850930; 
        const double grtalmandi_radTC_jmax =  1.4479836950;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3(MgSi)(SiO4)3]
        const double grtmajorit_radTC_c0  =  -11.782000000;
        const double grtmajorit_radTC_d1 =    2.0718000000;
        const double grtmajorit_radTC_jmin = -23.025850930;
        const double grtmajorit_radTC_jmax =  1.4479836950;

        // Coefficients for quartz
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double quartzpure_radTC_c0  =   0;
        const double quartzpure_radTC_d1  =   0;
        const double quartzpure_radTC_jmin = -23.025850930; 
        const double quartzpure_radTC_jmax = -23.050000000;

        // Coefficients for coesite
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double coesitSiO2_radTC_c0  =   0;
        const double coesitSiO2_radTC_d1  =   0;
        const double coesitSiO2_radTC_jmin = -23.025850930; 
        const double coesitSiO2_radTC_jmax = -23.050000000;

        // Coefficients for stishovite
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        const double stishovite_radTC_c0  =   0;
        const double stishovite_radTC_d1  =   0;
        const double stishovite_radTC_jmin = -23.025850930; 
        const double stishovite_radTC_jmax = -23.050000000;

        // Coefficients for Al-stishovite (5 vol%)
        // assumed 0 for now - no data available
        // mineral composition [(Al,Si)O2]
        const double stisho05Al_radTC_c0  =   0;
        const double stisho05Al_radTC_d1  =   0;
        const double stisho05Al_radTC_jmin = -23.025850930; 
        const double stisho05Al_radTC_jmax = -23.050000000;

        // Coefficients for antigorite (serpentine)
        // assumed 0 for now - no data available
        // mineral composition [///]
        // 010 direction
        const double antigor010_radTC_c0  =   0;
        const double antigor010_radTC_d1  =   0;
        const double antigor010_radTC_jmin = -23.025850930; 
        const double antigor010_radTC_jmax = -23.050000000;
        // 001 direction
        const double antigor001_radTC_c0  =   0;
        const double antigor001_radTC_d1  =   0;
        const double antigor001_radTC_jmin = -23.025850930; 
        const double antigor001_radTC_jmax = -23.050000000;

        // Coefficients for Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
        // assumed 0 for now - no data available
        // mineral composition [Mg1.19Fe0.12Al0.174Si1.71H2.02O6]
        const double phaseDFeAl_radTC_c0  =   0;
        const double phaseDFeAl_radTC_d1  =   0;
        const double phaseDFeAl_radTC_jmin = -23.025850930; 
        const double phaseDFeAl_radTC_jmax = -23.050000000;

        // Coefficients for Al-phase D (Dense Hydrous Magnesium Silicate)
        // assumed 0 for now - no data available
        // mineral composition [Mg1.29Al0.17Si1.73H1.98O6]
        const double phaseD02Al_radTC_c0  =   0;
        const double phaseD02Al_radTC_d1  =   0;
        const double phaseD02Al_radTC_jmin = -23.025850930; 
        const double phaseD02Al_radTC_jmax = -23.050000000;

        // Coefficients for ferropericlase (Mg1-xFexO)
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Mg0.92Fe0.08O] - (8% Iron)
        const double ferroper08_radTC_c0 =  66.278;
        const double ferroper08_radTC_d1 = -8.2756; 
        const double ferroper08_radTC_jmin = -7.2568958208;         
        const double ferroper08_radTC_jmax = -0.3403920329;  
        // mineral composition [Mg0.90Fe0.10O] - (10% Iron)
        const double ferroper10_radTC_c0 =  66.278;
        const double ferroper10_radTC_d1 = -8.2756; 
        const double ferroper10_radTC_jmin = -7.2568958208;         
        const double ferroper10_radTC_jmax = -0.3403920329; 
        // mineral composition [Mg0.80Fe0.20O] - (20% Iron)
        const double ferroper20_radTC_c0 =  66.278;
        const double ferroper20_radTC_d1 = -8.2756; 
        const double ferroper20_radTC_jmin = -7.2568958208;         
        const double ferroper20_radTC_jmax = -0.3403920329;  
        // mineral composition [Mg0.44Fe0.56O] (56% Iron)
        const double ferroper56_radTC_c0 =  66.278;
        const double ferroper56_radTC_d1 = -8.2756; 
        const double ferroper56_radTC_jmin = -7.2568958208;         
        const double ferroper56_radTC_jmax = -0.3403920329;  

        // Coefficients for davemaoite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [CaSiO3] 
        const double davemaoite_radTC_c0 =  66.278;
        const double davemaoite_radTC_d1 = -8.2756; 
        const double davemaoite_radTC_jmin = -7.2568958208;         
        const double davemaoite_radTC_jmax = -0.3403920329;  

        // Coefficients for new-hexagonal-alluminium-phase (FeNAL) 
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [Na0.71Mg2.05Al4.62Si1.16Fe(2+)0.09Fe(3+)0.17O12] 
        const double newhexAlph_radTC_c0 =  66.278;
        const double newhexAlph_radTC_d1 = -8.2756; 
        const double newhexAlph_radTC_jmin = -7.2568958208;         
        const double newhexAlph_radTC_jmax = -0.3403920329;  

        // Coefficients for akimotoite
        // retreived from fitting dataset of
        // [Lobanov et al., 2020, EPSL, vol. 537, 116176]
        // https://doi.org/10.1016/j.epsl.2020.116176
        // mineral composition [MgSiO3] 
        const double akimotoite_radTC_c0 =  66.278;
        const double akimotoite_radTC_d1 = -8.2756; 
        const double akimotoite_radTC_jmin = -7.2568958208;         
        const double akimotoite_radTC_jmax = -0.3403920329;  

        // Preallocate a vector for mineral fractions of different rocks and a vector for mineral indices

        // pyrolite Upper Mantle (58% olivine, 13% pyrope, 17% ensatite, 12% diopside)
        std::vector<double> minfract_pyrolite_UM = {0.58,0.13,0.17,0.12}; 
        std::vector<unsigned int> minindex_pyrolite_UM = {olivinedry_index, grtpyropes_index, opxenstati_index, cpxdiopsid_index}; 
        // pyrolite Upper Mantle Transition Zone (58% wadsleyite, 28% majorite, 14% diopside)
        std::vector<double> minfract_pyrolite_UMTZ = {0.58,0.28,0.14}; 
        std::vector<unsigned int> minindex_pyrolite_UMTZ = {wadsleydry_index, grtmajorit_index, cpxdiopsid_index}; 
        // pyrolite Lower Mantle Transition Zone (58% ringwoodite, 42% majorite)
        std::vector<double> minfract_pyrolite_LMTZ = {0.58,0.42}; 
        std::vector<unsigned int> minindex_pyrolite_LMTZ = {ringwoodry_index, grtmajorit_index}; 
        // pyrolite Lower Mantle (80% bridgmanite, 14% ferropericlase, 6% davemaoite)
        std::vector<double> minfract_pyrolite_LM = {0.80,0.14,0.06}; 
        std::vector<unsigned int> minindex_pyrolite_LM = {brigmaAlMg_index, ferroper10_index, davemaoite_index}; 

        // harzburgite Upper Mantle (80% olivine, 20% ensatite)
        std::vector<double> minfract_harzburg_UM = {0.80,0.20}; 
        std::vector<unsigned int> minindex_harzburg_UM = {olivinedry_index, opxenstati_index};
        // harzburgite Upper Mantle Transition Zone (80% wadsleyite, 13% diopside, 7% majorite)
        std::vector<double> minfract_harzburg_UMTZ = {0.80,0.13,0.07}; 
        std::vector<unsigned int> minindex_harzburg_UMTZ = {wadsleydry_index, cpxdiopsid_index, grtmajorit_index};
        // harzburgite Lower Mantle Transition Zone (80% ringwoodite, 20% majorite)
        std::vector<double> minfract_harzburg_LMTZ = {0.80,0.20}; 
        std::vector<unsigned int> minindex_harzburg_LMTZ = {ringwoodry_index, grtmajorit_index};
        // harzburgite Lower Mantle (76% bridgmanite, 24% ferropericlase)
        std::vector<double> minfract_harzburg_LM = {0.76,0.24}; 
        std::vector<unsigned int> minindex_harzburg_LM = {brigmaAlMg_index, ferroper10_index};

        // Meta-basaltic crust MORB Upper Mantle (80% diopside, 20% pyrope)
        std::vector<double> minfract_metaMORB_UM = {0.80,0.20}; 
        std::vector<unsigned int> minindex_metaMORB_UM = {cpxdiopsid_index, grtpyropes_index};
        // Meta-basaltic crust MORB Upper Mantle Transition Zone (50% majorite, 4% stishovite, 46% diopside)
        std::vector<double> minfract_metaMORB_UMTZ = {0.50,0.04,0.46}; 
        std::vector<unsigned int> minindex_metaMORB_UMTZ = {grtmajorit_index, stisho05Al_index, cpxdiopsid_index};
        // Meta-basaltic crust MORB Lower Mantle Transition Zone (92% majorite, 8% stishovite)
        std::vector<double> minfract_metaMORB_LMTZ = {0.92,0.08}; 
        std::vector<unsigned int> minindex_metaMORB_LMTZ = {grtmajorit_index, stisho05Al_index};
        // Meta-basaltic crust MORB Lower Mantle (35% bridgmanite, 28% davemaoite, 19% Fe-NAL, 18% stishovite)
        std::vector<double> minfract_metaMORB_LM = {0.35,0.28,0.19,0.18}; 
        std::vector<unsigned int> minindex_metaMORB_LM = {brigmaFeAl_index, davemaoite_index, newhexAlph_index, stisho05Al_index};

        // Dunite Upper Mantle (100% olivine)
        std::vector<double> minfract_duniteOl_UM = {1.00}; 
        std::vector<unsigned int> minindex_duniteOl_UM = {olivinedry_index};
        // Dunite Upper Mantle Transition Zone (100% wadsleyite)
        std::vector<double> minfract_duniteOl_UMTZ = {1.00};
        std::vector<unsigned int> minindex_duniteOl_UMTZ = {wadsleydry_index};
        // Dunite Lower Mantle Transition Zone (100% ringwoodite)
        std::vector<double> minfract_duniteOl_LMTZ = {1.00};
        std::vector<unsigned int> minindex_duniteOl_LMTZ = {ringwoodry_index};
        // Dunite Lower Mantle (100% bridgmanite)
        std::vector<double> minfract_duniteOl_LM = {1.00};
        std::vector<unsigned int> minindex_duniteOl_LM = {brigma90Mg_index};

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
        double sum_min_fract_pyrolite_UM = std::accumulate(minfract_pyrolite_UM.begin(), minfract_pyrolite_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_pyrolite_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_pyrolite_UM must be equal to 1."));
        double sum_min_fract_pyrolite_UMTZ = std::accumulate(minfract_pyrolite_UMTZ.begin(), minfract_pyrolite_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_pyrolite_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_pyrolite_UMTZ must be equal to 1."));
        double sum_min_fract_pyrolite_LMTZ = std::accumulate(minfract_pyrolite_LMTZ.begin(), minfract_pyrolite_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_pyrolite_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_pyrolite_LMTZ must be equal to 1."));
        double sum_min_fract_pyrolite_LM = std::accumulate(minfract_pyrolite_LM.begin(), minfract_pyrolite_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_pyrolite_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_pyrolite_LM must be equal to 1."));

        double sum_min_fract_harzburg_UM = std::accumulate(minfract_harzburg_UM.begin(), minfract_harzburg_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_harzburg_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_harzburg_UM must be equal to 1."));
        double sum_min_fract_harzburg_UMTZ = std::accumulate(minfract_harzburg_UMTZ.begin(), minfract_harzburg_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_harzburg_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_harzburg_UMTZ must be equal to 1."));
        double sum_min_fract_harzburg_LMTZ = std::accumulate(minfract_harzburg_LMTZ.begin(), minfract_harzburg_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_harzburg_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_harzburg_LMTZ must be equal to 1."));
        double sum_min_fract_harzburg_LM = std::accumulate(minfract_harzburg_LM.begin(), minfract_harzburg_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_harzburg_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_harzburg_LM must be equal to 1."));

        double sum_min_fract_metaMORB_UM = std::accumulate(minfract_metaMORB_UM.begin(), minfract_metaMORB_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_metaMORB_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_metaMORB_UM must be equal to 1."));
        double sum_min_fract_metaMORB_UMTZ = std::accumulate(minfract_metaMORB_UMTZ.begin(), minfract_metaMORB_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_metaMORB_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_metaMORB_UMTZ must be equal to 1."));
        double sum_min_fract_metaMORB_LMTZ = std::accumulate(minfract_metaMORB_LMTZ.begin(), minfract_metaMORB_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_metaMORB_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_metaMORB_LMTZ must be equal to 1."));
        double sum_min_fract_metaMORB_LM = std::accumulate(minfract_metaMORB_LM.begin(), minfract_metaMORB_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_metaMORB_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_metaMORB_LM must be equal to 1."));

        double sum_min_fract_duniteOl_UM = std::accumulate(minfract_duniteOl_UM.begin(), minfract_duniteOl_UM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_duniteOl_UM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_duniteOl_UM must be equal to 1."));
        double sum_min_fract_duniteOl_UMTZ = std::accumulate(minfract_duniteOl_UMTZ.begin(), minfract_duniteOl_UMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_duniteOl_UMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_duniteOl_UMTZ must be equal to 1."));
        double sum_min_fract_duniteOl_LMTZ = std::accumulate(minfract_duniteOl_LMTZ.begin(), minfract_duniteOl_LMTZ.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_duniteOl_LMTZ - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_duniteOl_LMTZ must be equal to 1."));
        double sum_min_fract_duniteOl_LM = std::accumulate(minfract_duniteOl_LM.begin(), minfract_duniteOl_LM.end(), 0.0);
        AssertThrow(std::abs(sum_min_fract_duniteOl_LM - 1.0) < 1e-6,
                    dealii::ExcMessage("Error: The sum of minfract_duniteOl_LM must be equal to 1."));

        // Define room temperature [K] 
        const double T_room = 298.15; 

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;

          // Compute natural logarithm of pressure and temperature
          AssertThrow(P_GPa > 0, dealii::ExcMessage("Pressure must be > 0 for log.")); 
          AssertThrow(in.temperature[i] > 0, dealii::ExcMessage("Temperature must be > 0 for log."));
          double P_log = std::log(P_GPa);
          double T_log = std::log(in.temperature[i]);

          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];

          // Take lithology of the model
          double lithology = in.composition[0][i];

          std::vector<double> mineral_fraction;    // Mineral fractions for the current lithology
          std::vector<unsigned int> mineral_index; // Mineral indexes for the current lithology

          if (P_GPa < 13.59280101) // upper mantle
          {
            if (lithology == 0) // pyrolite
            {
              mineral_fraction = minfract_pyrolite_UM; 
              mineral_index = minindex_pyrolite_UM;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_fraction = minfract_harzburg_UM; 
              mineral_index = minindex_harzburg_UM;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_fraction = minfract_metaMORB_UM; 
              mineral_index = minindex_metaMORB_UM;
            }
            else if (lithology == 3) // dunite
            {
              mineral_fraction = minfract_duniteOl_UM; 
              mineral_index = minindex_duniteOl_UM;
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
          else if (P_GPa >= 13.59280101 && P_GPa <= 17.69264984) // upper transition zone
          {
            if (lithology == 0) // pyrolite
            {
              mineral_fraction = minfract_pyrolite_UMTZ; 
              mineral_index = minindex_pyrolite_UMTZ;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_fraction = minfract_harzburg_UMTZ; 
              mineral_index = minindex_harzburg_UMTZ;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_fraction = minfract_metaMORB_UMTZ; 
              mineral_index = minindex_metaMORB_UMTZ;
            }
            else if (lithology == 3) // dunite
            {
              mineral_fraction = minfract_duniteOl_UMTZ; 
              mineral_index = minindex_duniteOl_UMTZ;
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
          else if (P_GPa > 17.6926498400000 && P_GPa <= 23.1122152) // lower transition zone
          {
            if (lithology == 0) // pyrolite
            {
              mineral_fraction = minfract_pyrolite_LMTZ; 
              mineral_index = minindex_pyrolite_LMTZ;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_fraction = minfract_harzburg_LMTZ; 
              mineral_index = minindex_harzburg_LMTZ;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_fraction = minfract_metaMORB_LMTZ; 
              mineral_index = minindex_metaMORB_LMTZ;
            }
            else if (lithology == 3) // dunite
            {
              mineral_fraction = minfract_duniteOl_LMTZ; 
              mineral_index = minindex_duniteOl_LMTZ;
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
          else if (P_GPa > 23.1122152) // lower mantle
          {
            if (lithology == 0) // pyrolite
            {
              mineral_fraction = minfract_pyrolite_LM; 
              mineral_index = minindex_pyrolite_LM;
            }
            else if (lithology == 1) // harzburgite
            {
              mineral_fraction = minfract_harzburg_LM; 
              mineral_index = minindex_harzburg_LM;
            }
            else if (lithology == 2) // meta-MORB
            {
              mineral_fraction = minfract_metaMORB_LM; 
              mineral_index = minindex_metaMORB_LM;
            }
            else if (lithology == 3) // dunite
            {
              mineral_fraction = minfract_duniteOl_LM; 
              mineral_index = minindex_duniteOl_LM;
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

          // unsigned int mineral_ID_nonused = in.Mineral_ID;

          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> mar25_minerals_latTcond(mineral_fraction.size(), 0.0); // Lattice thermal conductivity
          std::vector<double> mar25_minerals_radTcond(mineral_fraction.size(), 0.0); // Radiative thermal conductivity
          std::vector<double> mar25_minerals_totTcond(mineral_fraction.size(), 0.0); // Total thermal conductivity

          // Preallocate total thermal conductivity of the aggregate rock
          double mar25_aggregate_rock_totTcond = 1;

          for (size_t col = 0; col < mineral_fraction.size(); ++col)
          {

           unsigned int mID = mineral_index[col];

           switch (mID) // Compute the lattice, radiative and total thermal conductivities of the given mineral
           {
             case olivinedry_index: // Dry Olivine
             {      
               double olivinedry_latTCon = compute_lattice_thermal_conductivity_mar2025(
               olivinedry_latTC_a0, olivinedry_latTC_b1, olivinedry_latTC_ymin, olivinedry_latTC_ymax,
               P_log, T_mod, T_room, olivinedry_Tdep_n_exp); 
               double olivinedry_radTCon = compute_radiative_thermal_conductivity_mar2025(
               olivinedry_radTC_c0, olivinedry_radTC_d1, olivinedry_radTC_jmin, olivinedry_radTC_jmax, T_log); 
               double olivinedry_TotTCon = compute_total_thermal_conductivity_mar2025(
               olivinedry_latTCon, olivinedry_radTCon); 
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = olivinedry_latTCon;
               mar25_minerals_radTcond[col] = olivinedry_radTCon;
               mar25_minerals_totTcond[col] = olivinedry_TotTCon;
               break;
              }
             case wadsleydry_index: // Dry Wadsleyite 
             { 
               double wadsleydry_latTCon = compute_lattice_thermal_conductivity_mar2025(
               wadsleydry_latTC_a0, wadsleydry_latTC_b1, wadsleydry_latTC_ymin, wadsleydry_latTC_ymax,
               P_log, T_mod, T_room, wadsleydry_Tdep_n_exp);
               double wadsleydry_radTCon = compute_radiative_thermal_conductivity_mar2025(
               wadsleydry_radTC_c0, wadsleydry_radTC_d1, wadsleydry_radTC_jmin, wadsleydry_radTC_jmax, T_log);
               double wadsleydry_TotTCon = compute_total_thermal_conductivity_mar2025(
               wadsleydry_latTCon, wadsleydry_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = wadsleydry_latTCon;
               mar25_minerals_radTcond[col] = wadsleydry_radTCon;       
               mar25_minerals_totTcond[col] = wadsleydry_TotTCon;
               break;
              }
             case ringwoodry_index: // Dry Ringwoodite
             { 
               double ringwoodry_latTCon = compute_lattice_thermal_conductivity_mar2025(
               ringwoodry_latTC_a0, ringwoodry_latTC_b1, ringwoodry_latTC_ymin, ringwoodry_latTC_ymax,
               P_log, T_mod, T_room, ringwoodry_Tdep_n_exp);
               double ringwoodry_radTCon = compute_radiative_thermal_conductivity_mar2025(
               ringwoodry_radTC_c0, ringwoodry_radTC_d1, ringwoodry_radTC_jmin, ringwoodry_radTC_jmax, T_log);    
               double ringwoodry_TotTCon = compute_total_thermal_conductivity_mar2025(
               ringwoodry_latTCon, ringwoodry_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ringwoodry_latTCon;
               mar25_minerals_radTcond[col] = ringwoodry_radTCon;
               mar25_minerals_totTcond[col] = ringwoodry_TotTCon;
               break;
              }
             case brigm100Mg_index: // Mg-Bridgmanite
             { 
               double brigm100Mg_latTCon = compute_lattice_thermal_conductivity_mar2025(
               brigm100Mg_latTC_a0, brigm100Mg_latTC_b1, brigm100Mg_latTC_ymin, brigm100Mg_latTC_ymax,
               P_log, T_mod, T_room, brigm100Mg_Tdep_n_exp);
               double brigm100Mg_radTCon = compute_radiative_thermal_conductivity_mar2025(
               brigm100Mg_radTC_c0, brigm100Mg_radTC_d1, brigm100Mg_radTC_jmin, brigm100Mg_radTC_jmax, T_log);
               double brigm100Mg_TotTCon = compute_total_thermal_conductivity_mar2025(
               brigm100Mg_latTCon, brigm100Mg_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigm100Mg_latTCon;
               mar25_minerals_totTcond[col] = brigm100Mg_TotTCon;
               mar25_minerals_radTcond[col] = brigm100Mg_radTCon;
               break;
              }
             case brigma97Mg_index: // Fe-Bridgmanite (3%)
             { 
               double brigma97Mg_latTCon = compute_lattice_thermal_conductivity_mar2025(
               brigma97Mg_latTC_a0, brigma97Mg_latTC_b1, brigma97Mg_latTC_ymin, brigma97Mg_latTC_ymax,
               P_log, T_mod, T_room, brigma97Mg_Tdep_n_exp);
               double brigma97Mg_radTCon = compute_radiative_thermal_conductivity_mar2025(
               brigma97Mg_radTC_c0, brigma97Mg_radTC_d1, brigma97Mg_radTC_jmin, brigma97Mg_radTC_jmax, T_log);
               double brigma97Mg_TotTCon = compute_total_thermal_conductivity_mar2025(
               brigma97Mg_latTCon, brigma97Mg_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigma97Mg_latTCon;
               mar25_minerals_radTcond[col] = brigma97Mg_radTCon;
               mar25_minerals_totTcond[col] = brigma97Mg_TotTCon;
               break;
              }
             case brigma90Mg_index: // Fe-Bridgmanite (10%)
             { 
               double brigma90Mg_latTCon = compute_lattice_thermal_conductivity_mar2025(
               brigma90Mg_latTC_a0, brigma90Mg_latTC_b1, brigma90Mg_latTC_ymin, brigma90Mg_latTC_ymax,
               P_log, T_mod, T_room, brigma90Mg_Tdep_n_exp);
               double brigma90Mg_radTCon = compute_radiative_thermal_conductivity_mar2025(
               brigma90Mg_radTC_c0, brigma90Mg_radTC_d1, brigma90Mg_radTC_jmin, brigma90Mg_radTC_jmax, T_log);
               double brigma90Mg_TotTCon = compute_total_thermal_conductivity_mar2025(
               brigma90Mg_latTCon, brigma90Mg_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigma90Mg_latTCon;
               mar25_minerals_radTcond[col] = brigma90Mg_radTCon;
               mar25_minerals_totTcond[col] = brigma90Mg_TotTCon;
               break;
              }
             case brigmaAlMg_index: // Al-Bridgmanite
             { 
               double brigmaAlMg_latTCon = compute_lattice_thermal_conductivity_mar2025(
               brigmaAlMg_latTC_a0, brigmaAlMg_latTC_b1, brigmaAlMg_latTC_ymin, brigmaAlMg_latTC_ymax,
               P_log, T_mod, T_room, brigmaAlMg_Tdep_n_exp);
               double brigmaAlMg_radTCon = compute_radiative_thermal_conductivity_mar2025(
               brigmaAlMg_radTC_c0, brigmaAlMg_radTC_d1, brigmaAlMg_radTC_jmin, brigmaAlMg_radTC_jmax, T_log);
               double brigmaAlMg_TotTCon = compute_total_thermal_conductivity_mar2025(
               brigmaAlMg_latTCon, brigmaAlMg_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigmaAlMg_latTCon;
               mar25_minerals_radTcond[col] = brigmaAlMg_radTCon;
               mar25_minerals_totTcond[col] = brigmaAlMg_TotTCon;
               break;
              }
             case brigmaFeAl_index: // Fe,Al-Bridgmanite
             { 
               double brigmaFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025(
               brigmaFeAl_latTC_a0, brigmaFeAl_latTC_b1, brigmaFeAl_latTC_ymin, brigmaFeAl_latTC_ymax,
               P_log, T_mod, T_room, brigmaFeAl_Tdep_n_exp);
               double brigmaFeAl_radTCon = compute_radiative_thermal_conductivity_mar2025(
               brigmaFeAl_radTC_c0, brigmaFeAl_radTC_d1, brigmaFeAl_radTC_jmin, brigmaFeAl_radTC_jmax, T_log);
               double brigmaFeAl_TotTCon = compute_total_thermal_conductivity_mar2025(
               brigmaFeAl_latTCon, brigmaFeAl_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigmaFeAl_latTCon;
               mar25_minerals_radTcond[col] = brigmaFeAl_radTCon;
               mar25_minerals_totTcond[col] = brigmaFeAl_TotTCon;
               break;
              }
             case opxenstati_index: // Orthopyroxene (Enstatite)
             { 
               double opxenstati_latTCon = compute_lattice_thermal_conductivity_mar2025(
               opxenstati_latTC_a0, opxenstati_latTC_b1, opxenstati_latTC_ymin, opxenstati_latTC_ymax,
               P_log, T_mod, T_room, opxenstati_Tdep_n_exp);  
               double opxenstati_radTCon = compute_radiative_thermal_conductivity_mar2025(
               opxenstati_radTC_c0, opxenstati_radTC_d1, opxenstati_radTC_jmin, opxenstati_radTC_jmax, T_log);
               double opxenstati_TotTCon = compute_total_thermal_conductivity_mar2025(
               opxenstati_latTCon, opxenstati_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = opxenstati_latTCon;
               mar25_minerals_radTcond[col] = opxenstati_radTCon;
               mar25_minerals_totTcond[col] = opxenstati_TotTCon;
               break;
              }
             case cpxdiopsid_index: // Clinopyroxene (Diopside)
             { 
               double cpxdiopsid_latTCon = compute_lattice_thermal_conductivity_mar2025(
               cpxdiopsid_latTC_a0, cpxdiopsid_latTC_b1, cpxdiopsid_latTC_ymin, cpxdiopsid_latTC_ymax,
               P_log, T_mod, T_room, cpxdiopsid_Tdep_n_exp);      
               double cpxdiopsid_radTCon = compute_radiative_thermal_conductivity_mar2025(
               cpxdiopsid_radTC_c0, cpxdiopsid_radTC_d1, cpxdiopsid_radTC_jmin, cpxdiopsid_radTC_jmax, T_log);
               double cpxdiopsid_TotTCon = compute_total_thermal_conductivity_mar2025(
               cpxdiopsid_latTCon, cpxdiopsid_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = cpxdiopsid_latTCon;
               mar25_minerals_radTcond[col] = cpxdiopsid_radTCon;
               mar25_minerals_totTcond[col] = cpxdiopsid_TotTCon;
               break;
              }
             case grtpyropes_index: // Garnet (Pyrope)
             { 
               double grtpyropes_latTCon = compute_lattice_thermal_conductivity_mar2025(
               grtpyropes_latTC_a0, grtpyropes_latTC_b1, grtpyropes_latTC_ymin, grtpyropes_latTC_ymax,
               P_log, T_mod, T_room, grtpyropes_Tdep_n_exp);
               double grtpyropes_radTCon = compute_radiative_thermal_conductivity_mar2025(
               grtpyropes_radTC_c0, grtpyropes_radTC_d1, grtpyropes_radTC_jmin, grtpyropes_radTC_jmax, T_log);
               double grtpyropes_TotTCon = compute_total_thermal_conductivity_mar2025(
               grtpyropes_latTCon, grtpyropes_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtpyropes_latTCon;
               mar25_minerals_radTcond[col] = grtpyropes_radTCon;
               mar25_minerals_totTcond[col] = grtpyropes_TotTCon;
               break;
             }
             case grtgrossul_index: // Garnet (Grossular)
             { 
               double grtgrossul_latTCon = compute_lattice_thermal_conductivity_mar2025(
               grtgrossul_latTC_a0, grtgrossul_latTC_b1, grtgrossul_latTC_ymin, grtgrossul_latTC_ymax,
               P_log, T_mod, T_room, grtgrossul_Tdep_n_exp);
               double grtgrossul_radTCon = compute_radiative_thermal_conductivity_mar2025(
               grtgrossul_radTC_c0, grtgrossul_radTC_d1, grtgrossul_radTC_jmin, grtgrossul_radTC_jmax, T_log);
               double grtgrossul_TotTCon = compute_total_thermal_conductivity_mar2025(
               grtgrossul_latTCon, grtgrossul_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtgrossul_latTCon;
               mar25_minerals_radTcond[col] = grtgrossul_radTCon;
               mar25_minerals_totTcond[col] = grtgrossul_TotTCon;
               break;
              }
             case grtalmandi_index: // Garnet (Almandine)
             { 
               double grtalmandi_latTCon = compute_lattice_thermal_conductivity_mar2025(
               grtalmandi_latTC_a0, grtalmandi_latTC_b1, grtalmandi_latTC_ymin, grtalmandi_latTC_ymax,
               P_log, T_mod, T_room, grtalmandi_Tdep_n_exp);
               double grtalmandi_radTCon = compute_radiative_thermal_conductivity_mar2025(
               grtalmandi_radTC_c0, grtalmandi_radTC_d1, grtalmandi_radTC_jmin, grtalmandi_radTC_jmax, T_log);  
               double grtalmandi_TotTCon = compute_total_thermal_conductivity_mar2025(
               grtalmandi_latTCon, grtalmandi_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtalmandi_latTCon;
               mar25_minerals_radTcond[col] = grtalmandi_radTCon;
               mar25_minerals_totTcond[col] = grtalmandi_TotTCon;
               break;
              }
             case grtmajorit_index: // Garnet (Majorite)
             { 
               double grtmajorit_latTCon = compute_lattice_thermal_conductivity_mar2025(
               grtmajorit_latTC_a0, grtmajorit_latTC_b1, grtmajorit_latTC_ymin, grtmajorit_latTC_ymax,
               P_log, T_mod, T_room, grtmajorit_Tdep_n_exp);
               double grtmajorit_radTCon = compute_radiative_thermal_conductivity_mar2025(
               grtmajorit_radTC_c0, grtmajorit_radTC_d1, grtmajorit_radTC_jmin, grtmajorit_radTC_jmax, T_log);
               double grtmajorit_TotTCon = compute_total_thermal_conductivity_mar2025(
               grtmajorit_latTCon, grtmajorit_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtmajorit_latTCon;
               mar25_minerals_radTcond[col] = grtmajorit_radTCon;
               mar25_minerals_totTcond[col] = grtmajorit_TotTCon;
               break;
              }
             case quartzpure_index: // Quartz
             { 
               double quartzpure_latTCon = compute_lattice_thermal_conductivity_mar2025(
               quartzpure_latTC_a0, quartzpure_latTC_b1, quartzpure_latTC_ymin, quartzpure_latTC_ymax,
               P_log, T_mod, T_room, quartzpure_Tdep_n_exp);
               double quartzpure_radTCon = compute_radiative_thermal_conductivity_mar2025(
               quartzpure_radTC_c0, quartzpure_radTC_d1, quartzpure_radTC_jmin, quartzpure_radTC_jmax, T_log);
               double quartzpure_TotTCon = compute_total_thermal_conductivity_mar2025(
               quartzpure_latTCon, quartzpure_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = quartzpure_latTCon;
               mar25_minerals_radTcond[col] = quartzpure_radTCon;
               mar25_minerals_totTcond[col] = quartzpure_TotTCon;
               break;
              }
             case coesitSiO2_index: // Coesite
             { 
               double coesitSiO2_latTCon = compute_lattice_thermal_conductivity_mar2025(
               coesitSiO2_latTC_a0, coesitSiO2_latTC_b1, coesitSiO2_latTC_ymin, coesitSiO2_latTC_ymax,
               P_log, T_mod, T_room, coesitSiO2_Tdep_n_exp);
               double coesitSiO2_radTCon = compute_radiative_thermal_conductivity_mar2025(
               coesitSiO2_radTC_c0, coesitSiO2_radTC_d1, coesitSiO2_radTC_jmin, coesitSiO2_radTC_jmax, T_log);
               double coesitSiO2_TotTCon = compute_total_thermal_conductivity_mar2025(
               coesitSiO2_latTCon, coesitSiO2_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = coesitSiO2_latTCon;
               mar25_minerals_radTcond[col] = coesitSiO2_radTCon;
               mar25_minerals_totTcond[col] = coesitSiO2_TotTCon;
               break;
              }
             case stishovite_index: // stishovite
             { 
               double stishovite_latTCon; // Declare the variable
               if (P_GPa < 52) // Pressure < 52 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 stishovite_latTC_1_a0, stishovite_latTC_1_b1, stishovite_latTC_1_ymin, stishovite_latTC_1_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else if (P_GPa >= 52 && P_GPa <= 56) // Pressure between 52 and 56 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 stishovite_latTC_2_a0, stishovite_latTC_2_b1, stishovite_latTC_2_ymin, stishovite_latTC_2_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else if (P_GPa > 56) // Pressure > 56 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 stishovite_latTC_3_a0, stishovite_latTC_3_b1, stishovite_latTC_3_ymin, stishovite_latTC_3_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else
               {
                 AssertThrow(false, dealii::ExcMessage("Invalid pressure value for stishovite_latTC coefficients."));
               } 
               double stishovite_radTCon = compute_radiative_thermal_conductivity_mar2025(
               stishovite_radTC_c0, stishovite_radTC_d1, stishovite_radTC_jmin, stishovite_radTC_jmax, T_log);
               double stishovite_TotTCon = compute_total_thermal_conductivity_mar2025(
               stishovite_latTCon, stishovite_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = stishovite_latTCon;
               mar25_minerals_radTcond[col] = stishovite_radTCon;
               mar25_minerals_totTcond[col] = stishovite_TotTCon;
               break;
             }
             case stisho05Al_index: // Al-stishovite (5 vol%)
             { 
               double stisho05Al_latTCon = compute_lattice_thermal_conductivity_mar2025(
               stisho05Al_latTC_a0, stisho05Al_latTC_b1, stisho05Al_latTC_ymin, stisho05Al_latTC_ymax,
               P_log, T_mod, T_room, stisho05Al_Tdep_n_exp); 
               double stisho05Al_radTCon = compute_radiative_thermal_conductivity_mar2025(
               stisho05Al_radTC_c0, stisho05Al_radTC_d1, stisho05Al_radTC_jmin, stisho05Al_radTC_jmax, T_log);
               double stisho05Al_TotTCon = compute_total_thermal_conductivity_mar2025(
               stisho05Al_latTCon, stisho05Al_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = stisho05Al_latTCon;
               mar25_minerals_radTcond[col] = stisho05Al_radTCon;
               mar25_minerals_totTcond[col] = stisho05Al_TotTCon;
               break;
             }
             case antigor010_index: // Antigorite (010)
             { 
               double antigor010_latTCon = compute_lattice_thermal_conductivity_mar2025(
               antigor010_latTC_a0, antigor010_latTC_b1, antigor010_latTC_ymin, antigor010_latTC_ymax,
               P_log, T_mod, T_room, antigor010_Tdep_n_exp);
               double antigor010_radTCon = compute_radiative_thermal_conductivity_mar2025(
               antigor010_radTC_c0, antigor010_radTC_d1, antigor010_radTC_jmin, antigor010_radTC_jmax, T_log);
               double antigor010_TotTCon = compute_total_thermal_conductivity_mar2025(
               antigor010_latTCon, antigor010_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = antigor010_latTCon;
               mar25_minerals_radTcond[col] = antigor010_radTCon;
               mar25_minerals_totTcond[col] = antigor010_TotTCon;
               break;
             }
             case antigor001_index: // Antigorite (001)
             { 
               double antigor001_latTCon = compute_lattice_thermal_conductivity_mar2025(
               antigor001_latTC_a0, antigor001_latTC_b1, antigor001_latTC_ymin, antigor001_latTC_ymax,
               P_log, T_mod, T_room, antigor001_Tdep_n_exp);
               double antigor001_radTCon = compute_radiative_thermal_conductivity_mar2025(
               antigor001_radTC_c0, antigor001_radTC_d1, antigor001_radTC_jmin, antigor001_radTC_jmax, T_log);
               double antigor001_TotTCon = compute_total_thermal_conductivity_mar2025(
               antigor001_latTCon, antigor001_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = antigor001_latTCon;
               mar25_minerals_radTcond[col] = antigor001_radTCon;
               mar25_minerals_totTcond[col] = antigor001_TotTCon;
               break;
              }
             case phaseDFeAl_index: // Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
             { 
               double phaseDFeAl_latTCon; // Declare the variable
               if (P_GPa < 24) // Pressure < 24 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 phaseDFeAl_latTC_1_a0, phaseDFeAl_latTC_1_b1, phaseDFeAl_latTC_1_ymin, phaseDFeAl_latTC_1_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa >= 24 && P_GPa <= 38) // Pressure between 24 and 38 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 phaseDFeAl_latTC_2_a0, phaseDFeAl_latTC_2_b1, phaseDFeAl_latTC_2_ymin, phaseDFeAl_latTC_2_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa >= 38 && P_GPa <= 48) // Pressure between 38 and 48 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 phaseDFeAl_latTC_3_a0, phaseDFeAl_latTC_3_b1, phaseDFeAl_latTC_3_ymin, phaseDFeAl_latTC_3_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa > 48) // Pressure > 48 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025(
                 phaseDFeAl_latTC_4_a0, phaseDFeAl_latTC_4_b1, phaseDFeAl_latTC_4_ymin, phaseDFeAl_latTC_4_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else
               {
                 AssertThrow(false, dealii::ExcMessage("Invalid pressure value for phaseDFeAl_latTC coefficients."));
               } 
               double phaseDFeAl_radTCon = compute_radiative_thermal_conductivity_mar2025(
               phaseDFeAl_radTC_c0, phaseDFeAl_radTC_d1, phaseDFeAl_radTC_jmin, phaseDFeAl_radTC_jmax, T_log);
               double phaseDFeAl_TotTCon = compute_total_thermal_conductivity_mar2025(
               phaseDFeAl_latTCon, phaseDFeAl_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = phaseDFeAl_latTCon;
               mar25_minerals_radTcond[col] = phaseDFeAl_radTCon;
               mar25_minerals_totTcond[col] = phaseDFeAl_TotTCon;
               break;
              }
             case phaseD02Al_index: // Al-phase D (Dense Hydrous Magnesium Silicate)
             { 
               double phaseD02Al_latTCon = compute_lattice_thermal_conductivity_mar2025(
               phaseD02Al_latTC_a0, phaseD02Al_latTC_b1, phaseD02Al_latTC_ymin, phaseD02Al_latTC_ymax,
               P_log, T_mod, T_room, phaseD02Al_Tdep_n_exp);
               double phaseD02Al_radTCon = compute_radiative_thermal_conductivity_mar2025(
               phaseD02Al_radTC_c0, phaseD02Al_radTC_d1, phaseD02Al_radTC_jmin, phaseD02Al_radTC_jmax, T_log);
               double phaseD02Al_TotTCon = compute_total_thermal_conductivity_mar2025(
               phaseD02Al_latTCon, phaseD02Al_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = phaseD02Al_latTCon;
               mar25_minerals_radTcond[col] = phaseD02Al_radTCon;
               mar25_minerals_totTcond[col] = phaseD02Al_TotTCon;
               break;
              }
             case ferroper08_index: // Ferropericlase (Mg92Fe8O)
             { 
               double ferroper08_latTCon = compute_lattice_thermal_conductivity_mar2025(
               ferroper08_latTC_a0, ferroper08_latTC_b1, ferroper08_latTC_ymin, ferroper08_latTC_ymax,
               P_log, T_mod, T_room, ferroper08_Tdep_n_exp);
               double ferroper08_radTCon = compute_radiative_thermal_conductivity_mar2025(
               ferroper08_radTC_c0, ferroper08_radTC_d1, ferroper08_radTC_jmin, ferroper08_radTC_jmax, T_log);
               double ferroper08_TotTCon = compute_total_thermal_conductivity_mar2025(
               ferroper08_latTCon, ferroper08_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper08_latTCon;
               mar25_minerals_radTcond[col] = ferroper08_radTCon;
               mar25_minerals_totTcond[col] = ferroper08_TotTCon;
               break;
              }
             case ferroper10_index: // Ferropericlase (Mg90Fe10O)
             { 
               double ferroper10_latTCon = compute_lattice_thermal_conductivity_mar2025(
               ferroper10_latTC_a0, ferroper10_latTC_b1, ferroper10_latTC_ymin, ferroper10_latTC_ymax,
               P_log, T_mod, T_room, ferroper10_Tdep_n_exp);
               double ferroper10_radTCon = compute_radiative_thermal_conductivity_mar2025(
               ferroper10_radTC_c0, ferroper10_radTC_d1, ferroper10_radTC_jmin, ferroper10_radTC_jmax, T_log);        
               double ferroper10_TotTCon = compute_total_thermal_conductivity_mar2025(
               ferroper10_latTCon, ferroper10_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper10_latTCon;
               mar25_minerals_radTcond[col] = ferroper10_radTCon;
               mar25_minerals_totTcond[col] = ferroper10_TotTCon;
               break;
             }
             case ferroper20_index: // Ferropericlase (Mg80Fe20O)
             { 
               double ferroper20_latTCon = compute_lattice_thermal_conductivity_mar2025(
               ferroper20_latTC_a0, ferroper20_latTC_b1, ferroper20_latTC_ymin, ferroper20_latTC_ymax,
               P_log, T_mod, T_room, ferroper20_Tdep_n_exp);
               double ferroper20_radTCon = compute_radiative_thermal_conductivity_mar2025(
               ferroper20_radTC_c0, ferroper20_radTC_d1, ferroper20_radTC_jmin, ferroper20_radTC_jmax, T_log);
               double ferroper20_TotTCon = compute_total_thermal_conductivity_mar2025(
               ferroper20_latTCon, ferroper20_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper20_latTCon;
               mar25_minerals_radTcond[col] = ferroper20_radTCon;
               mar25_minerals_totTcond[col] = ferroper20_TotTCon;
               break;
             }
             case ferroper56_index: // Ferropericlase (Mg56Fe44O)
             { 
               double ferroper56_latTCon = compute_lattice_thermal_conductivity_mar2025(
               ferroper56_latTC_a0, ferroper56_latTC_b1, ferroper56_latTC_ymin, ferroper56_latTC_ymax,
               P_log, T_mod, T_room, ferroper56_Tdep_n_exp);
               double ferroper56_radTCon = compute_radiative_thermal_conductivity_mar2025(
               ferroper56_radTC_c0, ferroper56_radTC_d1, ferroper56_radTC_jmin, ferroper56_radTC_jmax, T_log);
               double ferroper56_TotTCon = compute_total_thermal_conductivity_mar2025(
               ferroper56_latTCon, ferroper56_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper56_latTCon;
               mar25_minerals_radTcond[col] = ferroper56_radTCon;
               mar25_minerals_totTcond[col] = ferroper56_TotTCon;
               break;
             }
             case davemaoite_index: // davemaoite
             { 
               double davemaoite_latTCon = compute_lattice_thermal_conductivity_mar2025(
               davemaoite_latTC_a0, davemaoite_latTC_b1, davemaoite_latTC_ymin, davemaoite_latTC_ymax,
               P_log, T_mod, T_room, davemaoite_Tdep_n_exp);
               double davemaoite_radTCon = compute_radiative_thermal_conductivity_mar2025(
               davemaoite_radTC_c0, davemaoite_radTC_d1, davemaoite_radTC_jmin, davemaoite_radTC_jmax, T_log);
               double davemaoite_TotTCon = compute_total_thermal_conductivity_mar2025(
               davemaoite_latTCon, davemaoite_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = davemaoite_latTCon;
               mar25_minerals_radTcond[col] = davemaoite_radTCon;
               mar25_minerals_totTcond[col] = davemaoite_TotTCon;
               break;
             }
             case newhexAlph_index: // New-hexagonal-alluminium-phase (FeNAL)
             { 
               double newhexAlph_latTCon = compute_lattice_thermal_conductivity_mar2025(
               newhexAlph_latTC_a0, newhexAlph_latTC_b1, newhexAlph_latTC_ymin, newhexAlph_latTC_ymax,
               P_log, T_mod, T_room, newhexAlph_Tdep_n_exp);
               double newhexAlph_radTCon = compute_radiative_thermal_conductivity_mar2025(
               newhexAlph_radTC_c0, newhexAlph_radTC_d1, newhexAlph_radTC_jmin, newhexAlph_radTC_jmax, T_log);   
               double newhexAlph_TotTCon = compute_total_thermal_conductivity_mar2025(
               newhexAlph_latTCon, newhexAlph_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = newhexAlph_latTCon;
               mar25_minerals_radTcond[col] = newhexAlph_radTCon; 
               mar25_minerals_totTcond[col] = newhexAlph_TotTCon;
               break;
             }
             case akimotoite_index: // akimotoite
             { 
               double akimotoite_latTCon = compute_lattice_thermal_conductivity_mar2025(
               akimotoite_latTC_a0, akimotoite_latTC_b1, akimotoite_latTC_ymin, akimotoite_latTC_ymax,
               P_log, T_mod, T_room, akimotoite_Tdep_n_exp);
               double akimotoite_radTCon = compute_radiative_thermal_conductivity_mar2025(
               akimotoite_radTC_c0, akimotoite_radTC_d1, akimotoite_radTC_jmin, akimotoite_radTC_jmax, T_log);   
               double akimotoite_TotTCon = compute_total_thermal_conductivity_mar2025(
               akimotoite_latTCon, akimotoite_radTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = akimotoite_latTCon;
               mar25_minerals_radTcond[col] = akimotoite_radTCon;
               mar25_minerals_totTcond[col] = akimotoite_TotTCon;
               break;
             }
           }

           // Thermal conductivity of the aggregate rock is computed as the
           // geometric mean of the total thermal conductivities of the minerals weighted by their fraction
           mar25_aggregate_rock_totTcond = mar25_aggregate_rock_totTcond * std::pow(mar25_minerals_totTcond[col], mineral_fraction[col]);
          
          }

          if (lithology != 99)
          {
             out.thermal_conductivities[i] = mar25_aggregate_rock_totTcond;
          }
          else if (lithology == 99)
             out.thermal_conductivities[i] = mar25_minerals_totTcond[0];
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
  template class marzotto_2025<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}

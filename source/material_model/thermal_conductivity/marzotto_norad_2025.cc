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

#include <aspect/material_model/thermal_conductivity/marzotto_norad_2025.h>
#include <deal.II/base/exceptions.h> // Ensure this is included for AssertThrow

// Helper functions in anonymous namespace to compute thermal conductivities using the Marzotto et al. (2025) formulations
namespace 
{
  // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) space considering the boundaries (ymin) and (ymax)
  double compute_lattice_thermal_conductivity_mar2025_norad(double a0, double b1, double ymin, double ymax, double P_log, double T_mod, double T_room, double n_exp)
  { 
    double zsimpl = a0 + b1 * P_log;
    double ysimpl = std::exp(zsimpl);
    double yprime = ysimpl / (1 + ysimpl);
    double yreals = ymin + (ymax - ymin) * yprime;
    double lattice_thermal_conductivity = std::exp(yreals);
    return lattice_thermal_conductivity * std::pow((T_room / T_mod), n_exp);
  }
   
  double compute_total_thermal_conductivity_mar2025_norad(double lattice_thermal_conductivity)
  {
    double total_thermal_conductivity = lattice_thermal_conductivity;
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
      marzotto_norad_2025<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        // Define coefficients for lattice thermal conductivity for different minerals 

        // Coefficients for dry olivine 
        // retreived from fitting TDTR dataset of
        // [Chang et al., 2017, PNAS, vol 114, p. 4078-4081]
        // https://doi.org/10.1073/pnas.1616216114
        // mineral chemical formula [Mg1.8 Fe0.2 SiO4]    
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
        // mineral chemical formula [(Mg1.8Fe0.2)SiO4]
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
        // mineral chemical formula [(Mg1.79Fe0.17)Si1.02O4]
        constexpr int ringwoodry_index = 2;
        const double ringwoodry_latTC_a0 =  -4.5433218986644; 
        const double ringwoodry_latTC_b1 =   0.95165619669065;
        const double ringwoodry_latTC_ymin = 2.1632750590560;
        const double ringwoodry_latTC_ymax = 5.7037824746562; 
        const double ringwoodry_Tdep_n_exp = 0.508; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [MgSiO3]
        constexpr int brigm100Mg_index = 3;
        const double brigm100Mg_latTC_a0 =  -4.3687;
        const double brigm100Mg_latTC_b1 =   1.0766; 
        const double brigm100Mg_latTC_ymin = 2.376025820; 
        const double brigm100Mg_latTC_ymax = 5.010635294;  
        const double brigm100Mg_Tdep_n_exp = 1.01000;

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [Fe0.03Mg0.97SiO3]
        constexpr int brigma97Mg_index = 4;
        const double brigma97Mg_latTC_a0 =  -4.520600000;
        const double brigma97Mg_latTC_b1 =   1.019900000; 
        const double brigma97Mg_latTC_ymin = 1.750524121; 
        const double brigma97Mg_latTC_ymax = 4.499809670;  
        const double brigma97Mg_Tdep_n_exp = 0.56605;

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [Fe0.1Mg0.9SiO3]
        constexpr int brigma90Mg_index = 5;
        const double brigma90Mg_latTC_a0 =  -4.883100000;
        const double brigma90Mg_latTC_b1 =   0.980900000; 
        const double brigma90Mg_latTC_ymin = 1.333739493; 
        const double brigma90Mg_latTC_ymax = 4.382026635;  
        const double brigma90Mg_Tdep_n_exp = 0.17054;

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [(Al,Mg)SiO3]
        constexpr int brigmaAlMg_index = 6;
        const double brigmaAlMg_latTC_a0 =  -4.331500000;
        const double brigmaAlMg_latTC_b1 =   1.027000000; 
        const double brigmaAlMg_latTC_ymin = 1.845020046; 
        const double brigmaAlMg_latTC_ymax = 4.605170186;  
        const double brigmaAlMg_Tdep_n_exp = 0.61983;

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [(Fe,Al,Mg)SiO3]
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
        // mineral chemical formula [Mg2Si2O6]
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
        // mineral chemical formula [CaMgSi2O6]
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
        // mineral chemical formula [Mg3Al2Si3O12]
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
        // mineral chemical formula [(Ca0.986Fe0.014)3Al2(SiO4)3]
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
        // mineral chemical formula [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3]
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
        // mineral chemical formula [Mg3(MgSi)(SiO4)3]
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
        // mineral chemical formula [SiO2]
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
        // mineral chemical formula [SiO2]
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
        // mineral chemical formula [SiO2]
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
        // mineral chemical formula [(Al,Si)O2]
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
        // mineral chemical formula [(Mg2.80Fe0.05)Si2.08O5(OH)3.77]
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
        // mineral chemical formula [Mg1.19Fe0.12Al0.174Si1.71H2.02O6]
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
        // mineral chemical formula [Mg1.29Al0.17Si1.73H1.98O6]
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
        // mineral chemical formula [Mg0.92Fe0.08O] - (8% Iron)
        constexpr int ferroper08_index = 22;
        const double ferroper08_latTC_a0 = -6.9942;
        const double ferroper08_latTC_b1 =  1.953;
        const double ferroper08_latTC_ymin = 1.629240539; 
        const double ferroper08_latTC_ymax = 4.118362306;
        const double ferroper08_Tdep_n_exp = 0.5;
        // mineral chemical formula [Mg0.90Fe0.10O] - (10% Iron)
        constexpr int ferroper10_index = 23;
        const double ferroper10_latTC_a0 = -7.0133;
        const double ferroper10_latTC_b1 =  1.9321;
        const double ferroper10_latTC_ymin = 1.5040773968; 
        const double ferroper10_latTC_ymax = 4.0250359042;
        const double ferroper10_Tdep_n_exp = 0.5;
        // mineral chemical formula [Mg0.80Fe0.20O] (20% Iron)
        constexpr int ferroper20_index = 24;
        const double ferroper20_latTC_a0 = -5.2408;
        const double ferroper20_latTC_b1 =  0.9649;
        const double ferroper20_latTC_ymin = 1.2490430868; 
        const double ferroper20_latTC_ymax = 3.9318256327;
        const double ferroper20_Tdep_n_exp = 0.025;
        // mineral chemical formula [Mg0.44Fe0.56O] (56% Iron)
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
        // mineral chemical formula [CaSiO3]
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
        // mineral chemical formula [Na0.71Mg2.05Al4.62Si1.16Fe(2+)0.09Fe(3+)0.17O12]
        constexpr int newhexAlph_index = 27;
        const double newhexAlph_latTC_a0 = -29.421;
        const double newhexAlph_latTC_b1 =  7.7792;
        const double newhexAlph_latTC_ymin = 2.363551955; 
        const double newhexAlph_latTC_ymax = 3.653998874;
        const double newhexAlph_Tdep_n_exp = 0.5;

        // Coefficients for akimotoite
        // assumed to be equal to En100-Bridgmanite
        // [Zhang & Marzotto 2025, in preparation]
        // mineral chemical formula [MgSiO3]
        constexpr int akimotoite_index = 28;
        const double akimotoite_latTC_a0 =  -4.368700000;
        const double akimotoite_latTC_b1 =   1.076600000; 
        const double akimotoite_latTC_ymin = 2.376025820; 
        const double akimotoite_latTC_ymax = 5.010635294;  
        const double akimotoite_Tdep_n_exp = 1.01000;

        unsigned int mineralpar_index = akimotoite_index+1; // Number of minerals

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

          double current_temperature = in.temperature[i];
          if (in.temperature[i] <= 0.0) // Avoid log(0) or negative temperature
          {
            current_temperature = 1e-10; 
            // AssertThrow(in.temperature[i] > 0, dealii::ExcMessage("Temperature must be > 0 for log."));
          }

          double current_pressure = in.pressure[i];
          if (in.pressure[i] <= 0.0) // Avoid log(0) or negative pressure
          {
            current_pressure = 1e-10; 
            // AssertThrow(P_GPa > 0, dealii::ExcMessage("Pressure must be > 0 for log."));
          }

          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = current_pressure/1e9;

          // Compute natural logarithm of pressure and temperature
          double P_log = std::log(P_GPa);

          // Take the temperature field of the model [K]
          double T_mod = current_temperature;

          // Take lithology of the model
          // double lithology = in.composition[0][i];

          double lithology = 0.0;
          // if there is a compositional field, use the first one as indicator for lithology
          if (in.composition[i].size() > 0)
          {
            lithology = in.composition[i][0];
          }

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

          // Preallocate a vector for storing thermal conductivities of minerals
          std::vector<double> mar25_minerals_latTcond(mineral_fraction.size(), 0.0); // Lattice thermal conductivity
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
               double olivinedry_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               olivinedry_latTC_a0, olivinedry_latTC_b1, olivinedry_latTC_ymin, olivinedry_latTC_ymax,
               P_log, T_mod, T_room, olivinedry_Tdep_n_exp); 
               double olivinedry_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               olivinedry_latTCon); 
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = olivinedry_latTCon;
               mar25_minerals_totTcond[col] = olivinedry_TotTCon;
               break;
              }
             case wadsleydry_index: // Dry Wadsleyite 
             { 
               double wadsleydry_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               wadsleydry_latTC_a0, wadsleydry_latTC_b1, wadsleydry_latTC_ymin, wadsleydry_latTC_ymax,
               P_log, T_mod, T_room, wadsleydry_Tdep_n_exp);
               double wadsleydry_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               wadsleydry_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = wadsleydry_latTCon;    
               mar25_minerals_totTcond[col] = wadsleydry_TotTCon;
               break;
              }
             case ringwoodry_index: // Dry Ringwoodite
             { 
               double ringwoodry_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               ringwoodry_latTC_a0, ringwoodry_latTC_b1, ringwoodry_latTC_ymin, ringwoodry_latTC_ymax,
               P_log, T_mod, T_room, ringwoodry_Tdep_n_exp);   
               double ringwoodry_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               ringwoodry_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ringwoodry_latTCon;
               mar25_minerals_totTcond[col] = ringwoodry_TotTCon;
               break;
              }
             case brigm100Mg_index: // Mg-Bridgmanite
             { 
               double brigm100Mg_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               brigm100Mg_latTC_a0, brigm100Mg_latTC_b1, brigm100Mg_latTC_ymin, brigm100Mg_latTC_ymax,
               P_log, T_mod, T_room, brigm100Mg_Tdep_n_exp);
               double brigm100Mg_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               brigm100Mg_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigm100Mg_latTCon;
               mar25_minerals_totTcond[col] = brigm100Mg_TotTCon;
               break;
              }
             case brigma97Mg_index: // Fe-Bridgmanite (3%)
             { 
               double brigma97Mg_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               brigma97Mg_latTC_a0, brigma97Mg_latTC_b1, brigma97Mg_latTC_ymin, brigma97Mg_latTC_ymax,
               P_log, T_mod, T_room, brigma97Mg_Tdep_n_exp);
               double brigma97Mg_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               brigma97Mg_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigma97Mg_latTCon;
               mar25_minerals_totTcond[col] = brigma97Mg_TotTCon;
               break;
              }
             case brigma90Mg_index: // Fe-Bridgmanite (10%)
             { 
               double brigma90Mg_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               brigma90Mg_latTC_a0, brigma90Mg_latTC_b1, brigma90Mg_latTC_ymin, brigma90Mg_latTC_ymax,
               P_log, T_mod, T_room, brigma90Mg_Tdep_n_exp);
               double brigma90Mg_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               brigma90Mg_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigma90Mg_latTCon;
               mar25_minerals_totTcond[col] = brigma90Mg_TotTCon;
               break;
              }
             case brigmaAlMg_index: // Al-Bridgmanite
             { 
               double brigmaAlMg_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               brigmaAlMg_latTC_a0, brigmaAlMg_latTC_b1, brigmaAlMg_latTC_ymin, brigmaAlMg_latTC_ymax,
               P_log, T_mod, T_room, brigmaAlMg_Tdep_n_exp);
               double brigmaAlMg_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               brigmaAlMg_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigmaAlMg_latTCon;
               mar25_minerals_totTcond[col] = brigmaAlMg_TotTCon;
               break;
              }
             case brigmaFeAl_index: // Fe,Al-Bridgmanite
             { 
               double brigmaFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               brigmaFeAl_latTC_a0, brigmaFeAl_latTC_b1, brigmaFeAl_latTC_ymin, brigmaFeAl_latTC_ymax,
               P_log, T_mod, T_room, brigmaFeAl_Tdep_n_exp);
               double brigmaFeAl_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               brigmaFeAl_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = brigmaFeAl_latTCon;
               mar25_minerals_totTcond[col] = brigmaFeAl_TotTCon;
               break;
              }
             case opxenstati_index: // Orthopyroxene (Enstatite)
             { 
               double opxenstati_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               opxenstati_latTC_a0, opxenstati_latTC_b1, opxenstati_latTC_ymin, opxenstati_latTC_ymax,
               P_log, T_mod, T_room, opxenstati_Tdep_n_exp);  
               double opxenstati_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               opxenstati_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = opxenstati_latTCon;
               mar25_minerals_totTcond[col] = opxenstati_TotTCon;
               break;
              }
             case cpxdiopsid_index: // Clinopyroxene (Diopside)
             { 
               double cpxdiopsid_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               cpxdiopsid_latTC_a0, cpxdiopsid_latTC_b1, cpxdiopsid_latTC_ymin, cpxdiopsid_latTC_ymax,
               P_log, T_mod, T_room, cpxdiopsid_Tdep_n_exp);      
               double cpxdiopsid_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               cpxdiopsid_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = cpxdiopsid_latTCon;
               mar25_minerals_totTcond[col] = cpxdiopsid_TotTCon;
               break;
              }
             case grtpyropes_index: // Garnet (Pyrope)
             { 
               double grtpyropes_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               grtpyropes_latTC_a0, grtpyropes_latTC_b1, grtpyropes_latTC_ymin, grtpyropes_latTC_ymax,
               P_log, T_mod, T_room, grtpyropes_Tdep_n_exp);
               double grtpyropes_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               grtpyropes_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtpyropes_latTCon;
               mar25_minerals_totTcond[col] = grtpyropes_TotTCon;
               break;
             }
             case grtgrossul_index: // Garnet (Grossular)
             { 
               double grtgrossul_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               grtgrossul_latTC_a0, grtgrossul_latTC_b1, grtgrossul_latTC_ymin, grtgrossul_latTC_ymax,
               P_log, T_mod, T_room, grtgrossul_Tdep_n_exp);
               double grtgrossul_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               grtgrossul_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtgrossul_latTCon;
               mar25_minerals_totTcond[col] = grtgrossul_TotTCon;
               break;
              }
             case grtalmandi_index: // Garnet (Almandine)
             { 
               double grtalmandi_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               grtalmandi_latTC_a0, grtalmandi_latTC_b1, grtalmandi_latTC_ymin, grtalmandi_latTC_ymax,
               P_log, T_mod, T_room, grtalmandi_Tdep_n_exp);
               double grtalmandi_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               grtalmandi_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtalmandi_latTCon;
               mar25_minerals_totTcond[col] = grtalmandi_TotTCon;
               break;
              }
             case grtmajorit_index: // Garnet (Majorite)
             { 
               double grtmajorit_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               grtmajorit_latTC_a0, grtmajorit_latTC_b1, grtmajorit_latTC_ymin, grtmajorit_latTC_ymax,
               P_log, T_mod, T_room, grtmajorit_Tdep_n_exp);
               double grtmajorit_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               grtmajorit_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = grtmajorit_latTCon;
               mar25_minerals_totTcond[col] = grtmajorit_TotTCon;
               break;
              }
             case quartzpure_index: // Quartz
             { 
               double quartzpure_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               quartzpure_latTC_a0, quartzpure_latTC_b1, quartzpure_latTC_ymin, quartzpure_latTC_ymax,
               P_log, T_mod, T_room, quartzpure_Tdep_n_exp);
               double quartzpure_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               quartzpure_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = quartzpure_latTCon;
               mar25_minerals_totTcond[col] = quartzpure_TotTCon;
               break;
              }
             case coesitSiO2_index: // Coesite
             { 
               double coesitSiO2_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               coesitSiO2_latTC_a0, coesitSiO2_latTC_b1, coesitSiO2_latTC_ymin, coesitSiO2_latTC_ymax,
               P_log, T_mod, T_room, coesitSiO2_Tdep_n_exp);
               double coesitSiO2_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               coesitSiO2_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = coesitSiO2_latTCon;
               mar25_minerals_totTcond[col] = coesitSiO2_TotTCon;
               break;
              }
             case stishovite_index: // stishovite
             { 
               double stishovite_latTCon; // Declare the variable
               if (P_GPa < 52) // Pressure < 52 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 stishovite_latTC_1_a0, stishovite_latTC_1_b1, stishovite_latTC_1_ymin, stishovite_latTC_1_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else if (P_GPa >= 52 && P_GPa <= 56) // Pressure between 52 and 56 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 stishovite_latTC_2_a0, stishovite_latTC_2_b1, stishovite_latTC_2_ymin, stishovite_latTC_2_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else if (P_GPa > 56) // Pressure > 56 [GPa]
               {
                 stishovite_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 stishovite_latTC_3_a0, stishovite_latTC_3_b1, stishovite_latTC_3_ymin, stishovite_latTC_3_ymax,
                 P_log, T_mod, T_room, stishovite_Tdep_n_exp);
               }
               else
               {
                 AssertThrow(false, dealii::ExcMessage("Invalid pressure value for stishovite_latTC coefficients."));
               } 
               double stishovite_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               stishovite_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = stishovite_latTCon;
               mar25_minerals_totTcond[col] = stishovite_TotTCon;
               break;
             }
             case stisho05Al_index: // Al-stishovite (5 vol%)
             { 
               double stisho05Al_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               stisho05Al_latTC_a0, stisho05Al_latTC_b1, stisho05Al_latTC_ymin, stisho05Al_latTC_ymax,
               P_log, T_mod, T_room, stisho05Al_Tdep_n_exp); 
               double stisho05Al_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               stisho05Al_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = stisho05Al_latTCon;
               mar25_minerals_totTcond[col] = stisho05Al_TotTCon;
               break;
             }
             case antigor010_index: // Antigorite (010)
             { 
               double antigor010_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               antigor010_latTC_a0, antigor010_latTC_b1, antigor010_latTC_ymin, antigor010_latTC_ymax,
               P_log, T_mod, T_room, antigor010_Tdep_n_exp);
               double antigor010_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               antigor010_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = antigor010_latTCon;
               mar25_minerals_totTcond[col] = antigor010_TotTCon;
               break;
             }
             case antigor001_index: // Antigorite (001)
             { 
               double antigor001_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               antigor001_latTC_a0, antigor001_latTC_b1, antigor001_latTC_ymin, antigor001_latTC_ymax,
               P_log, T_mod, T_room, antigor001_Tdep_n_exp);
               double antigor001_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               antigor001_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = antigor001_latTCon;
               mar25_minerals_totTcond[col] = antigor001_TotTCon;
               break;
              }
             case phaseDFeAl_index: // Fe,Al-phase D (Dense Hydrous Magnesium Silicate)
             { 
               double phaseDFeAl_latTCon; // Declare the variable
               if (P_GPa < 24) // Pressure < 24 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 phaseDFeAl_latTC_1_a0, phaseDFeAl_latTC_1_b1, phaseDFeAl_latTC_1_ymin, phaseDFeAl_latTC_1_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa >= 24 && P_GPa <= 38) // Pressure between 24 and 38 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 phaseDFeAl_latTC_2_a0, phaseDFeAl_latTC_2_b1, phaseDFeAl_latTC_2_ymin, phaseDFeAl_latTC_2_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa >= 38 && P_GPa <= 48) // Pressure between 38 and 48 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 phaseDFeAl_latTC_3_a0, phaseDFeAl_latTC_3_b1, phaseDFeAl_latTC_3_ymin, phaseDFeAl_latTC_3_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else if (P_GPa > 48) // Pressure > 48 [GPa]
               {
                 phaseDFeAl_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
                 phaseDFeAl_latTC_4_a0, phaseDFeAl_latTC_4_b1, phaseDFeAl_latTC_4_ymin, phaseDFeAl_latTC_4_ymax,
                 P_log, T_mod, T_room, phaseDFeAl_Tdep_n_exp);
               }
               else
               {
                 AssertThrow(false, dealii::ExcMessage("Invalid pressure value for phaseDFeAl_latTC coefficients."));
               } 
               double phaseDFeAl_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               phaseDFeAl_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = phaseDFeAl_latTCon;
               mar25_minerals_totTcond[col] = phaseDFeAl_TotTCon;
               break;
              }
             case phaseD02Al_index: // Al-phase D (Dense Hydrous Magnesium Silicate)
             { 
               double phaseD02Al_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               phaseD02Al_latTC_a0, phaseD02Al_latTC_b1, phaseD02Al_latTC_ymin, phaseD02Al_latTC_ymax,
               P_log, T_mod, T_room, phaseD02Al_Tdep_n_exp);
               double phaseD02Al_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               phaseD02Al_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = phaseD02Al_latTCon;
               mar25_minerals_totTcond[col] = phaseD02Al_TotTCon;
               break;
              }
             case ferroper08_index: // Ferropericlase (Mg92Fe8O)
             { 
               double ferroper08_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               ferroper08_latTC_a0, ferroper08_latTC_b1, ferroper08_latTC_ymin, ferroper08_latTC_ymax,
               P_log, T_mod, T_room, ferroper08_Tdep_n_exp);
               double ferroper08_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               ferroper08_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper08_latTCon;
               mar25_minerals_totTcond[col] = ferroper08_TotTCon;
               break;
              }
             case ferroper10_index: // Ferropericlase (Mg90Fe10O)
             { 
               double ferroper10_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               ferroper10_latTC_a0, ferroper10_latTC_b1, ferroper10_latTC_ymin, ferroper10_latTC_ymax,
               P_log, T_mod, T_room, ferroper10_Tdep_n_exp);       
               double ferroper10_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               ferroper10_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper10_latTCon;
               mar25_minerals_totTcond[col] = ferroper10_TotTCon;
               break;
             }
             case ferroper20_index: // Ferropericlase (Mg80Fe20O)
             { 
               double ferroper20_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               ferroper20_latTC_a0, ferroper20_latTC_b1, ferroper20_latTC_ymin, ferroper20_latTC_ymax,
               P_log, T_mod, T_room, ferroper20_Tdep_n_exp);
               double ferroper20_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               ferroper20_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper20_latTCon;
               mar25_minerals_totTcond[col] = ferroper20_TotTCon;
               break;
             }
             case ferroper56_index: // Ferropericlase (Mg56Fe44O)
             { 
               double ferroper56_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               ferroper56_latTC_a0, ferroper56_latTC_b1, ferroper56_latTC_ymin, ferroper56_latTC_ymax,
               P_log, T_mod, T_room, ferroper56_Tdep_n_exp);
               double ferroper56_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               ferroper56_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = ferroper56_latTCon;
               mar25_minerals_totTcond[col] = ferroper56_TotTCon;
               break;
             }
             case davemaoite_index: // davemaoite
             { 
               double davemaoite_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               davemaoite_latTC_a0, davemaoite_latTC_b1, davemaoite_latTC_ymin, davemaoite_latTC_ymax,
               P_log, T_mod, T_room, davemaoite_Tdep_n_exp);
               double davemaoite_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               davemaoite_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = davemaoite_latTCon;
               mar25_minerals_totTcond[col] = davemaoite_TotTCon;
               break;
             }
             case newhexAlph_index: // New-hexagonal-alluminium-phase (FeNAL)
             { 
               double newhexAlph_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               newhexAlph_latTC_a0, newhexAlph_latTC_b1, newhexAlph_latTC_ymin, newhexAlph_latTC_ymax,
               P_log, T_mod, T_room, newhexAlph_Tdep_n_exp);
               double newhexAlph_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               newhexAlph_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = newhexAlph_latTCon;
               mar25_minerals_totTcond[col] = newhexAlph_TotTCon;
               break;
             }
             case akimotoite_index: // akimotoite
             { 
               double akimotoite_latTCon = compute_lattice_thermal_conductivity_mar2025_norad(
               akimotoite_latTC_a0, akimotoite_latTC_b1, akimotoite_latTC_ymin, akimotoite_latTC_ymax,
               P_log, T_mod, T_room, akimotoite_Tdep_n_exp);  
               double akimotoite_TotTCon = compute_total_thermal_conductivity_mar2025_norad(
               akimotoite_latTCon);
               // Store the thermal conductivities in the vector
               mar25_minerals_latTcond[col] = akimotoite_latTCon;
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
  template class marzotto_norad_2025<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}


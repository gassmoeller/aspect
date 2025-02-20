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
      template <int dim>
      void
      PTdepRbounded<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        // Define coefficients for lattice thermal conductivity for different minerals 

        // Coefficients for dry olivine 
        // retreived from fitting TDTR dataset of
        // [Chang et al., 2017, PNAS, vol 114, p. 4078-4081]
        // https://doi.org/10.1073/pnas.1616216114
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        const double OlivineDry_LatTC_a0 =   -4.124100000;
        const double OlivineDry_LatTC_b1 =    2.146900000;
        const double OlivineDry_LatTC_ymin =  1.28093384543429;
        const double OlivineDry_LatTC_ymax =  2.60726381956037;
        const double OlivineDry_TDep_n_Exp =  0.5;     
        
        // Coefficients for dry wadsleyite
        // retreived from fitting dataset of
        // [Xu et al., 2004, PEPI, vol 143, pp. 321-336]
        // https://doi.org/10.1016/j.pepi.2004.03.005
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        // const double WadsleyDry_LatTC_a0 =  -0.665600000;
        // const double WadsleyDry_LatTC_b1 =   0.308200000;
        // const double WadsleyDry_LatTC_ymin = 1.75735804249439;
        // const double WadsleyDry_LatTC_ymax = 2.37090451537473; 
        // const double WadsleyDry_TDep_n_Exp = 0.5; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Marzotto et al., 2020, GRL, vol 47, issue 13]
        // https://doi.org/10.1029/2020GL087607
        // mineral composition [(Mg1.79Fe0.17)Si1.02O4]
        // const double RingwooDry_LatTC_a0 =  -5.462400000;
        // const double RingwooDry_LatTC_b1 =   2.079100000;
        // const double RingwooDry_LatTC_ymin = 1.60943791241410;
        // const double RingwooDry_LatTC_ymax = 2.94939766245070; 
        // const double RingwooDry_TDep_n_Exp = 0.5; 

        // Coefficients for Mg-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [MgSiO3]
        // const double En100Brigm_LatTC_a0 =  -4.368700000;
        // const double En100Brigm_LatTC_b1 =   1.076600000; 
        // const double En100Brigm_LatTC_ymin = 2.376025820; 
        // const double En100Brigm_LatTC_ymax = 5.010635294;  
        // const double En100Brigm_TDep_n_Exp = 1.01000;

        // Coefficients for Fe-bridgmanite (3%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.03Mg0.97SiO3]
        // const double En97Brigma_LatTC_a0 =  -4.520600000;
        // const double En97Brigma_LatTC_b1 =   1.019900000; 
        // const double En97Brigma_LatTC_ymin = 1.750524121; 
        // const double En97Brigma_LatTC_ymax = 4.499809670;  
        // const double En97Brigma_TDep_n_Exp = 0.56605;

        // Coefficients for Fe-bridgmanite (10%)
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [Fe0.1Mg0.9SiO3]
        // const double En90Brigma_LatTC_a0 =  -4.883100000;
        // const double En90Brigma_LatTC_b1 =   0.980900000; 
        // const double En90Brigma_LatTC_ymin = 1.333739493; 
        // const double En90Brigma_LatTC_ymax = 4.382026635;  
        // const double En90Brigma_TDep_n_Exp = 0.17054;

        // Coefficients for Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Al,Mg)SiO3]
        // const double En90Brigma_LatTC_a0 =  -4.331500000;
        // const double En90Brigma_LatTC_b1 =   1.027000000; 
        // const double En90Brigma_LatTC_ymin = 1.845020046; 
        // const double En90Brigma_LatTC_ymax = 4.605170186;  
        // const double En90Brigma_TDep_n_Exp = 0.61983;

        // Coefficients for Fe,Al-bridgmanite
        // retreived from fitting dataset of
        // [Zhang & Marzotto 2025, in preparation]
        // mineral composition [(Fe,Al,Mg)SiO3]
        // const double En90Brigma_LatTC_a0 =  -4.510600000;
        // const double En90Brigma_LatTC_b1 =   1.066800000; 
        // const double En90Brigma_LatTC_ymin = 1.389093953; 
        // const double En90Brigma_LatTC_ymax = 3.912023005;  
        // const double En90Brigma_TDep_n_Exp = 0.46815;
         
        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of 
        // [Schloessin & Dvorak, 1972, GJI, 27(5), 499-516]
        // https://doi.org/10.1111/j.1365-246X.1972.tb06105.x
        // mineral composition [Mg2Si2O6]
        // const double OpxEnstati_LatTC_a0 =   -3.004700000;
        // const double OpxEnstati_LatTC_b1 =    2.600000000;
        // const double OpxEnstati_LatTC_ymin =  1.760865151; 
        // const double OpxEnstati_LatTC_ymax =  2.096937429;
        // const double OpxEnstati_TDep_n_Exp =  0.5;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of 
        // [Wang et al., 2014, JGR: Solid Earth, 119(8), 6277-6287]
        // https://doi.org/10.1002/2014JB011208
        // mineral composition [CaMgSi2O6]
        // const double CpxDiopsid_LatTC_a0 =   -3.251100000;
        // const double CpxDiopsid_LatTC_b1 =    1.689100000;
        // const double CpxDiopsid_LatTC_ymin =  1.793640135; 
        // const double CpxDiopsid_LatTC_ymax =  2.389462023;
        // const double CpxDiopsid_TDep_n_Exp =  0.5;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [Mg3Al2Si3O12]
        // const double GrtPyropes_LatTC_a0 =   -4.363700000;
        // const double GrtPyropes_LatTC_b1 =    2.036800000;
        // const double GrtPyropes_LatTC_ymin =  1.481604541; 
        // const double GrtPyropes_LatTC_ymax =  2.443131606;
        // const double GrtPyropes_TDep_n_Exp =  0.4314;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Ca0.986Fe0.014)3Al2(SiO4)3]
        // const double GrtGrossul_LatTC_a0 =  -4.7584;
        // const double GrtGrossul_LatTC_b1 =   2.0816;
        // const double GrtGrossul_LatTC_ymin = 1.410986974; 
        // const double GrtGrossul_LatTC_ymax = 2.457625992;
        // const double GrtGrossul_TDep_n_Exp = 0.4589;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // https://doi.org/10.2138/am-2023-8953
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3]
        // const double GrtAlmandi_LatTC_a0 =  -4.5047;
        // const double GrtAlmandi_LatTC_b1 =   2.0988;
        // const double GrtAlmandi_LatTC_ymin = 1.223775432; 
        // const double GrtAlmandi_LatTC_ymax = 2.374762159;
        // const double GrtAlmandi_TDep_n_Exp = 0.4172;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Giesting et al.2004  EPSL, 218(1-2), 45-56]
        // https://doi.org/10.1016/S0012-821X(03)00630-7
        // mineral composition [Mg3(MgSi)(SiO4)3]
        // const double GrtMajorit_LatTC_a0 =  -4.3637;
        // const double GrtMajorit_LatTC_b1 =   2.0368;
        // const double GrtMajorit_LatTC_ymin = 2.279316466; 
        // const double GrtMajorit_LatTC_ymax = 2.718047842;
        // const double GrtMajorit_TDep_n_Exp = 0.5;

        // Coefficients for quartz 
        // retreived from fitting dataset of
        // [Xiong et al., 2019 - Journal of Applied Physics, 126(21)]
        // https://doi.org/10.1063/1.5114992
        // mineral composition [SiO2]
        // const double QuartzPure_LatTC_a0 =   -2.0203;
        // const double QuartzPure_LatTC_b1 =    2.4456;
        // const double QuartzPure_LatTC_ymin =  2.260981081; 
        // const double QuartzPure_LatTC_ymax =  2.745391462;
        // const double QuartzPure_TDep_n_Exp = -1.015433333;

        // Coefficients for coesite
        // retreived from fitting dataset of
        // [Yukutake & Shimada, 1978, PEPI, 17(3), 193-200]
        // https://doi.org/10.1016/0031-9201(78)90036-5
        // mineral composition [SiO2]
        // const double CoesitSiO2_LatTC_a0 =   -12.728;
        // const double CoesitSiO2_LatTC_b1 =    2.9998;
        // const double CoesitSiO2_LatTC_ymin =  1.982022416; 
        // const double CoesitSiO2_LatTC_ymax =  2.249036030;
        // const double CoesitSiO2_TDep_n_Exp = -1.015433333;

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
        // const double WadsleyDry_RadTC_c0 =   -21.717000000;
        // const double WadsleyDry_RadTC_d1 =    3.4271000000;
        // const double WadsleyDry_RadTC_jmin = -23.025850930;
        // const double WadsleyDry_RadTC_jmax =  1.0773006810; 

        // Coefficients for dry ringwoodite
        // retreived from fitting dataset of
        // [Thomas et al., 2012, EPSL, vol. 357, p. 130-136.]
        // https://doi.org/10.1016/j.epsl.2012.09.035
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        // const double RingwooDry_RadTC_c0 =   -23.067000000;
        // const double RingwooDry_RadTC_d1 =    3.5985000000;
        // const double RingwooDry_RadTC_jmin = -23.025850930;
        // const double RingwooDry_RadTC_jmax =  0.4027750000; 

        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg2Si2O6]
        // const double OpxEnstati_RadTC_c0 =   -13.532000000;
        // const double OpxEnstati_RadTC_d1 =    2.4004000000;
        // const double OpxEnstati_RadTC_jmin = -23.025850930; 
        // const double OpxEnstati_RadTC_jmax =  1.4456685920;

        // Coefficients for clinopyroxene (diopside)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [CaMgSi2O6]
        // const double CpxDiopsid_RadTC_c0 =   -13.972000000;
        // const double CpxDiopsid_RadTC_d1 =    2.4892000000;
        // const double CpxDiopsid_RadTC_jmin = -23.025850930; 
        // const double CpxDiopsid_RadTC_jmax =  1.4446390700;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3Al2(SiO4)3]
        // const double GrtPyropes_RadTC_c0 =   -11.782000000;
        // const double GrtPyropes_RadTC_d1 =    2.0718000000;
        // const double GrtPyropes_RadTC_jmin = -23.025850930; 
        // const double GrtPyropes_RadTC_jmax =  1.4479836950;

        // Coefficients for garnet (grossular)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Ca3Al2(SiO4)3]
        // const double GrtGrossul_RadTC_c0 =   -11.261000000;
        // const double GrtGrossul_RadTC_d1 =    2.0132000000;
        // const double GrtGrossul_RadTC_jmin = -23.025850930; 
        // const double GrtGrossul_RadTC_jmax =  1.4486839760;

        // Coefficients for garnet (almandine)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [(Mg0.44Fe0.45Ca0.1Mn0.01)3Al2(SiO4)3 ]
        // const double GrtAlmandi_RadTC_c0  =  -11.261000000;
        // const double GrtAlmandi_RadTC_dj1 =   2.0132000000;
        // const double GrtAlmandi_RadTC_jmin = -23.025850930; 
        // const double GrtAlmandi_RadTC_jmax =  1.4486839760;

        // Coefficients for garnet (majorite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, GCubed, 20(5), 2378-2394]
        // https://doi.org/10.1029/2019GC008187
        // mineral composition [Mg3(MgSi)(SiO4)3]
        // const double GrtMajorit_RadTC_c0  =  -11.261000000;
        // const double GrtMajorit_RadTC_dj1 =   2.0132000000;
        // const double GrtMajorit_RadTC_jmin = -23.025850930; 
        // const double GrtMajorit_RadTC_jmax =  1.4486839760;

        // Coefficients for quartz
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        // const double QuartzPure_RadTC_c0  =   0;
        // const double QuartzPure_RadTC_dj1 =   0;
        // const double QuartzPure_RadTC_jmin = -23.025850930; 
        // const double QuartzPure_RadTC_jmax = -23.050000000;

        // Coefficients for coesite
        // assumed 0 for now - no data available
        // mineral composition [SiO2]
        // const double CoesitSiO2_RadTC_c0  =   0;
        // const double CoesitSiO2_RadTC_dj1 =   0;
        // const double CoesitSiO2_RadTC_jmin = -23.025850930; 
        // const double CoesitSiO2_RadTC_jmax = -23.050000000;

        // Define room temperature [K] 
        const double T_room = 298.15; 

        const unsigned int n_points = in.n_evaluation_points();

        for (unsigned int i = 0; i < n_points; ++i) 
        {
          // Convert pressure unit from [Pa] to [GPa]
          double P_GPa = in.pressure[i]/1e9;

          // Compute natural logarithm of pressure and temperature 
          double P_log = std::log(P_GPa);
          double T_log = std::log(in.temperature[i]);

          // Take the temperature field of the model [K]
          double T_mod = in.temperature[i];

          // Take the mineral fraction of the model
          double min_frac = in.composition[0][i];

          // Compute the lattice thermal conductivity in real (+,-) and simplex (0->1) 
          // space considering the boundaries ymin and ymax
          // Dry Olivine
          double OlivineDry_LatTC_zSimpl = OlivineDry_LatTC_a0 + OlivineDry_LatTC_b1*P_log;
          double OlivineDry_LatTC_ySimpl = std::exp(OlivineDry_LatTC_zSimpl);
          double OlivineDry_LatTC_yPrime = OlivineDry_LatTC_ySimpl/(1+OlivineDry_LatTC_ySimpl);
          double OlivineDry_LatTC_yRealS = OlivineDry_LatTC_ymin+(OlivineDry_LatTC_ymax-OlivineDry_LatTC_ymin)*OlivineDry_LatTC_yPrime;
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_LatTC_zSimpl = OpxEnstati_LatTC_a0 + OpxEnstati_LatTC_b1*P_log;
          // double OpxEnstati_LatTC_ySimpl = std::exp(OpxEnstati_LatTC_zSimpl);
          // double OpxEnstati_LatTC_yPrime = OpxEnstati_LatTC_ySimpl/(1+OpxEnstati_LatTC_ySimpl);
          // double OpxEnstati_LatTC_yRealS = OpxEnstati_LatTC_ymin+(OpxEnstati_LatTC_ymax-OpxEnstati_LatTC_ymin)*OpxEnstati_LatTC_yPrime;
          // Garnet (Pyrope)
          // double GrtPyropes_LatTC_zSimpl = GrtPyropes_LatTC_a0 + GrtPyropes_LatTC_b1*P_log;
          // double GrtPyropes_LatTC_ySimpl = std::exp(GrtPyropes_LatTC_zSimpl);
          // double GrtPyropes_LatTC_yPrime = GrtPyropes_LatTC_ySimpl/(1+GrtPyropes_LatTC_ySimpl);
          // double GrtPyropes_LatTC_yRealS = GrtPyropes_LatTC_ymin+(GrtPyropes_LatTC_ymax-GrtPyropes_LatTC_ymin)*GrtPyropes_LatTC_yPrime;

          // Compute the P-dependent lattice thermal conductivity of minerals 
          // Dry Olivine
          double OlivineDry_PDep_LatTCon = std::exp(OlivineDry_LatTC_yRealS);
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_PDep_LatTCon = std::exp(OpxEnstati_LatTC_yRealS);
          // Garnet (Pyrope)
          // double GrtPyropes_PDep_LatTCon = std::exp(GrtPyropes_LatTC_yRealS);

          // Compute T-dependent lattice thermal conductivity of minerals 
          // Dry Olivine
          double OlivineDry_TDep_LatTCon = OlivineDry_PDep_LatTCon*std::pow((T_room/T_mod),OlivineDry_TDep_n_Exp);
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_TDep_LatTCon = OpxEnstati_PDep_LatTCon*std::pow((T_room/T_mod),OpxEnstati_TDep_n_Exp);
          // Garnet (Pyrope)
          // double GrtPyropes_TDep_LatTCon = GrtPyropes_PDep_LatTCon*std::pow((T_room/T_mod),GrtPyropes_TDep_n_Exp);

          // Compute P,T-dependent lattice thermal conductivity of minerals 
          // Dry Olivine
          double OlivineDry_PTDep_LatTCo = OlivineDry_TDep_LatTCon;
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_PTDep_LatTCo = OpxEnstati_TDep_LatTCon;
          // Garnet (Pyrope)
          // const double GrtPyropes_PTDep_LatTCo = GrtPyropes_TDep_LatTCon;

          // Compute the radiative thermal conductivity in real (+,-) and simplex (0->1) 
          // space considering the boundaries ymin and ymax
          // Dry Olivine
          double OlivineDry_RadTC_zSimpl = OlivineDry_RadTC_c0 + OlivineDry_RadTC_d1*T_log;
          double OlivineDry_RadTC_ySimpl = std::exp(OlivineDry_RadTC_zSimpl);
          double OlivineDry_RadTC_yPrime = OlivineDry_RadTC_ySimpl/(1+OlivineDry_RadTC_ySimpl);
          double OlivineDry_RadTC_yRealS = OlivineDry_RadTC_jmin+(OlivineDry_RadTC_jmax-OlivineDry_RadTC_jmin)*OlivineDry_RadTC_yPrime;
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_RadTC_zSimpl = OpxEnstati_RadTC_c0 + OpxEnstati_RadTC_d1*T_log;
          // double OpxEnstati_RadTC_ySimpl = std::exp(OpxEnstati_RadTC_zSimpl);
          // double OpxEnstati_RadTC_yPrime = OpxEnstati_RadTC_ySimpl/(1+OpxEnstati_RadTC_ySimpl);
          // double OpxEnstati_RadTC_yRealS = OpxEnstati_RadTC_jmin+(OpxEnstati_RadTC_jmax-OpxEnstati_RadTC_jmin)*OpxEnstati_RadTC_yPrime;
          // Garnet (Pyrope)
          // double GrtPyropes_RadTC_zSimpl = GrtPyropes_RadTC_c0 + GrtPyropes_RadTC_d1*T_log;
          // double GrtPyropes_RadTC_ySimpl = std::exp(GrtPyropes_RadTC_zSimpl);
          // double GrtPyropes_RadTC_yPrime = GrtPyropes_RadTC_ySimpl/(1+GrtPyropes_RadTC_ySimpl);
          // double GrtPyropes_RadTC_yRealS = GrtPyropes_RadTC_jmin+(GrtPyropes_RadTC_jmax-GrtPyropes_RadTC_jmin)*GrtPyropes_RadTC_yPrime;

          // Compute T-dependent lattice thermal conductivity of minerals
          // Dry Olivine
          double OlivineDry_TDep_RadTCon = std::exp(OlivineDry_RadTC_yRealS);
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_TDep_RadTCon = std::exp(OpxEnstati_RadTC_yRealS);
          // Garnet (Pyrope)
          // double GrtPyropes_TDep_RadTCon = std::exp(GrtPyropes_RadTC_yRealS);
    
          // Compute P,T-dependent total thermal conductivity of minerals 
          // Dry Olivine
          double OlivineDry_PTDep_TotTCo = OlivineDry_PTDep_LatTCo+OlivineDry_TDep_RadTCon;
          // Orthopyroxene (Enstatite)
          // double OpxEnstati_PTDep_TotTCo = OpxEnstati_PTDep_LatTCo+OpxEnstati_TDep_RadTCon;
          // Garnet (Pyrope)
          // double GrtPyropes_PTDep_TotTCo = GrtPyropes_PTDep_LatTCo+GrtPyropes_TDep_RadTCon;

          // Compute P,T-dependent thermal conductivities of aggregate rocks 
          // double AggRock_PTDep_LatTCo = std::pow(OlivineDry_PTDep_LatTCo,in.composition[i])*std::pow(OpxEnstati_PTDep_LatTCo,in.composition[i])*std::pow(GrtPyropes_PTDep_LatTCo,in.composition[i]);
          // double AggRock_PTDep_RadTCo = std::pow(OlivineDry_TDep_RadTCon,in.composition[i])*std::pow(OpxEnstati_TDep_RadTCon,in.composition[i])*std::pow(GrtPyropes_TDep_RadTCon,in.composition[i]);
          double AggRock_PTDep_LatTCo = std::pow(OlivineDry_PTDep_LatTCo, min_frac);
          double AggRock_PTDep_RadTCo = std::pow(OlivineDry_TDep_RadTCon, min_frac);
          double AggRock_PTDep_TotTCo = AggRock_PTDep_LatTCo+AggRock_PTDep_RadTCo;

          // out.lat_thermal_conductivi[i] = OlivineDry_PTDep_LatTCo;
          // out.rad_thermal_conductivi[i] = OlivineDry_TDep_RadTCon;
          out.thermal_conductivities[i] = AggRock_PTDep_TotTCo;
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
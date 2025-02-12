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
        // mineral composition [(Mg1.8Fe0.2)SiO4]
        // const double WadsleyDry_LatTC_a0 =  -0.665600000;
        // const double WadsleyDry_LatTC_b1 =   0.308200000;
        // const double WadsleyDry_LatTC_ymin = 1.75735804249439;
        // const double WadsleyDry_LatTC_ymax = 2.37090451537473; 
        // const double WadsleyDry_TDep_n_Exp =  0.5; 

        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of 
        // [Schloessin & Dvorak, 1972, GJI, 27(5), 499-516]
        // mineral composition [Mg2Si2O6]
        // const double OpxEnstati_LatTC_a0 =   -3.004700000;
        // const double OpxEnstati_LatTC_b1 =    2.600000000;
        // const double OpxEnstati_LatTC_ymin =  1.760865151; 
        // const double OpxEnstati_LatTC_ymax =  2.096937429;
        // const double OpxEnstati_TDep_n_Exp = -0.5;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Hung et al. 2024, American Mineralogist, 109(3), 482-487]
        // mineral composition [Mg3Al2Si3O12]
        // const double GrtPyropes_LatTC_a0 =   -4.363700000;
        // const double GrtPyropes_LatTC_b1 =    2.036800000;
        // const double GrtPyropes_LatTC_ymin =  1.481604541; 
        // const double GrtPyropes_LatTC_ymax =  2.443131606;
        // const double GrtPyropes_TDep_n_Exp = -0.4314;

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
        // mineral composition [Mg1.8 Fe0.2 SiO4]
        // const double WadsleyDry_RadTC_a0 =   -21.717000000;
        // const double WadsleyDry_RadTC_b1 =    3.4271000000;
        // const double WadsleyDry_RadTC_ymin = -23.025850930;
        // const double WadsleyDry_RadTC_ymax =  1.0773006810; 

        // Coefficients for orthopyroxene (enstatite)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, Journal of Petrology, vol. 60(4), p. 755-790]
        // mineral composition [Mg2Si2O6]
        // const double OpxEnstati_RadTC_c0 =   -13.532000000;
        // const double OpxEnstati_RadTC_d1 =    2.4004000000;
        // const double OpxEnstati_RadTC_jmin = -23.025850930; 
        // const double OpxEnstati_RadTC_jmax =  1.4456685920;

        // Coefficients for garnet (pyrope)
        // retreived from fitting dataset of
        // [Grose & Afonso, 2019, Journal of Petrology, vol. 60(4), p. 755-790]
        // mineral composition [Mg3Al2(SiO4)3]
        // const double GrtPyropes_RadTC_c0 =   -11.782000000;
        // const double GrtPyropes_RadTC_d1 =    2.0718000000;
        // const double GrtPyropes_RadTC_jmin = -23.025850930; 
        // const double GrtPyropes_RadTC_jmax =  1.4479836950;

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
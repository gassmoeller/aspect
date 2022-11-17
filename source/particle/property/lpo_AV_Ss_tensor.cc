/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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

//#include <cstdlib>
#include <aspect/particle/property/lpo_AV_Ss_tensor.h>
#include <aspect/particle/property/lpo.h>
#include <aspect/particle/world.h>

#include <aspect/utilities.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {


      template <int dim>
      LpoSsTensor<dim>::LpoSsTensor ()
      {
        permutation_operator_3d[0][1][2]  = 1;
        permutation_operator_3d[1][2][0]  = 1;
        permutation_operator_3d[2][0][1]  = 1;
        permutation_operator_3d[0][2][1]  = -1;
        permutation_operator_3d[1][0][2]  = -1;
        permutation_operator_3d[2][1][0]  = -1;



        // tensors of indices
        indices_tensor[0][0] = 0;
        indices_tensor[0][1] = 5;
        indices_tensor[0][2] = 4;
        indices_tensor[1][0] = 5;
        indices_tensor[1][1] = 1;
        indices_tensor[1][2] = 3;
        indices_tensor[2][0] = 4;
        indices_tensor[2][1] = 3;
        indices_tensor[2][2] = 2;

        // vectors of indices
        indices_vector_1.resize(6);
        indices_vector_1[0] = 0;
        indices_vector_1[1] = 1;
        indices_vector_1[2] = 2;
        indices_vector_1[3] = 1;
        indices_vector_1[4] = 2;
        indices_vector_1[5] = 0;

        indices_vector_2.resize(6);
        indices_vector_2[0] = 0;
        indices_vector_2[1] = 1;
        indices_vector_2[2] = 2;
        indices_vector_2[3] = 2;
        indices_vector_2[4] = 0;
        indices_vector_2[5] = 1;
      }

      template <int dim>
      void
      LpoSsTensor<dim>::initialize ()
      {
        // todo: check wheter this works correctly. Since the get_random_number function takes a reference
        // to the random_number_generator function, changing the function should mean that I have to update the
        // get_random_number function as well. But I will need to test this.
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        this->random_number_generator.seed(random_number_seed+my_rank);

        const auto &manager = this->get_particle_world().get_property_manager();
        AssertThrow(manager.plugin_name_exists("lpo"),
                    ExcMessage("No lpo property plugin found."));
        Assert(manager.plugin_name_exists("lpo Ss tensor"),
               ExcMessage("No LPO aniso stress property plugin found."));

        AssertThrow(manager.check_plugin_order("lpo","lpo Ss tensor"),
                    ExcMessage("To use the lpo Ss tensor plugin, the lpo plugin need to be defined before this plugin."));

        lpo_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("lpo"));

        if (dim == 2)
          {
            Assert(false,ExcMessage("This PROPERTY is not implemented for 2D."));

          }


      }


      template <int dim>
      SymmetricTensor<2,dim>
      LpoSsTensor<dim>::compute_S_tensor (const SymmetricTensor<2,dim> &strain_rate,
                                          const double grain_size,
                                          std::vector<double> volume_fraction_mineral,
                                          std::vector<std::vector<double>> volume_fractions_grains,
                                          const std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains,
                                          const std::vector<unsigned int> &deformation_type,
                                          const double &temperature) const
      {
        Assert(false,ExcMessage("This PROPERTY is not implemented for 2D."));
      }
      /*template <>
      SymmetricTensor<2,2>
      LpoSsTensor<2>::compute_S_tensor (const SymmetricTensor<2,2> &strain_rate,
                                          const double &grain_size,
                                          const std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains,
                                          const std::vector<unsigned int> &deformation_type,
                                          const double &temperature) const;

      {
          Assert(false,ExcMessage("This PROPERTY is not implemented for 2D."));
      }*/
      template <>
      SymmetricTensor<2,3>
      LpoSsTensor<3>::compute_S_tensor (const SymmetricTensor<2,3> &strain_rate,
                                        const double grain_size,
                                        std::vector<double> volume_fraction_mineral,
                                        std::vector<std::vector<double>> volume_fractions_grains,
                                        const std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_grains,
                                        const std::vector<unsigned int> &deformation_type,
                                        const double &temperature) const

      {
        const size_t n_minerals_local = a_cosine_matrices_grains.size();
        //std::cout<<"n_minerals_local: "<< n_minerals_local<<  std::endl;
        const size_t n_grains_local = a_cosine_matrices_grains[0].size();
        //std::cout<<"n_grains_local: "<< n_grains_local<<  std::endl;
        SymmetricTensor<2,3, double> S_sum;
        const int dim=3;
        double nFo = 4.1;
        double A0 = 1.1e5*std::exp(-530000/8.314/temperature);
        //std::cout<<"T: "<<temperature<<std::endl;
        //std::cout<<"A0: "<<A0<<std::endl;
        FullMatrix<double> Schm(6,3); //Schmid tensor, 6x3 matrix
        FullMatrix<double> pinvschm(3,6); //pseudoinverse of Schmid tensor, 3x6 matrix


        Schm[3][2] = 1;
        Schm[4][1] = 1;
        Schm[5][0] = 1;
        pinvschm[0][5] = 1;
        pinvschm[1][4] = 1;
        pinvschm[2][3] = 1;

        for (size_t mineral_i = 0; mineral_i < n_minerals_local; mineral_i++)
          {
            //std::cout<<"Def style: "<<deformation_type[mineral_i]<<std::endl;
            Tensor<1,3> A_ss; //A_ss is the invers of the minimum resolved stress on the slip systems on the nth power
            //std::cout<<"Volume fractions minerals: "<< volume_fraction_mineral[mineral_i]<<  std::endl;
            if (deformation_type[mineral_i] == (unsigned int)DeformationTypeSelector::Enstatite)
              {
                A_ss[0] = 1.;
                A_ss[1] = 1.;
                A_ss[2] = 1.;
                //std::cout<<"Setting Enstatite A_ss: "<<A_ss<<std::endl;
              }
            else
              {
                A_ss[0] = 139.2525;
                A_ss[1] = 214.4907;
                A_ss[2] = 0.3520;
                //std::cout<<"Setting Olivine A_ss: "<<A_ss<<std::endl;
              }

            //std::cout<<"A_ss: "<<A_ss<<std::endl;
            for (size_t i = 0; i < n_grains_local; i++) //NOTE:  need to use volume fractions for the final stress per particle for enstatite and olivine
              {
                //std::cout<<"strain rate: "<<strain_rate<<  "A_ss: "<<A_ss<<  std::endl;
                //std::cout<<"A_ss: "<<A_ss<<  std::endl;
                //std::cout<<"Volume fraction grains: "<< volume_fractions_grains[mineral_i][i]<<  std::endl;
                Tensor<2,3> R = a_cosine_matrices_grains[mineral_i][i];
                //std::cout<<"Rotation matrix: "<<R<<  std::endl;
                SymmetricTensor<2,3> Rate_grain=symmetrize(R*strain_rate*transpose(R));
                //std::cout<<"Rate_grain "<<Rate_grain<<  std::endl;
                std::array<std::pair<double, Tensor<1, 3>>, 3> Rate_gr_eig = eigenvectors(Rate_grain,SymmetricTensorEigenvectorMethod::jacobi);
                //std::cout<<"Rate_grain eigen values: "<<Rate_gr_eig[0].first<<"; "<<Rate_gr_eig[1].first<<"; "<<Rate_gr_eig[2].first<<std::endl;
                double inv2=std::pow(Rate_gr_eig[0].first-Rate_gr_eig[1].first,2)
                            +std::pow(Rate_gr_eig[1].first-Rate_gr_eig[2].first,2)
                            +std::pow(Rate_gr_eig[2].first-Rate_gr_eig[0].first,2);
                //std::cout<<"inv2: "<<inv2<<  std::endl;

                FullMatrix<double> Rate_grain_voigt(6,1);
                Rate_grain_voigt[0][0]=Rate_grain[0][0];
                Rate_grain_voigt[1][0]=Rate_grain[1][1];
                Rate_grain_voigt[2][0]=Rate_grain[2][2];
                Rate_grain_voigt[3][0]=2*Rate_grain[1][2];
                Rate_grain_voigt[4][0]=2*Rate_grain[0][2];
                Rate_grain_voigt[5][0]=2*Rate_grain[0][1];

                //std::cout<<"Rate_grain_voigt: "<<Rate_grain_voigt[0][0]<<"; "<<Rate_grain_voigt[1][0]<<"; "<<Rate_grain_voigt[2][0]<<"; "<<Rate_grain_voigt[3][0]<<"; "<<Rate_grain_voigt[4][0]<<"; "<<Rate_grain_voigt[5][0]<<std::endl;

                FullMatrix<double> r_ss(3,1); //Optimazition to find shear strain rate on slip system
                FullMatrix<double> r_gc_v(6,1); //strain rate tensor for grain in Voigt notation and crystal reference frame
                pinvschm.mmult(r_ss, Rate_grain_voigt);
                //std::cout<<"r_ss"<<r_ss[0][0]<<"; "<<r_ss[1][0]<<"; "<<r_ss[2][0]<<std::endl;
                Schm.mmult(r_gc_v,r_ss);

                SymmetricTensor<2,3> r_gc;
                r_gc[0][0]=r_gc_v[0][0];
                r_gc[1][1]=r_gc_v[1][0];
                r_gc[2][2]=r_gc_v[2][0];
                r_gc[1][2]=0.5*r_gc_v[3][0];
                r_gc[0][2]=0.5*r_gc_v[4][0];
                r_gc[0][1]=0.5*r_gc_v[5][0];
                //std::cout<<"r_gc:  "<<r_gc<<  std::endl;
                std::array<std::pair<double, Tensor<1, dim>>, dim> r_gc_eig = eigenvectors(r_gc, SymmetricTensorEigenvectorMethod::jacobi);

                double inv2best =std::pow(r_gc_eig[0].first-r_gc_eig[1].first,2)
                                 +std::pow(r_gc_eig[1].first-r_gc_eig[2].first,2)
                                 +std::pow(r_gc_eig[0].first-r_gc_eig[2].first,2);

                for (unsigned int i=0; i<dim; ++i)
                  {
                    r_ss[i][0]=r_ss[i][0]*std::pow(inv2/inv2best,0);
                  }
                FullMatrix<double> tau_ss(3,1);
                //std::cout<<"std::copysignf(1.0,r_ss[0][0]): "<<std::copysignf(1.0,r_ss[0][0])<<  std::endl;
                //std::cout<<"1.0/A_ss[0]: "<<1.0/A_ss[0]<<  std::endl;
                //std::cout<<"1.0/A0: "<<1.0/A0<<std::endl;
                AssertThrow(isfinite(1./A0),
                            ExcMessage("1/A0 is infinite"))
                tau_ss[0][0]= std::copysignf(1.0,r_ss[0][0])*std::pow(1.0/A_ss[0]*1.0/A0*std::pow(grain_size,0.73)*std::fabs(r_ss[0][0]/2),1.0/nFo);
                tau_ss[1][0]= std::copysignf(1.0,r_ss[1][0])*std::pow(1.0/A_ss[1]*1.0/A0*std::pow(grain_size,0.73)*std::fabs(r_ss[1][0]/2),1.0/nFo);
                tau_ss[2][0]= std::copysignf(1.0,r_ss[2][0])*std::pow(1.0/A_ss[2]*1.0/A0*std::pow(grain_size,0.73)*std::fabs(r_ss[2][0]/2),1.0/nFo);


                FullMatrix<double>  S_gc_v(6,1);
                Schm.mmult(S_gc_v,tau_ss); //Voigt notation of the resolved stress on the grain
                SymmetricTensor<2,3> S_gc;
                S_gc[0][0] = S_gc_v[0][0];
                S_gc[1][1] = S_gc_v[1][0];
                S_gc[2][2] = S_gc_v[2][0];
                S_gc[1][2] = S_gc_v[3][0];
                S_gc[0][2] = S_gc_v[4][0];
                S_gc[0][1] = S_gc_v[5][0];

                SymmetricTensor<2,3> S_g= symmetrize(transpose(R)*S_gc*R); //Here instead of making a multidimensional array what I sum at the end, I create S_g and add it to S_sum
                //SymmetricTensor<2,3> S_sum;
                //std::cout<<"Stress on grain: "<<S_g<<  std::endl;
                SymmetricTensor<2,3> S_g_contrib = S_g*volume_fraction_mineral[mineral_i]*volume_fractions_grains[mineral_i][i]; //Each particle has n grains that have multiple mineral fractions (with some volume fraction)
                //std::cout<<"gtains contribution of total stress MPa: "<<S_g_contrib<<  std::endl;
                S_sum += S_g_contrib;
                //std::cout<<"S_sum: "<<S_sum<<  std::endl;

              }


          }

        S_sum *= 1e6;
        //std::cout<<"S_sum final Pa: "<<S_sum<<  std::endl;

        return S_sum;
      }



      template <int dim>
      void
      LpoSsTensor<dim>::initialize_one_particle_property(const Point<dim> &,
                                                         std::vector<double> &data) const
      {
        std::vector<unsigned int> deformation_type;
        std::vector<double> volume_fraction_mineral;
        std::vector<std::vector<double>> volume_fractions_grains;
        std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains;

        Particle::Property::LPO<dim>::load_particle_data(lpo_data_position,
                                                         data,
                                                         deformation_type,
                                                         volume_fraction_mineral,
                                                         volume_fractions_grains,
                                                         a_cosine_matrices_grains);



        Tensor<2,6> Ss_tensor; //The Ss tensor is a compilation of the stresses needed for the calculation of the viscosity tensor

        // There is a bug up to dealii 9.3.0, so we have to work around it.
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)

          {
            data.push_back(Ss_tensor[Tensor<2,6>::unrolled_to_component_indices(i)]);
          }


      }

      template <int dim>
      void
      LpoSsTensor<dim>::update_one_particle_property(const unsigned int data_position,
                                                     const Point<dim> &,
                                                     const Vector<double> &solution,
                                                     const std::vector<Tensor<1,dim> > &gradients,
                                                     const ArrayView<double> &data) const
      {
        std::vector<unsigned int> deformation_type;
        std::vector<double> volume_fraction_mineral;
        std::vector<std::vector<double>> volume_fractions_grains;
        std::vector<std::vector<Tensor<2,3> > > a_cosine_matrices_grains;

        Particle::Property::LPO<dim>::load_particle_data(lpo_data_position,
                                                         data,
                                                         deformation_type,
                                                         volume_fraction_mineral,
                                                         volume_fractions_grains,
                                                         a_cosine_matrices_grains);

        //std::cout<<"a_cos_matrices -size1: "<< a_cosine_matrices_grains.size()<< "a_cos_matrices -size2: "<< a_cosine_matrices_grains[0].size()<< std::endl;
        //std::cout<<"deformation style -size: "<< deformation_type.size()<<  std::endl;

        const double grain_size=1000.0; //micron --> should be an input?
        Tensor<2,6> Ss_tensor; //Initial value 0, because at initial timestep we don't have strain rate
        double temperature = solution[this->introspection().component_indices.temperature];

        if  (this->get_timestep_number() > 0 && temperature > 1000)
          {
            Tensor<2,dim> velocity_gradient;
            for (unsigned int d=0; d<dim; ++d)
              {
                velocity_gradient[d] = gradients[d];
              }

            const SymmetricTensor<2,dim> strain_rate = symmetrize (velocity_gradient);
            //std::cout<<"strain rate: "<< strain_rate<< std::endl;
            double E_eq;
            SymmetricTensor<2,dim> e1, e2, e3, e4, e5, E;
            E=strain_rate;
            E_eq=(1.0/6.0*(std::pow(double (E[0][0] - E[1][1]),2) + std::pow(double (E[1][1] - E[2][2]),2)+std::pow(double (E[2][2] - E[0][0]),2)))+(std::pow(E[0][1],2)+std::pow(E[1][2],2)+std::pow(E[2][0],2));//J2
            E_eq= std::sqrt((4./3.)*E_eq);// Second invariant of strain-rate

            AssertThrow(isfinite(1/E.norm()),
                        ExcMessage("Strain rate should be finite"));

            //We define 5 independent strainrates, of which E is the linear combination
            e1[0][0]=E_eq;
            e1[1][1]=E_eq;
            e1[2][2]=-2*E_eq;
            e2[0][0]=E_eq;
            e2[1][1]=-2*E_eq;
            e2[2][2]=E_eq;
            e3[0][1]=E_eq;
            e3[1][0]=E_eq;
            e4[0][2]=E_eq;
            e4[2][0]=E_eq;
            e5[1][2]=E_eq;
            e5[2][1]=E_eq;

            //We calculate the stress response for each strain rate with the micromechanical model
            // AssertThrow(in.temperature[q] != 0,
            //     ExcMessage("Temperature is 0"));
            //std::cout<<"e1: "<< e1<<std::endl;
            //std::cout<<"grain size: "<< grain_size<<std::endl;
            //std::cout<<"T: "<< temperature<<std::endl;
            SymmetricTensor<2,dim> stress1, stress2, stress3, stress4, stress5, Stress;
            stress1=compute_S_tensor(e1, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            stress2=compute_S_tensor(e2, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            stress3=compute_S_tensor(e3, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            stress4=compute_S_tensor(e4, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            stress5=compute_S_tensor(e5, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            Stress =compute_S_tensor(E, grain_size, volume_fraction_mineral, volume_fractions_grains, a_cosine_matrices_grains, deformation_type, temperature);
            //std::cout << "Strain rate particle " << E << std::endl;
            //std::cout << "Stress tensor particle " << Stress << std::endl;


            for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
              {
                Ss_tensor[0][i] = Stress[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
                Ss_tensor[1][i] = stress1[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
                Ss_tensor[2][i] = stress2[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
                Ss_tensor[3][i] = stress3[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
                Ss_tensor[4][i] = stress4[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
                Ss_tensor[5][i] = stress5[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];

              }
          }
        //std::cout << "Ss tensor " << Ss_tensor << std::endl;
        Particle::Property::LpoSsTensor<dim>::store_particle_data(data_position,
                                                                  data,
                                                                  Ss_tensor);


      }


      template <int dim>
      void
      LpoSsTensor<dim>::load_particle_data(unsigned int lpo_data_position,
                                           const ArrayView<double> &data,
                                           Tensor<2,6> &Ss_tensor)
      {

        // There is a bug up to dealii 9.3.0, so we have to work around it.
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)
          {
            Ss_tensor[Tensor<2,6>::unrolled_to_component_indices(i)] = data[lpo_data_position + i];
          }

        //for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
        //elastic_tensor[SymmetricTensor<2,6>::unrolled_to_component_indices(i)] = data[lpo_data_position + i];
      }


      template <int dim>
      void
      LpoSsTensor<dim>::store_particle_data(unsigned int lpo_data_position,
                                            const ArrayView<double> &data,
                                            Tensor<2,6> &Ss_tensor)
      {
        // There is a bug up to dealii 9.3.0, so we have to work around it.
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)
          {
            data[lpo_data_position + i] = Ss_tensor[Tensor<2,6>::unrolled_to_component_indices(i)];
          }

      }



      template <int dim>
      UpdateTimeFlags
      LpoSsTensor<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      UpdateFlags
      LpoSsTensor<dim>::get_needed_update_flags () const
      {
        return update_default;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      LpoSsTensor<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information;

        property_information.push_back(std::make_pair("lpo_Ss_tensor",Tensor<2,6>::n_independent_components));

        return property_information;
      }

      template <int dim>
      void
      LpoSsTensor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LpoSsTensor");
            {
              prm.declare_entry ("Random number seed", "1",
                                 Patterns::Integer (0),
                                 "The seed used to generate random numbers. This will make sure that "
                                 "results are reproducable as long as the problem is run with the "
                                 "same amount of MPI processes. It is implemented as final seed = "
                                 "user seed + MPI Rank. ");


              prm.declare_entry ("Volume fraction olivine", "0.5",
                                 Patterns::Double(0),
                                 "The volume fraction of the olivine phase (0 is no olivine, 1 is fully olivine). "
                                 "The rest of the volume fraction is set to be entstatite. "
                                 "Todo: if full olivine make not enstite grains and vice-versa.");

              prm.declare_entry ("Number of samples", "0",
                                 Patterns::Double(0),
                                 "This determines how many samples are taken when using the random "
                                 "draw volume averaging. Setting it to zero means that the number of "
                                 "samples is set to be equal to the number of grains.");
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }


      template <int dim>
      void
      LpoSsTensor<dim>::parse_parameters (ParameterHandler &prm)
      {

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LpoSsTensor");
            {

              random_number_seed = prm.get_integer ("Random number seed"); // 2
              n_grains = LPO<dim>::get_number_of_grains();
              n_minerals = LPO<dim>::get_number_of_minerals();
              n_samples = prm.get_integer("Number of samples"); // 0
              if (n_samples == 0)
                n_samples = n_grains;
            }
            prm.leave_subsection ();
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();


      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(LpoSsTensor,
                                        "lpo Ss tensor",
                                        "A plugin in which the particle property tensor is "
                                        "defined as he collection of stresses resulted from "
                                        "the micromechanical model for olivine aggregate deformation "
                                        "with the current strain rate and 5 independent strain rates "
                                        "with the same amplitude. These stresses can be used to construct "
                                        "the full rank4 viscosity tensor.")
    }
  }
}

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
#include <aspect/particle/property/lpo_visco_tensor.h>
#include <aspect/particle/property/lpo.h>
#include <aspect/particle/world.h>

#include <aspect/utilities.h>
#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace internal
  {
    
  }
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      LpoViscoTensor<dim>::LpoViscoTensor ()
      {
       
      }

      template <int dim>
      void
      LpoViscoTensor<dim>::initialize ()
      {
        // todo: check wheter this works correctly. Since the get_random_number function takes a reference
        // to the random_number_generator function, changing the function should mean that I have to update the
        // get_random_number function as well. But I will need to test this.
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
        this->random_number_generator.seed(random_number_seed+my_rank);

        const auto &manager = this->get_particle_world().get_property_manager();
        AssertThrow(manager.plugin_name_exists("lpo"),
                    ExcMessage("No lpo property plugin found."));
        Assert(manager.plugin_name_exists("lpo visco tensor"),
               ExcMessage("No viscosity anisotropy property plugin found."));

        AssertThrow(manager.check_plugin_order("lpo","lpo visco tensor"),
                    ExcMessage("To use the lpo visco tensor plugin, the lpo plugin need to be defined before this plugin."));

        lpo_data_position = manager.get_data_info().get_position_by_plugin_index(manager.get_plugin_index_by_name("lpo"));


      }

      template <int dim>
      SymmetricTensor<2,dim>
      LpoViscoTensor<dim>::Stress_strain_aggregate((const SymmetricTensor<2,dim,double> rate,const std::vector<std::vector<Tensor<2,dim, double> > > R_matrix, const double temperature, const std::vector<std::vector<double> > grain_size, const std::vector<unsigned int> deformation_type));

      {
        if(dim == 2)
        {
          Assert(false,ExcMessage("This PROPERTY is not implemented for 2D."));
        }
        else
        {


    
      
          /*Micromechanical model for olivine deformation by Hansen et al., (2016,
          JGR) using a pseudo-Taylor method, assuming that each grain experiences
          the same strain rate. It results in the best fitting stress in MPa that we convert at the end to unit in Pa.*/
          const int n_grains=R_matrix[0].size();
          const int m_minerals=R_matrix.size();
          const int dim=3;
          double nFo = 4.1;
          double A0 = 1.1e5*std::exp(-530000/8.314/temperature);
          FullMatrix<double> Schm(6,3); //Schmid tensor, 6x3 matrix
          FullMatrix<double> pinvschm(3,6); //pseudoinverse of Schmid tensor, 3x6 matrix
          Tensor<1,3> A_ss; //A_ss is the invers of the minimum resolved stress on the slip systems on the nth power
          
          Schm[3][2] = 1;
          Schm[4][1] = 1;
          Schm[5][0] = 1;
          pinvschm[0][5] = 1;
          pinvschm[1][4] = 1;
          pinvschm[2][3] = 1;
          SymmetricTensor<2,3> S_sum;
          for (int i_mineral =0; i_mineral < m_minerals; i_mineral++)
          {
            if (deformation_type[i_mineral] == (unsigned int)DeformationTypeSelector::Enstatite) //Decide about what to do with Ensatite
            {
              A_ss[0] = 1;
              A_ss[1] = 1;
              A_ss[2] = 1;
            }
            else
            {
              A_ss[0] = 139.2525;
              A_ss[1] = 214.4907;
              A_ss[2] = 0.3520;
            } 
            for (int i = 0; i < n_grains; i++) {
              Tensor<2,3> R = R_matrix[i_mineral][i];
              SymmetricTensor<2,3> Rate_grain=symmetrize(R*rate*transpose(R));
              std::array<std::pair<double, Tensor<1, 3>>, 3> Rate_gr_eig = eigenvectors(Rate_grain,SymmetricTensorEigenvectorMethod::jacobi);
		          double inv2=std::pow(Rate_gr_eig[0].first-Rate_gr_eig[1].first,2)
                          +std::pow(Rate_gr_eig[1].first-Rate_gr_eig[2].first,2)
                          +std::pow(Rate_gr_eig[2].first-Rate_gr_eig[0].first,2);
              FullMatrix<double> Rate_grain_voigt;
              Rate_grain_voigt[0][0]=Rate_grain[0][0];
              Rate_grain_voigt[1][0]=Rate_grain[1][1];
              Rate_grain_voigt[2][0]=Rate_grain[2][2];
              Rate_grain_voigt[3][0]=2*Rate_grain[1][2];
              Rate_grain_voigt[4][0]=2*Rate_grain[0][2];
              Rate_grain_voigt[5][0]=2*Rate_grain[0][1];
              
              FullMatrix<double> r_ss(3,1); //Optimazition to find shear strain rate on slip system
              FullMatrix<double> r_gc_v(6,1); //strain rate tensor for grain in Voigt notation and crystal reference frame
              pinvschm.mmult(r_ss, Rate_grain_voigt);
              Schm.mmult(r_gc_v,r_ss);
        
              SymmetricTensor<2,3> r_gc; 
              r_gc[0][0]=r_gc_v[0][0];
              r_gc[1][1]=r_gc_v[1][0];
              r_gc[2][2]=r_gc_v[2][0];
              r_gc[1][2]=0.5*r_gc_v[3][0];
              r_gc[0][2]=0.5*r_gc_v[4][0];
              r_gc[0][1]=0.5*r_gc_v[5][0];
        
		          std::array<std::pair<double, Tensor<1, dim>>, dim> r_gc_eig = eigenvectors(r_gc, SymmetricTensorEigenvectorMethod::jacobi);
		      
              double inv2best =std::pow(r_gc_eig[0].first-r_gc_eig[1].first,2)
                              +std::pow(r_gc_eig[1].first-r_gc_eig[2].first,2)
                              +std::pow(r_gc_eig[0].first-r_gc_eig[2].first,2);

		          for (unsigned int i=0; i<dim; ++i)
              {
                r_ss[i][0]=r_ss[i][0]*std::pow(inv2/inv2best,0);
              }
		      
              FullMatrix<double> tau_ss(3,1);
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
		          
              S_sum += S_g; 
            }//end loop  for number of grains
            S_sum=S_sum/n_grains; //Stress for mineralphase
          } //end loop  for number of mineral phases
          S_sum=S_sum/m_minerals; //stress for particle
          S_sum *= 1e6; //convert from MPa to Pa 
        }
      }  
        

        
     
      
      template <int dim>
      void
      LpoViscoTensor<dim>::initialize_one_particle_property(const Point<dim> &,
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

        
        Tensor<2,6> SsTensor; //The Ss tensor is a compilation of the stresses needed for the calculation of the viscosity tensor

        
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)

          {
            data.push_back(SsTensor[Tensor<2,6>::unrolled_to_component_indices(i)]);
          }


      }

      template <int dim>
      void
      LpoViscoTensor<dim>::update_one_particle_property(const unsigned int data_position,
                                                          const Point<dim> &position,
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


        const size_t n_minerals_local = volume_fractions_grains.size();
        const size_t n_grains_local = volume_fractions_grains[0].size();
        const double grainsize_init=1000.0; //micron --> should be an input?
        Tensor<2,dim> velocity_gradient;
        for (unsigned int d=0; d<dim; ++d)
        {
          velocity_gradient[d] = gradients[d]; 
        }
        double temperature = solution[this->introspection().component_indices.temperature];
        const SymmetricTensor<2,dim> strain_rate = symmetrize (velocity_gradient);
        std::vector<std::vector<double> > grain_size;
        for (size_t mineral_i = 0; mineral_i < n_minerals_local; mineral_i++)
        {
          for (unsigned int grain_i = 0; grain_i < n_grains_local; ++grain_i)
          {
            grain_size[mineral_i][grain_i]=n_minerals_local*n_grains_local*volume_fractions_grains[mineral_i][grain_i]*grainsize_init;
          }
        }
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
        // AssertThrow(grain_size == 1000,
        //     ExcMessage("Something is wrong with the grain_size"));
        SymmetricTensor<2,dim> stress1, stress2, stress3, stress4, stress5, Stress;
        stress1=LpoViscoTensor::Stress_strain_aggregate(e1, a_cosine_matrices_grains, temperature, grain_size);
        stress2=LpoViscoTensor::Stress_strain_aggregate(e2, a_cosine_matrices_grains, temperature, grain_size);
        stress3=LpoViscoTensor::Stress_strain_aggregate(e3, a_cosine_matrices_grains, temperature, grain_size);
        stress4=LpoViscoTensor::Stress_strain_aggregate(e4, a_cosine_matrices_grains, temperature, grain_size);
        stress5=LpoViscoTensor::Stress_strain_aggregate(e5, a_cosine_matrices_grains, temperature, grain_size);
        Stress =LpoViscoTensor::Stress_strain_aggregate(E, a_cosine_matrices_grains, temperature, grain_size);
        Tensor<2,6> SsTensor;
        for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
        {
          SsTensor[0][i] = Stress[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
          SsTensor[1][i] = stress1[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
          SsTensor[2][i] = stress2[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
          SsTensor[3][i] = stress3[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
          SsTensor[4][i] = stress4[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];
          SsTensor[5][i] = stress5[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)];

        }  

        Particle::Property::LpoViscoTensor<dim>::store_particle_data(data_position,
                                                                       data,
                                                                       SsTensor);


      }


      template <int dim>
      void
      LpoViscoTensor<dim>::load_particle_data(unsigned int lpo_data_position,
                                                const ArrayView<double> &data,
                                                SymmetricTensor<2,6> &SsTensor)
      {

        // There is a bug up to dealii 9.3.0, so we have to work around it.
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)

          SsTensor[Tensor<2,6>::unrolled_to_component_indices(i)] = data[lpo_data_position + i];

        //for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
        //elastic_tensor[SymmetricTensor<2,6>::unrolled_to_component_indices(i)] = data[lpo_data_position + i];
      }


      template <int dim>
      void
      LpoViscoTensor<dim>::store_particle_data(unsigned int lpo_data_position,
                                                 const ArrayView<double> &data,
                                                 SymmetricTensor<2,6> &SsTensor)
      {
        // There is a bug up to dealii 9.3.0, so we have to work around it.
        for (unsigned int i = 0; i < Tensor<2,6>::n_independent_components ; ++i)

          data[lpo_data_position + i] = SsTensor[Tensor<2,6>::unrolled_to_component_indices(i)];

        //for (unsigned int i = 0; i < SymmetricTensor<2,6>::n_independent_components ; ++i)
        //  data[lpo_data_position + i] = elastic_tensor[SymmetricTensor<2,6>::unrolled_to_component_indices(i)];
      }


      template <int dim>
      UpdateTimeFlags
      LpoViscoTensor<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      UpdateFlags
      LpoViscoTensor<dim>::get_needed_update_flags () const
      {
        return update_default;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      LpoViscoTensor<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information;

        property_information.push_back(std::make_pair("LPO_AV_stresses",Tensor<2,6>::n_independent_components));

        return property_information;
      }

      template <int dim>
      void
      LpoViscoTensor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LpoViscoTensor");
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
      LpoViscoTensor<dim>::parse_parameters (ParameterHandler &prm)
      {

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("LpoViscoTensor");
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(LpoViscoTensor,
                                        "LPO_AV_stresses",
                                        "A plugin in which the particle property tensor is "
                                        "defined as a compilation of anisotropic stresses "
                                        "that are required for calculating the full viscosity tensor")
    }
  }
}

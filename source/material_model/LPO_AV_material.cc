
/*
 Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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

#include <aspect/material_model/LPO_AV_material.h>
#include <aspect/introspection.h>
#include <aspect/material_model/interface.h>
#include <aspect/plugins.h>
#include <aspect/simulator/assemblers/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

//#include <world_builder/utilities.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/geometry_info.h>
#include <aspect/simulator_access.h>

#include <aspect/material_model/simple.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator_signals.h>
#include <aspect/postprocess/particles.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
//#include <aspect/postprocess/particle_lpo.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>



namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for anisotropic viscosities to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class AnisotropicViscosity : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        AnisotropicViscosity(const unsigned int n_points);

        static double J2_second_invariant(const SymmetricTensor<2,dim> t, const double min_strain_rate);

        static FullMatrix<double> transform_Symmetric3x3_matrix_to_6D_vector(const SymmetricTensor<2,3> &tensor);

        static SymmetricTensor<2,3> transform_6D_vector_to_Symmetric3x3_matrix(const FullMatrix<double> &matrix); 

        static Tensor <2,3> euler_angles_to_rotation_matrix(double phi1_d, double theta_d, double phi2_d);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Stress-strain "director" tensors at the given positions. This
         * variable is used to implement anisotropic viscosity.
         *
         * @note The strain rate term in equation (1) of the manual will be
         * multiplied by this tensor *and* the viscosity scalar ($\eta$ /i.e. effective viscosity), as
         * described in the manual section titled "Constitutive laws". This
         * variable is assigned the rank-four identity tensor by default.
         * This leaves the isotropic constitutive law unchanged if the material
         * model does not explicitly assign a value.
         */
        std::vector<SymmetricTensor<4,dim> > stress_strain_directors;
    };

    namespace
    {

      template <int dim>
      std::vector<std::string> make_AnisotropicViscosity_additional_outputs_names()
      {
        std::vector<std::string> names;

        for (unsigned int i = 0; i < Tensor<4,dim>::n_independent_components ; ++i)
          {
            TableIndices<4> indices(Tensor<4,dim>::unrolled_to_component_indices(i));
            names.push_back("anisotropic_viscosity"+std::to_string(indices[0])+std::to_string(indices[1])+std::to_string(indices[2])+std::to_string(indices[3]));
          }
        return names;
      }
    }



    template <int dim>
    AnisotropicViscosity<dim>::AnisotropicViscosity (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_AnisotropicViscosity_additional_outputs_names<dim>()),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ())
    {}



    template <int dim>
    std::vector<double>
    AnisotropicViscosity<dim>::get_nth_output(const unsigned int idx) const
    {
      std::vector<double> output(stress_strain_directors.size());
      for (unsigned int i = 0; i < stress_strain_directors.size() ; ++i)
        {
          output[i]= stress_strain_directors[i][Tensor<4,dim>::unrolled_to_component_indices(idx)];
        }
      return output;
    }
  }
}

namespace aspect
{
  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the Stokes preconditioner.
     */
    template <int dim>
    class StokesPreconditionerAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AnisotropicViscosities.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Stokes equation for the current cell.
     */
    template <int dim>
    class StokesIncompressibleTermsAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AdditionalMaterialOutputsStokesRHS if we need to do so.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };



    template <int dim>
    void
    StokesPreconditionerAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      // First loop over all dofs and find those that are in the Stokes system
      // save the component (pressure and dim velocities) each belongs to.
      for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
            {
              scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
              ++i_stokes;
            }
          ++i;
        }

      // Loop over all quadrature points and assemble their contributions to
      // the preconditioner matrix
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] =
                    scratch.finite_element_values[introspection.extractors
                                                  .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;
          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += (2.0 * eta * (scratch.grads_phi_u[i]
                                                         * stress_strain_director
                                                         * scratch.grads_phi_u[j])
                                            + one_over_eta * pressure_scaling
                                            * pressure_scaling
                                            * (scratch.phi_p[i]
                                               * scratch.phi_p[j]))
                                           * JxW;
        }
    }



    template <int dim>
    void
    StokesPreconditionerAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  if (scratch.rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }
          // Viscosity scalar
          const double eta = (scratch.rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());

          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                   * JxW;

              if (force != nullptr)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (scratch.rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                // assemble \nabla p as -(p, div v):
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * JxW;
                  }
            }
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }
      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
             == n_points, ExcInternalError());
    }
  }

  namespace HeatingModel
  {
    template <int dim>
    class ShearHeatingAnisotropicViscosity : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const;
    };



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.strain_rate.size(),
             ExcMessage ("The shear heating plugin needs the strain rate!"));

      // Some material models provide dislocation viscosities and boundary area work fractions
      // as additional material outputs. If they are attached, use them.
      const MaterialModel::DislocationViscosityOutputs<dim> *disl_viscosities_out =
        material_model_outputs.template get_additional_output<MaterialModel::DislocationViscosityOutputs<dim> >();

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // If there is an anisotropic viscosity, use it to compute the correct stress
          const SymmetricTensor<2,dim> &directed_strain_rate = ((anisotropic_viscosity != nullptr)
                                                                ?
                                                                anisotropic_viscosity->stress_strain_directors[q]
                                                                * material_model_inputs.strain_rate[q]
                                                                :
                                                                material_model_inputs.strain_rate[q]);

          const SymmetricTensor<2,dim> stress =
            2 * material_model_outputs.viscosities[q] *
            (this->get_material_model().is_compressible()
             ?
             directed_strain_rate - 1./3. * trace(directed_strain_rate) * unit_symmetric_tensor<dim>()
             :
             directed_strain_rate);

          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            (this->get_material_model().is_compressible()
             ?
             material_model_inputs.strain_rate[q]
             - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
             :
             material_model_inputs.strain_rate[q]);

          heating_model_outputs.heating_source_terms[q] = stress * deviatoric_strain_rate;

          // If dislocation viscosities and boundary area work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to increase surface energy)
          if (disl_viscosities_out != 0)
            {
              heating_model_outputs.heating_source_terms[q] *= 1 - disl_viscosities_out->boundary_area_change_work_fractions[q] *
                                                               material_model_outputs.viscosities[q] / disl_viscosities_out->dislocation_viscosities[q];
            }

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      const unsigned int n_points = material_model_outputs.viscosities.size();

      if (material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          material_model_outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }
  }
  
  /*
  namespace MaterialModel
  {
    /* The LPO_AV material model calculates an anisotropic viscosity tensor from the orientation of the olivine grains. 
    We first calculate a V tensor that is stress independent, and from that we create the stress_strain_directors, 
    that will be used in the Assemblers (2*eta_eff*s-s-d*strain rate). To get the V tensor we will use the micromechanical model 
    by Hansen et al., 2016, similarly as it was described in Kiraly et al., 2020 *
    template <int dim>
    class LPO_AV : public MaterialModel::Simple<dim>
    {
      public:
        virtual void initialize();
        virtual void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters (ParameterHandler &prm);
        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;
        virtual double reference_density () const;
        virtual void create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const;
      private:
        double eta; //reference viscosity
        double min_strain_rate;
        double grain_size; 
        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;
    };
  }
  */
}

namespace aspect
{

//Next session is a more evolved implementation of anisotropic viscosity in the material model based on Hansen et al 2016 and Kiraly et al 2020
  namespace MaterialModel
  {
    template <int dim>
    void
    LPO_AV<dim>::set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesPreconditioner<dim>>(*(assemblers.stokes_preconditioner[i])))
            assemblers.stokes_preconditioner[i] = std_cxx14::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim> > ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesIncompressibleTerms<dim>>(*(assemblers.stokes_system[i])))
            assemblers.stokes_system[i] = std_cxx14::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim> > ();
        }
    }



    template <int dim>
    void
    LPO_AV<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&LPO_AV<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
      AssertThrow((dim==3),
                  ExcMessage("Olivine has 3 independent slip systems, allowing for deformation in 3 independent directions, hence these models only work in 3D"));

    }

    namespace internal
    {
      template <int dim>
      SymmetricTensor<2,dim>
      Stress_strain_aggregate(const SymmetricTensor<2,dim,double> rate,const Tensor<1,dim> euler, const double temperature, const double grain_size)
      {
        Assert(false,ExcMessage("This material model is not implemented for 2D."));
      }

      //grainsize - either input, or eventually need to get from particles grain volume fraction; 
      //ALSO LATER - get directly the rotaion matrix from particles instead of euler angles. ?Updating euler angles?
      template <> //Specialization of the above function for 3D
      SymmetricTensor<2,3>
      Stress_strain_aggregate(SymmetricTensor<2,3,double> rate,Tensor<1,3,double> euler,double temperature,double grain_size)
      {
        /*Micromechanical model for olivine deformation by Hansen et al., (2016,
        JGR) using a pseudo-Taylor method, assuming that each grain experiences
        the same strain rate. It results in the best fitting stress in MPa that we convert at the end to unit in Pa.*/
        const int dim=3;
        double nFo = 4.1;
        double A0 = 1.1e5*std::exp(-530000/8.314/temperature);
        FullMatrix<double> Schm(6,3); //Schmid tensor, 6x3 matrix
        FullMatrix<double> pinvschm(3,6); //pseudoinverse of Schmid tensor, 3x6 matrix
        Tensor<1,3> A_ss; //A_ss is the invers of the minimum resolved stress on the slip systems on the nth power
        A_ss[0] = 139.2525;
        A_ss[1] = 214.4907;
        A_ss[2] = 0.3520;
        Schm[3][2] = 1;
        Schm[4][1] = 1;
        Schm[5][0] = 1;
        pinvschm[0][5] = 1;
        pinvschm[1][4] = 1;
        pinvschm[2][3] = 1;
        //for (int i = 0; i < ngrains; i++) {
          Tensor<2,3> R = AnisotropicViscosity<dim>::euler_angles_to_rotation_matrix(euler[0],euler[1],euler[2]);
          SymmetricTensor<2,3> Rate_grain=symmetrize(R*rate*transpose(R));
          std::array<std::pair<double, Tensor<1, 3>>, 3> Rate_gr_eig = eigenvectors(Rate_grain,SymmetricTensorEigenvectorMethod::jacobi);
		      double inv2=std::pow(Rate_gr_eig[0].first-Rate_gr_eig[1].first,2)
          +std::pow(Rate_gr_eig[1].first-Rate_gr_eig[2].first,2)
          +std::pow(Rate_gr_eig[0].first-Rate_gr_eig[2].first,2);
          FullMatrix<double> Rate_grain_voigt = AnisotropicViscosity<dim>::transform_Symmetric3x3_matrix_to_6D_vector(Rate_grain);

          // Either like this
          //Vector<double> r_ss(3);
          //Vector<double> rate_grain_vector(6) = ...;
          //pinvschm.Tvmult(r_ss,Rate_grain_voigt);

          // Or like this
          FullMatrix<double> r_ss(3,1);
          FullMatrix<double> r_gc_v(6,1);
          pinvschm.mmult(r_ss, Rate_grain_voigt);
          Schm.mmult(r_gc_v,r_ss);
            
          //FullMatrix<double> r_ss(3,1) = pinvschm * Rate_grain_voigt; //!!Problem: Rate_grain_voigt is not a column vector
          //Tensor<1,6, double> r_gc_v = Schm * r_ss;
          auto r_gc = AnisotropicViscosity<dim>::transform_6D_vector_to_Symmetric3x3_matrix(r_gc_v);
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
		      //error: no matching function for call to 'copysign', and/or 'fabs'
          FullMatrix<double>  S_gc_v(6,1);
          Schm.mmult(S_gc_v,tau_ss); //Voigt notation of the resolved stress on the grain
		      SymmetricTensor<2,3> S_gc = AnisotropicViscosity<dim>::transform_6D_vector_to_Symmetric3x3_matrix(S_gc_v);
		      SymmetricTensor<2,3> S_g= symmetrize(transpose(R)*S_gc*R); //Here instead of making a multidimensional array what I sum at the end, I create S_g and add it to S_sum
		      SymmetricTensor<2,3> S_sum;
          S_sum += S_g; //error: no viable overloaded '+='
        //}//end loop  for number of grains
        S_sum *= 1e6; //convert from MPa to Pa

        return S_sum;
      }
    }
    

    template <int dim> 
    double
    AnisotropicViscosity<dim>::J2_second_invariant(const SymmetricTensor<2,dim> t, const double min_strain_rate)
    {
      const double J2_strict = 1/6*(std::pow(t[0][0] - t[1][1],2) + std::pow(t[1][1] - t[2][2],2)+std::pow(t[2][2] - t[0][0],2))+std::pow(t[0][1],2)+std::pow(t[1][2],2)+std::pow(t[2][0],2); 
      const double J2 = std::max(J2_strict, std::pow(min_strain_rate,2)); //prevents having too small values (also used in compute_second_invariant for strain rate)
      return J2;
    } 

    template<int dim> 
    FullMatrix<double>
    AnisotropicViscosity<dim>::transform_Symmetric3x3_matrix_to_6D_vector(const SymmetricTensor<2,3> &input) //error: cannot define or redeclare 'transform_Symmetric3x3_matrix_to_6D_vector' here because namespace 'MaterialModel' does not enclose namespace 'AnisotropicViscosity'
    {
      FullMatrix<double> result(6,1); 
      
      result[0][0]= input[0][0];                  // 0 - 1 
      result[1][0]= input[1][1];                  // 1  // 2
      result[2][0]= input[2][2];                  // 2  // 3
      result[3][0]= std::sqrt(2.0)*input[1][2];   // 3  // 4
      result[4][0]= std::sqrt(2.0)*input[0][2];   // 4  // 5
      result[5][0]= std::sqrt(2.0)*input[0][1];   // 5  // 6
      return result;
    }

    template<int dim> //Where do I have to declare this function?
    SymmetricTensor<2,3>
    AnisotropicViscosity<dim>::transform_6D_vector_to_Symmetric3x3_matrix(const FullMatrix<double> &input) 
    {
      SymmetricTensor<2,3> result;

      const double sqrt_2_inv = 1/std::sqrt(2.0);

      result[0][0] = input[0][0];
      result[1][1] = input[1][0];
      result[2][2] = input[2][0];
      result[1][2] = sqrt_2_inv * input[3][0];
      result[0][2] = sqrt_2_inv * input[4][0];
      result[0][1] = sqrt_2_inv * input[5][0];

      return result;

    }
    template<int dim>  
    Tensor<2,3>
    AnisotropicViscosity<dim>::euler_angles_to_rotation_matrix(double phi1_d, double theta_d, double phi2_d)
    {
      constexpr double const_pi = 3.141592653589793238462643383279502884;
      const double degree_to_rad = const_pi/180.0;
      const double phi1 = phi1_d * degree_to_rad;
      const double theta = theta_d * degree_to_rad;
      const double phi2 = phi2_d * degree_to_rad;
      Tensor<2, 3> rot_matrix;


      rot_matrix[0][0] = cos(phi2)*cos(phi1) - cos(theta)*sin(phi1)*sin(phi2);
      rot_matrix[0][1] = -cos(phi2)*sin(phi1) - cos(theta)*cos(phi1)*sin(phi2);
      rot_matrix[0][2] = -sin(phi2)*sin(theta);

      rot_matrix[1][0] = sin(phi2)*cos(phi1) + cos(theta)*sin(phi1)*cos(phi2);
      rot_matrix[1][1] = -sin(phi2)*sin(phi1) + cos(theta)*cos(phi1)*cos(phi2);
      rot_matrix[1][2] = cos(phi2)*sin(theta);

      rot_matrix[2][0] = -sin(theta)*sin(phi1);
      rot_matrix[2][1] = -sin(theta)*cos(phi1);
      rot_matrix[2][2] = cos(theta);
      return rot_matrix;
    }

    template <int dim> 
    void
    LPO_AV<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity;
        anisotropic_viscosity = out.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      AssertThrow((this->introspection().compositional_name_exists("euler1")),
                  ExcMessage("LPO_AV material model only works if there is a compositional field called euler1 (first euler angle)."));
      AssertThrow(this->introspection().compositional_name_exists("euler2"),
                  ExcMessage("LPO_AV material model only works if there is a compositional field called euler2 (second euler angle)."));
      AssertThrow(this->introspection().compositional_name_exists("euler3"),
                  ExcMessage("LPO_AV material model only works if there is a compositional field called euler3 (third euler angle)."));

      std::vector<unsigned int> c_idx_euler;
      c_idx_euler.push_back (this->introspection().compositional_index_for_name("euler1"));
      c_idx_euler.push_back (this->introspection().compositional_index_for_name("euler2"));
      c_idx_euler.push_back (this->introspection().compositional_index_for_name("euler3"));

      // Get the grad_u tensor, at the center of this cell, if possible. DO I NEED THIS? I NEED THE CURRENT STRAIN RATE

      std::vector<Tensor<2,dim> > velocity_gradients (in.n_evaluation_points());
      if (in.current_cell.state() == IteratorState::valid)
        {
          std::vector<Point<dim> > quadrature_positions(in.n_evaluation_points());
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

          // FEValues requires a quadrature and we provide the default quadrature
          // as we only need to evaluate the gradients of the solution.
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   Quadrature<dim>(quadrature_positions),
                                   update_gradients);
          fe_values.reinit (in.current_cell);
          fe_values[this->introspection().extractors.velocities]
          .get_function_gradients(this->get_solution(), velocity_gradients);
        }

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          out.densities[q] = 0;//Change this to 0 for the simple shear box test
          out.viscosities[q] = eta; //Later it is going to be overwritten by the effective viscosity
          out.thermal_expansion_coefficients[q] = 1e-10;
          out.specific_heat[q] = 1;
          out.thermal_conductivities[q] = 1;
          out.compressibilities[q] = 0.0;
          out.entropy_derivative_pressure[q] = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          for (unsigned int c=0; c<in.composition[q].size(); ++c)
            out.reaction_terms[q][c] = 0.0;

          Tensor<1,dim> euler;
          for (unsigned int i=0; i<dim; ++i){
            euler[i] = in.composition[q][c_idx_euler[i]];}

          // The computation of the viscosity tensor is only
          // necessary after the simulator has been initialized
          if  ((this->simulator_is_past_initialization()) && (this->get_timestep_number() > 0) && (in.temperature[q]>1000))
            {
              double E_eq;
              SymmetricTensor<2,dim> e1, e2, e3, e4, e5, E;
              
                E_eq= std::sqrt(4/3*AnisotropicViscosity<dim>::J2_second_invariant(in.strain_rate[q], min_strain_rate)); //Second invariant of strain-rate
                E=in.strain_rate[q];
                
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
              SymmetricTensor<2,dim> stress1, stress2, stress3, stress4, stress5, Stress, s1, s2, s3, s4, s5, S;
              stress1=internal::Stress_strain_aggregate(e1, euler, in.temperature[q],grain_size);
              stress2=internal::Stress_strain_aggregate(e2, euler, in.temperature[q],grain_size);
              stress3=internal::Stress_strain_aggregate(e3, euler, in.temperature[q],grain_size);
              stress4=internal::Stress_strain_aggregate(e4, euler, in.temperature[q],grain_size);
              stress5=internal::Stress_strain_aggregate(e5, euler, in.temperature[q],grain_size);
              Stress=internal::Stress_strain_aggregate(E, euler, in.temperature[q],grain_size);
              const double Stress_eq= std::sqrt(3*AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate));
              AssertThrow(Stress_eq != 0,
                  ExcMessage("Equivalent stress should not be 0"));
              AssertThrow(isfinite(Stress_eq),
                  ExcMessage("Stress should be finite")); 
              //To calculate a "stress independent viscosity" (i.e. inverse of fluidty) 
              //we have to convert the stress to a "non-Newtonian modified stress" which is the stress *second invariant on the power of (n-1)/2
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  s1[i][j]= stress1[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress1, min_strain_rate),1.25);//assuming n=3.5, so (n-1)/2=1.25
                  s2[i][j]= stress2[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress2, min_strain_rate),1.25);
                  s3[i][j]= stress3[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress3, min_strain_rate),1.25);
                  s4[i][j]= stress4[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress4, min_strain_rate),1.25);
                  s5[i][j]= stress5[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(stress5, min_strain_rate),1.25);
                  S[i][j]= Stress[i][j]*std::pow(AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate),1.25);
                }
              }
              
              //Build the stress independent V tensor 
              SymmetricTensor<4,dim> V, ViscoTensor_r4;
              for (int i = 0; i < dim; i++)
               {
                for (int j = 0; j < dim; j++)
                {
                  V[i][j][2][2] = (S[i][j] - ((2/3*s1[i][j]/E_eq)+(1/3*s2[i][j]/E_eq)) *E[0][0] 
                                         - ((1/3*s1[i][j]/E_eq)-(1/3*s2[i][j]/E_eq)) *E[1][1]
                                         - s3[i][j]*E[0][1]/E_eq - s4[i][j]*E[0][2]/E_eq - s5[i][j]*E[1][2]/E_eq)/
                                         (E[0][0]/3 + 2*E[1][1]/3 + E[2][2]);
                  V[i][j][0][0] = (2/3* s1[i][j] + 1/3*s2[i][j])/E_eq + 1/3*V[i][j][2][2];
                  V[i][j][1][1] = (1/3* s1[i][j] - 1/3*s2[i][j])/E_eq + 2/3*V[i][j][2][2];
                  V[i][j][0][1] = 0.5*s3[i][j]/E_eq;
                  V[i][j][0][2] = 0.5*s4[i][j]/E_eq;
                  V[i][j][1][2] = 0.5*s5[i][j]/E_eq; 
                  for (int k = 0; k<dim; k++)
                  {
                    for (int l = 0; l<dim; l++)
                    {
                      ViscoTensor_r4[i][j][k][l]= V[i][j][k][l]/std::pow(AnisotropicViscosity<dim>::J2_second_invariant(Stress, min_strain_rate),1.25);
                    }
                  }
                }
              }
              
              // Overwrite the scalar viscosity with an effective viscosity
              out.viscosities[q] = std::abs(Stress_eq/E_eq);
              AssertThrow(out.viscosities[q] != 0,
                  ExcMessage("EViscosity should not be 0")); 
              AssertThrow(isfinite(out.viscosities[q]),
                  ExcMessage("EViscosity should not be finite")); 
              if (anisotropic_viscosity != nullptr)
              {
                anisotropic_viscosity->stress_strain_directors[q] = ViscoTensor_r4/(Stress_eq/E_eq);
                SymmetricTensor<2,6> ViscoTensor;
                ViscoTensor[0][0]=anisotropic_viscosity->stress_strain_directors[q][0][0][0][0] * out.viscosities[q];
                ViscoTensor[0][1]=anisotropic_viscosity->stress_strain_directors[q][0][0][1][1] * out.viscosities[q];
                ViscoTensor[0][2]=anisotropic_viscosity->stress_strain_directors[q][0][0][2][2] * out.viscosities[q];
                ViscoTensor[0][3]=anisotropic_viscosity->stress_strain_directors[q][0][0][1][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[0][4]=anisotropic_viscosity->stress_strain_directors[q][0][0][0][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[0][5]=anisotropic_viscosity->stress_strain_directors[q][0][0][0][1] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[1][1]=anisotropic_viscosity->stress_strain_directors[q][1][1][1][1] * out.viscosities[q];
                ViscoTensor[1][2]=anisotropic_viscosity->stress_strain_directors[q][1][1][2][2] * out.viscosities[q];
                ViscoTensor[1][3]=anisotropic_viscosity->stress_strain_directors[q][1][1][1][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[1][4]=anisotropic_viscosity->stress_strain_directors[q][1][1][0][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[1][5]=anisotropic_viscosity->stress_strain_directors[q][1][1][0][1] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[2][2]=anisotropic_viscosity->stress_strain_directors[q][2][2][2][2] * out.viscosities[q];
                ViscoTensor[2][3]=anisotropic_viscosity->stress_strain_directors[q][2][2][1][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[2][4]=anisotropic_viscosity->stress_strain_directors[q][2][2][0][2] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[2][5]=anisotropic_viscosity->stress_strain_directors[q][2][2][0][1] * std::sqrt(2) * out.viscosities[q];
                ViscoTensor[3][3]=anisotropic_viscosity->stress_strain_directors[q][1][2][1][2] * 2 * out.viscosities[q];
                ViscoTensor[3][4]=anisotropic_viscosity->stress_strain_directors[q][1][2][0][2] * 2 * out.viscosities[q];
                ViscoTensor[3][5]=anisotropic_viscosity->stress_strain_directors[q][1][2][0][1] * 2 * out.viscosities[q];
                ViscoTensor[4][4]=anisotropic_viscosity->stress_strain_directors[q][0][2][0][2] * 2 * out.viscosities[q];
                ViscoTensor[4][5]=anisotropic_viscosity->stress_strain_directors[q][0][2][0][1] * 2 * out.viscosities[q];
                ViscoTensor[5][5]=anisotropic_viscosity->stress_strain_directors[q][0][1][0][1] * 2 * out.viscosities[q];


                }
            }

      }      
    }




    template <int dim>
    bool
    LPO_AV<dim>::is_compressible () const
    {
      return false;
    }



    template <int dim>
    double
    LPO_AV<dim>::reference_density () const
    {
      return 1.0;
    }



    template <int dim>
    double
    LPO_AV<dim>::reference_viscosity () const
    {
      return 1e20;
    }



    template <int dim>
    void
    LPO_AV<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AnisotropicViscosity");
        {
          eta = prm.get_double("Reference viscosity");
          min_strain_rate = prm.get_double("Minimum strain rate");
          grain_size = prm.get_double("Grain size");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LPO_AV<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AnisotropicViscosity");
        {
          prm.declare_entry ("Reference viscosity", "1e20",
                             Patterns::Double(),
                             "Magnitude of reference viscosity.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}"); 
          prm.declare_entry ("Grain size", "1000",
                             Patterns::Double(),
                             "Olivine anisotropic viscosity is dependent of grain size. Value is given in microns");                  
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LPO_AV<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<AnisotropicViscosity<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class StokesPreconditioner<dim>; \
  template class StokesIncompressibleTerms<dim>; \
  template class StokesBoundaryTraction<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingAnisotropicViscosity,
                                  "anisotropic shear heating",
                                  "Implementation of a standard model for shear heating. "
                                  "Adds the term: "
                                  "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                  "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                  "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                  "right-hand side of the temperature equation.")
  }



  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(LPO_AV,
                                   "AnisotropicViscosity material",
                                   "Olivine LPO related viscous anisotropy")
  }
}

/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/heating_model/shear_heating.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ShearHeating<dim>::
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
      const ShearHeatingOutputs<dim> *shear_heating_out =
        material_model_outputs.template get_additional_output<ShearHeatingOutputs<dim>>();

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            (this->get_material_model().is_compressible()
             ?
             material_model_inputs.strain_rate[q]
             - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
             :
             material_model_inputs.strain_rate[q]);

          SymmetricTensor<2,dim> stress =
            2 * material_model_outputs.viscosities[q] *
            deviatoric_strain_rate;

          heating_model_outputs.heating_source_terms[q] = stress * deviatoric_strain_rate;

          // If elasticity is used, the stress should include the elastic stresses
          // and only the visco-plastic (non-elastic) strain rate should contribute.
          if (this->get_parameters().enable_elasticity == true)
            {
              // Visco-elastic stresses are stored on the fields
              SymmetricTensor<2, dim> stress_0, stress_old;
              stress_0[0][0] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
              stress_0[1][1] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
              stress_0[0][1] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

              stress_old[0][0] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx_old")];
              stress_old[1][1] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy_old")];
              stress_old[0][1] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy_old")];

              if (dim == 3)
                {
                  stress_0[2][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                  stress_0[0][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                  stress_0[1][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];

                  stress_old[2][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz_old")];
                  stress_old[0][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz_old")];
                  stress_old[1][2] = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz_old")];
                }

              const MaterialModel::ElasticAdditionalOutputs<dim> *elastic_out = material_model_outputs.template get_additional_output<MaterialModel::ElasticAdditionalOutputs<dim>>();
              AssertThrow(elastic_out != nullptr, ExcMessage("ElasticAdditionalOutputs are requested, but have not been computed."));

              const double shear_modulus = elastic_out->elastic_shear_moduli[q];

              // Retrieve the elastic timestep and viscosity, only two material models
              // support elasticity.
              double elastic_timestep = 0.;
              double elastic_viscosity = 0.;
              if (Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim>>(this->get_material_model()))
                {
                  const MaterialModel::ViscoPlastic<dim> &vp = Plugins::get_plugin_as_type<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model());
                  elastic_viscosity = vp.get_elastic_viscosity(shear_modulus);
                  elastic_timestep = vp.get_elastic_timestep();
                }
              else if (Plugins::plugin_type_matches<MaterialModel::Viscoelastic<dim>>(this->get_material_model()))
                {
                  const MaterialModel::Viscoelastic<dim> &ve = Plugins::get_plugin_as_type<const MaterialModel::Viscoelastic<dim>>(this->get_material_model());
                  elastic_viscosity = ve.get_elastic_viscosity(shear_modulus);
                  elastic_timestep = ve.get_elastic_timestep();
                }
              else
                AssertThrow(false, ExcMessage("Shear heating cannot be used with elasticity for material models other than ViscoPlastic and Viscoelastic."));

              double dtc = this->get_timestep();
              if (dtc == 0 && this->get_timestep_number() == 0)
                dtc = std::min(std::min(this->get_parameters().maximum_time_step, this->get_parameters().maximum_first_time_step), elastic_timestep);
              const double timestep_ratio = dtc / elastic_timestep;
              // Scale the elastic viscosity with the timestep ratio, eta is already scaled.
              elastic_viscosity *= timestep_ratio;

              // Apply the stress update to get the total stress of timestep t.
              stress = 2. * material_model_outputs.viscosities[q] * deviatoric_strain_rate + material_model_outputs.viscosities[q] / elastic_viscosity * stress_0 +
                       (1. - timestep_ratio) * (1. - material_model_outputs.viscosities[q] / elastic_viscosity) * stress_old;

              SymmetricTensor<2, dim> visco_plastic_strain_rate = material_model_inputs.strain_rate[q] - ((stress - stress_0) / (2. * dtc * shear_modulus));

              if (this->get_material_model().is_compressible())
                visco_plastic_strain_rate = visco_plastic_strain_rate -
                                            1. / 3. * trace(visco_plastic_strain_rate) * unit_symmetric_tensor<dim>();

              heating_model_outputs.heating_source_terms[q] = stress * visco_plastic_strain_rate;
            }

          // If shear heating work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to be converted into other forms of energy)
          if (shear_heating_out != nullptr)
            {
              heating_model_outputs.heating_source_terms[q] *= shear_heating_out->shear_heating_work_fractions[q];
            }

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeating<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }



    template <int dim>
    ShearHeatingOutputs<dim>::ShearHeatingOutputs (const unsigned int n_points)
      :
      MaterialModel::NamedAdditionalMaterialOutputs<dim>({"shear_heating_work_fraction"}),
                  shear_heating_work_fractions(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    ShearHeatingOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void) idx;
      AssertIndexRange (idx, 1);

      return shear_heating_work_fractions;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeating,
                                  "shear heating",
                                  "Implementation of a standard model for shear heating. "
                                  "Adds the term: "
                                  "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                  "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                  "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                  "right-hand side of the temperature equation.")

#define INSTANTIATE(dim) \
  template class ShearHeatingOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

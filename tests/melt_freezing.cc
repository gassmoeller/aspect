#include <aspect/material_model/interface.h>
#include <aspect/melt.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

using namespace dealii;


namespace aspect
{
  template <int dim>
  class MeltingRate:
    public MaterialModel::MeltInterface<dim>, public ::aspect::SimulatorAccess<dim>
  {
      bool is_compressible () const override
      {
        return false;
      }

      double reference_viscosity () const override
      {
        return 5e20;
      }

      double reference_darcy_coefficient () const override
      {
        return 1e-8 * std::pow(0.01, 3.0) / 10.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0; i<in.position.size(); ++i)
          {
            out.viscosities[i] = 5e20;
            out.densities[i] = 3000.0;
            out.thermal_expansion_coefficients[i] = 2e-5;
            out.specific_heat[i] = 1250.0;
            out.thermal_conductivities[i] = 4.7;
            out.compressibilities[i] = 0.0;

            if (this->simulator_is_past_initialization() && this->get_timestep() > 0.0)
{
            const double Gamma = (in.position[i][1] > 0.5) ? 0.1 * (in.position[i][1]-0.5) : 0.0;
            // Porosity needs to be a mass reaction term
            out.reaction_terms[i][0] = -Gamma * out.densities[i] / this->get_timestep();
            out.reaction_terms[i][1] = Gamma * in.composition[i][2];
            out.reaction_terms[i][2] = 0.0;
}
else
{
            out.reaction_terms[i][0] = 0.0;
            out.reaction_terms[i][1] = 0.0;
            out.reaction_terms[i][2] = 0.0;
}
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != nullptr)
          {
            for (unsigned int i=0; i<in.position.size(); ++i)
              {
                melt_out->compaction_viscosities[i] = 5e20;
                melt_out->fluid_viscosities[i]= 10.0;
                melt_out->permeabilities[i]= 1e-8;
                melt_out->fluid_densities[i]= 3000.0;
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
              }
          }
      }

  };


}



// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(MeltingRate,
                                 "melting rate",
                                 "")
}

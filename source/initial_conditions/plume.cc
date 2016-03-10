/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/initial_conditions/plume.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <cmath>

namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    Plume<dim>::Plume ()
    {}

    template <int dim>
    void
    Plume<dim>::initialize ()
    {
      if (use_lithosphere_file)
        {
          AssertThrow ((dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model())) != 0,
                       ExcMessage ("This initial plume plugin can only be used with a lithosphere "
                                   "file when using a box geometry."));

          lithosphere_thickness.reset(new Utilities::AsciiDataLookup<dim-1> (1,1));

          const std::string filename = data_directory + data_file_name;

          this->get_pcout() << std::endl << "   Loading Ascii data initial file "
                            << filename << "." << std::endl << std::endl;

          if (Utilities::fexists(filename))
            lithosphere_thickness->load_file(filename);
          else
            AssertThrow(false,
                        ExcMessage (std::string("Ascii data file <")
                                    +
                                    filename
                                    +
                                    "> not found!"));
        }
    }

    template <int dim>
    double
    Plume<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // convert input ages to seconds
      const double age_top =    (this->convert_output_to_years() ? age_top_boundary_layer * year_in_seconds
                                 : age_top_boundary_layer);
      const double age_bottom = (this->convert_output_to_years() ? age_bottom_boundary_layer * year_in_seconds
                                 : age_bottom_boundary_layer);

      // First, get the temperature of the adiabatic profile at a representative
      // point at the top and bottom boundary of the model
      // if adiabatic heating is switched off, assume a constant profile
      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const Point<dim> bottom_point = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const double adiabatic_surface_temperature = this->get_adiabatic_conditions().temperature(surface_point);
      const double adiabatic_bottom_temperature = (this->include_adiabatic_heating())
                                                  ?
                                                  this->get_adiabatic_conditions().temperature(bottom_point)
                                                  :
                                                  adiabatic_surface_temperature;

      // then, get the temperature at the top and bottom boundary of the model
      // if no boundary temperature is prescribed simply use the adiabatic.
      // This implementation assumes that the top and bottom boundaries have
      // prescribed temperatures and minimal_temperature() returns the value
      // at the surface and maximal_temperature() the value at the bottom.
      const double T_surface = (&this->get_boundary_temperature() != 0)
                               ?
                               this->get_boundary_temperature().minimal_temperature(
                                 this->get_fixed_temperature_boundary_indicators())
                               :
                               adiabatic_surface_temperature;
      const double T_bottom = (&this->get_boundary_temperature() != 0)
                              ?
                              this->get_boundary_temperature().maximal_temperature(
                                this->get_fixed_temperature_boundary_indicators())
                              :
                              adiabatic_bottom_temperature;

      // get a representative profile of the compositional fields as an input
      // for the material model
      const double depth = this->get_geometry_model().depth(position);

      // look up material properties
      typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, this->n_compositional_fields());
      in.position[0]=position;
      in.temperature[0]=this->get_adiabatic_conditions().temperature(position);
      in.pressure[0]=this->get_adiabatic_conditions().pressure(position);
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        in.composition[0][c] = function->value(Point<1>(depth),c);
      in.strain_rate[0] = SymmetricTensor<2,dim>(); // adiabat has strain=0.
      this->get_material_model().evaluate(in, out);

      const double kappa = out.thermal_conductivities[0] / (out.densities[0] * out.specific_heat[0]);

      double surface_cooling_temperature;

      if (!use_lithosphere_file)
        {
          // analytical solution for the thermal boundary layer from half-space cooling model
          surface_cooling_temperature = age_top > 0.0 ?
                                        (T_surface - adiabatic_surface_temperature) *
                                        erfc(this->get_geometry_model().depth(position) /
                                             (2 * sqrt(kappa * age_top)))
                                        : 0.0;
        }
      else
        {
          Point<dim-1> data_position;
          for (unsigned int i = 0; i < dim - 1; i++)
            data_position[i] = position[i];

          const double lith_thickness = lithosphere_thickness->get_data(data_position,0);
          const double depth = this->get_geometry_model().depth(position);

          surface_cooling_temperature =  (T_surface - adiabatic_surface_temperature) *
                                         erfc(depth / lith_thickness);
        }

      const double bottom_heating_temperature = age_bottom > 0.0 ?
                                                (T_bottom - adiabatic_bottom_temperature + subadiabaticity)
                                                * erfc((this->get_geometry_model().maximal_depth()
                                                        - this->get_geometry_model().depth(position)) /
                                                       (2 * sqrt(kappa * age_bottom)))
                                                : 0.0;

      // set the initial temperature perturbation
      // first: get the center of the perturbation, then check the distance to the
      // evaluation point. the center is supposed to lie at the center of the bottom
      // surface.
      Point<dim> mid_point;
      if (perturbation_position == "center")
        {
          if (const GeometryModel::SphericalShell<dim> *
              shell_geometry_model = dynamic_cast <const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()))
            {
              const double inner_radius = shell_geometry_model->inner_radius();
              const double half_opening_angle = numbers::PI/180.0 * 0.5 * shell_geometry_model->opening_angle();
              if (dim==2)
                {
                  // choose the center of the perturbation at half angle along the inner radius
                  mid_point(0) = inner_radius * std::sin(half_opening_angle),
                  mid_point(1) = inner_radius * std::cos(half_opening_angle);
                }
              else if (dim==3)
                {
                  // if the opening angle is 90 degrees (an eighth of a full spherical
                  // shell, then choose the point on the inner surface along the first
                  // diagonal
                  if (shell_geometry_model->opening_angle() == 90)
                    {
                      mid_point(0) = inner_radius*std::sqrt(1./3),
                      mid_point(1) = inner_radius*std::sqrt(1./3),
                      mid_point(2) = inner_radius*std::sqrt(1./3);
                    }
                  else
                    {
                      // otherwise do the same as in 2d
                      mid_point(0) = inner_radius * std::sin(half_opening_angle) * std::cos(half_opening_angle),
                      mid_point(1) = inner_radius * std::sin(half_opening_angle) * std::sin(half_opening_angle),
                      mid_point(2) = inner_radius * std::cos(half_opening_angle);
                    }
                }
            }
          else if (const GeometryModel::Box<dim> *
                   box_geometry_model = dynamic_cast <const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
            {
              // for the box geometry, choose a point at the center of the bottom face.
              // (note that the loop only runs over the first dim-1 coordinates, leaving
              // the depth variable at zero)
              mid_point = box_geometry_model->get_origin();
              for (unsigned int i=0; i<dim-1; ++i)
                mid_point(i) += 0.5 * box_geometry_model->get_extents()[i];
            }
          else
            AssertThrow (false,
                         ExcMessage ("Not a valid geometry model for the initial conditions model"
                                     "adiabatic."));
        }
      else if (perturbation_position == "coordinates")
        {
          mid_point = midpoint_coordinates;
        }
      else
        AssertThrow(false, ExcNotImplemented());

      const double perturbation = (mid_point.distance(position) < radius) ? amplitude
                                  : 0.0;


      // add the subadiabaticity
      const double zero_depth = 0.174;
      const double nondimesional_depth = (this->get_geometry_model().depth(position) / this->get_geometry_model().maximal_depth() - zero_depth)
                                         / (1.0 - zero_depth);
      double subadiabatic_T = 0.0;
      if (nondimesional_depth > 0)
        subadiabatic_T = -subadiabaticity * nondimesional_depth * nondimesional_depth;

      // If adiabatic heating is disabled, apply all perturbations to
      // constant adiabatic surface temperature instead of adiabatic profile.
      const double temperature_profile = (this->include_adiabatic_heating())
                                         ?
                                         this->get_adiabatic_conditions().temperature(position)
                                         :
                                         adiabatic_surface_temperature;

      // return sum of the adiabatic profile, the boundary layer temperatures and the initial
      // temperature perturbation.
      return temperature_profile + surface_cooling_temperature
             + (perturbation > 0.0 ? std::max(bottom_heating_temperature + subadiabatic_T,perturbation)
                : bottom_heating_temperature + subadiabatic_T);
    }


    template <int dim>
    void
    Plume<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plume");
        {
          prm.declare_entry ("Age top boundary layer", "0e0",
                             Patterns::Double (0),
                             "The age of the upper thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Age bottom boundary layer", "0e0",
                             Patterns::Double (0),
                             "The age of the lower thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Radius", "0e0",
                             Patterns::Double (0),
                             "The Radius (in m) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain.");
          prm.declare_entry ("Amplitude", "0e0",
                             Patterns::Double (0),
                             "The amplitude (in K) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain. This perturbation will be added to "
                             "the adiabatic temperature profile, but not to the bottom thermal "
                             "boundary layer. Instead, the maximum of the perturbation and the bottom "
                             "boundary layer temperature will be used.");
          prm.declare_entry ("Position", "center",
                             Patterns::Selection ("center|coordinates"),
                             "Where the initial temperature perturbation should be placed. If 'center' is "
                             "given, then the perturbation will be centered along a 'midpoint' of some "
                             "sort of the bottom boundary. For example, in the case of a box geometry, "
                             "this is the center of the bottom face; in the case of a spherical shell "
                             "geometry, it is along the inner surface halfway between the bounding "
                             "radial lines. If 'coordinates' is given, then the perturbation "
                             "will be centered around the given coordinates in 'Midpoint coordinates'.");
          prm.declare_entry ("Midpoint coordinates", "",
                             Patterns::List(Patterns::Double()),
                             "Where the initial temperature perturbation should be placed. The list is "
                             "enforced to have dim entries if 'coordinates is selected for 'Position' "
                             "nd is interpreted to be the cartesian coordinates of the center of the "
                             "anomaly.");
          prm.declare_entry ("Subadiabaticity", "0e0",
                             Patterns::Double (0),
                             "If this value is larger than 0, the initial temperature profile will "
                             "not be adiabatic, but subadiabatic. This value gives the maximal "
                             "deviation from adiabaticity. Set to 0 for an adiabatic temperature "
                             "profile. Units: K.\n\n"
                             "The function object in the Function subsection "
                             "represents the compositional fields that will be used as a reference "
                             "profile for calculating the thermal diffusivity. "
                             "This function is one-dimensional and depends only on depth. The format of this "
                             "functions follows the syntax understood by the "
                             "muparser library, see Section~\\ref{sec:muparser-format}.");
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/initial-conditions/plume/test/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Data file name",
                             "lithosphere_test",
                             Patterns::Anything (),
                             "The file name of the lithosphere thickness data.");
          prm.declare_entry ("Use lithosphere thickness file",
                             "false",
                             Patterns::Bool (),
                             "Determines whether the plugin attempts to load a lithosphere "
                             "thickness file with the same file format that the ascii data plugin "
                             "expects.");

          prm.enter_subsection("Function");
          {
            Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      // we need to get the number of compositional fields here to
      // initialize the function parser. unfortunately, we can't get it
      // via SimulatorAccess from the simulator itself because at the
      // current point the SimulatorAccess hasn't been initialized
      // yet. so get it from the parameter file directly.
      prm.enter_subsection ("Compositional fields");
      const unsigned int n_compositional_fields = prm.get_integer ("Number of fields");
      prm.leave_subsection ();

      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Plume");
        {
          age_top_boundary_layer = prm.get_double ("Age top boundary layer");
          age_bottom_boundary_layer = prm.get_double ("Age bottom boundary layer");
          radius = prm.get_double ("Radius");
          amplitude = prm.get_double ("Amplitude");
          perturbation_position = prm.get("Position");
          if (perturbation_position == "coordinates")
            {
              std::vector<double> coordinates = dealii::Utilities::string_to_double(
                                                  dealii::Utilities::split_string_list(prm.get("Midpoint coordinates")));

              for (unsigned int i = 0; i < dim; i++)
                midpoint_coordinates(i) = coordinates[i];
            }

          subadiabaticity = prm.get_double ("Subadiabaticity");

          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory    = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = data_directory.find (subst_text),  position!=std::string::npos)
              data_directory.replace (data_directory.begin()+position,
                                      data_directory.begin()+position+subst_text.size(),
                                      ASPECT_SOURCE_DIR);
          }

          data_file_name    = prm.get ("Data file name");
          use_lithosphere_file    = prm.get_bool ("Use lithosphere thickness file");

          if (n_compositional_fields > 0)
            {
              prm.enter_subsection("Function");
              try
                {
                  function.reset (new Functions::ParsedFunction<1>(n_compositional_fields));
                  function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t<Initial conditions/Plume/Function>\n"
                            << "with expression\n"
                            << "\t<" << prm.get("Function expression") << ">";
                  throw;
                }

              prm.leave_subsection();
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(Plume,
                                       "plume",
                                       "Temperature is prescribed as an adiabatic "
                                       "profile with upper and lower thermal boundary layers, "
                                       "whose ages are given as input parameters or as a data file.")
  }
}

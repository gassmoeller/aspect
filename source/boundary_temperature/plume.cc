/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/boundary_temperature/plume.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/parameter_handler.h>
#include <fstream>
#include <iostream>
#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {

    namespace internal
    {
      template <int dim>
      PlumeLookup<dim>::PlumeLookup(const std::string &filename,
                                    const ConditionalOStream &pcout)
      {
        pcout << std::endl << "   Loading Plume position file for boundary temperature: "
              << filename << "." << std::endl << std::endl;

        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (fexists(filename),
                     ExcMessage (std::string("Plume position file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        const double Myr_in_seconds = 1e6 * year_in_seconds;
        const double km_in_m = 1e3;

        std::string temp;
        std::ifstream in(filename.c_str(), std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't open file <") + filename));

        double time,start_time,x,y;

        while (in >> time >> x >> y)
          {
            Point<dim> position;
            switch (dim)
              {
                case 2:
                  position(0) = x;
                  break;
                case 3:
                  position(0) = x;
                  position(1) = y;
                  break;
                default:
                  AssertThrow(false,ExcNotImplemented());
                  break;
              }
            plume_positions.push_back(position*km_in_m);

            if (times.size() == 0)
              start_time = time * Myr_in_seconds;
            times.push_back(start_time - time * Myr_in_seconds);
          }
      }

      template <int dim>
      bool
      PlumeLookup<dim>::fexists(const std::string &filename) const
      {
        std::ifstream ifile(filename.c_str());
        return !(!ifile); // only in c++11 you can convert to bool directly
      }

      template <int dim>
      Point<dim>
      PlumeLookup<dim>::plume_position(const double time) const
      {
        if (time <= times.front())
          return plume_positions.front();
        else if (time >= times.back())
          return plume_positions.back();
        else
          {
            for (unsigned int i = 0; i < times.size() - 1; i++)
              {
                if ((time >= times[i])
                    && (time < times[i+1]))
                  {
                    const double timestep = times[i+1]-times[i];
                    const double time_weight = (time - times[i]) / timestep;

                    return Point<dim> ((1.0 - time_weight) * plume_positions[i]
                                       + time_weight * plume_positions[i+1]);
                  }
              }
            AssertThrow(false,
                        ExcMessage("Did not find time interval for plume location."))
          }
        return Point<dim>();
      }
    }

// ------------------------------ Box -------------------

    template <int dim>
    void
    Plume<dim>::
    initialize ()
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model())
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      lookup.reset(new internal::PlumeLookup<dim>(data_directory+plume_file_name,
                                                  this->get_pcout()));

      if (use_lithosphere_thickness_files)
        {
          const std::set<types::boundary_id> boundary_ids = this->get_fixed_temperature_boundary_indicators();
          lithosphere_thickness_file->initialize(boundary_ids, 1);
        }
    }

    template <int dim>
    void
    Plume<dim>::
    update ()
    {
      plume_position = lookup->plume_position(this->get_time()
                                              - model_time_to_start_plume_tail);

      if (use_lithosphere_thickness_files)
        {
          lithosphere_thickness_file->update();
        }
    }

    template <int dim>
    Point<dim>
    Plume<dim>::
    get_plume_position () const
    {
      return plume_position;
    }

    template <int dim>
    double
    Plume<dim>::
    adiabatic_temperature (const Point<dim> &position,
                           const unsigned int boundary_indicator) const
    {
      const unsigned int bottom = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");
      const unsigned int top = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      // then, get the temperature of the adiabatic profile at a representative
      // point at the top and bottom boundary of the model
      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const Point<dim> bottom_point = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const double adiabatic_surface_temperature = this->get_adiabatic_conditions().temperature(surface_point);
      const double adiabatic_bottom_temperature = this->get_adiabatic_conditions().temperature(bottom_point);

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
      in.strain_rate.resize(0); // adiabat has strain=0.
      this->get_material_model().evaluate(in, out);

      const double kappa = out.thermal_conductivities[0] / (out.densities[0] * out.specific_heat[0]);

      double lithosphere_age = age_top_boundary_layer;

      if (use_lithosphere_thickness_files)
        {
          const types::boundary_id boundary_id(boundary_indicator);
          const double lithosphere_thickness = lithosphere_thickness_file->get_data_component(boundary_id, position, 0);
          // conversion from lithosphere thickness in lithosphere age as required for half-space cooling model equation
          // see Turcotte & Schubert (2002) "Geodynamics", p. 287, equation 4.126
          lithosphere_age = lithosphere_thickness * lithosphere_thickness / (2.32 * 2.32 * kappa);
        }

      // analytical solution for the thermal boundary layer from half-space cooling model
      const double surface_cooling_temperature = lithosphere_age > 0.0 ?
                                                 (temperature_[top] - adiabatic_surface_temperature) *
                                                 erfc(this->get_geometry_model().depth(position) /
                                                      (2 * sqrt(kappa * lithosphere_age)))
                                                 : 0.0;

      const double bottom_heating_temperature = age_bottom_boundary_layer > 0.0 ?
                                                (temperature_[bottom] - adiabatic_bottom_temperature + subadiabaticity)
                                                * erfc((this->get_geometry_model().maximal_depth()
                                                        - this->get_geometry_model().depth(position)) /
                                                       (2 * sqrt(kappa * age_bottom_boundary_layer)))
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
             + bottom_heating_temperature + subadiabatic_T;
    }

    template <int dim>
    double
    Plume<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &position) const
    {

      Assert (boundary_indicator<2*dim, ExcMessage ("Unknown boundary indicator."));
      Assert (this->get_adiabatic_conditions().is_initialized(),
              ExcMessage ("The adiabatic conditions are not yet initialized,"
                          "but they are necessary for the plume boundary temperature plugin."));

      double boundary_temperature(0);
      const types::boundary_id boundary_id(boundary_indicator);

      if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "bottom")
        {
          const Point<dim> distance = position - plume_position;

          double distance_head_to_boundary,current_head_radius(0);
          if ((head_radius > 0) && (head_velocity > 0) && (head_amplitude > 0))
            {
              distance_head_to_boundary = fabs(head_velocity * (this->get_time() - model_time_to_start_plume_tail));

              // If the plume is not yet there, perturbation will not be set
              if (distance_head_to_boundary < head_radius)
                {
                  current_head_radius = sqrt(head_radius * head_radius
                                             - distance_head_to_boundary * distance_head_to_boundary);
                }
            }

          //Normal plume tail
          if ((this->get_time() >= model_time_to_start_plume_tail)
              && (current_head_radius < tail_radius))
            {
              // T=T_0*exp-(r/r_0)**2
              boundary_temperature += tail_amplitude * std::exp(-std::pow(distance.norm()/tail_radius,2));
            }
          // Plume head
          else
            {
              if (distance.norm() < current_head_radius)
                {
                  boundary_temperature += head_amplitude;
                  //const double current_head_amplitude = head_amplitude * std::exp(-std::pow(distance_head_to_boundary/head_radius,2));
                  //boundary_temperature += current_head_amplitude * std::exp(-std::pow(distance.norm()/current_head_radius,2));
                }
            }

          boundary_temperature += this->get_adiabatic_conditions().temperature(position);
        }
      else if (this->get_geometry_model().translate_id_to_symbol_name(boundary_id) == "top")
        boundary_temperature = temperature_[boundary_indicator];
      else
        {
          //TODO: initial is not working properly right now, it returns 0 at the surface
          if (side_boundary_type == initial)
            boundary_temperature = this->get_initial_conditions().initial_temperature(position);
          else if (side_boundary_type == constant)
            boundary_temperature = temperature_[boundary_indicator];
          else if (side_boundary_type == adiabatic)
            boundary_temperature = adiabatic_temperature(position, boundary_indicator);
          else
            AssertThrow (false, ExcNotImplemented());
        }

      return boundary_temperature;
    }


    template <int dim>
    double
    Plume<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::min_element(temperature_, temperature_+2*dim);
      else
        {
          double min = maximal_temperature(fixed_boundary_ids);
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              min = std::min(min,temperature_[*p]);
          return min;
        }
    }



    template <int dim>
    double
    Plume<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::max_element(temperature_, temperature_+2*dim);
      else
        {
          double max = -std::numeric_limits<double>::max();
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              max = std::max(max,temperature_[*p]);
          return max;
        }
    }

    template <int dim>
    void
    Plume<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Plume");
      {
        prm.declare_entry ("Data directory",
                           "$ASPECT_SOURCE_DIR/data/boundary-temperature/plume/",
                           Patterns::DirectoryName (),
                           "The name of a directory that contains the model data. This path "
                           "may either be absolute (if starting with a '/') or relative to "
                           "the current directory. The path may also include the special "
                           "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the 'data/' subdirectory of ASPECT. ");
        prm.declare_entry ("Plume position file name", "Tristan.sur",
                           Patterns::Anything (),
                           "The file name of the plume position data.");
        prm.declare_entry ("Amplitude", "0",
                           Patterns::Double (),
                           "Amplitude of the plume tail temperature anomaly. Units: K.");
        prm.declare_entry ("Inflow velocity", "0",
                           Patterns::Double (),
                           "Magnitude of the velocity inflow. Units: K.");
        prm.declare_entry ("Radius", "0",
                           Patterns::Double (),
                           "Radius of the plume tail temperature anomaly. Units: m.");
        prm.declare_entry ("Head amplitude", "0",
                           Patterns::Double (),
                           "Amplitude of the plume head temperature anomaly. Units: K.");
        prm.declare_entry ("Head radius", "0",
                           Patterns::Double (),
                           "Radius of the plume head temperature anomaly. Units: m.");
        prm.declare_entry ("Head velocity", "0",
                           Patterns::Double (),
                           "Magnitude of the plume head velocity inflow. Units: m/s or m/yr.");
        prm.declare_entry ("Model time to start plume tail", "0",
                           Patterns::Double (),
                           "Time before the start of the plume position data at which "
                           "the head starts to flow into the model. Units: years or "
                           "seconds.");
        /**
         * Choose the type of side boundary temperature applied. Current choices
         * are to apply temperatures from an adiabatic profile, the initial
         * temperature or constant boundary temperatures.
         */
        prm.declare_entry ("Side boundary type", "initial",
                           Patterns::Selection ("initial|adiabatic|constant"),
                           "A selection that determines the assumed temperatures "
                           "at the side boundaries. Choices are initial, adiabatic "
                           "and constant to prescribe the initial temperature, "
                           "temperatures from an adiabatic profile or a constant "
                           "temperature.");
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
        prm.declare_entry ("Left temperature", "1613",
                           Patterns::Double (),
                           "Temperature at the left boundary (at minimal x-value). Units: K.");
        prm.declare_entry ("Right temperature", "1613",
                           Patterns::Double (),
                           "Temperature at the right boundary (at maximal x-value). Units: K.");
        prm.declare_entry ("Bottom temperature", "4000",
                           Patterns::Double (),
                           "Temperature at the bottom boundary (at minimal z-value). Units: K.");
        prm.declare_entry ("Top temperature", "273",
                           Patterns::Double (),
                           "Temperature at the top boundary (at maximal x-value). Units: K.");
        if (dim==3)
          {
            prm.declare_entry ("Front temperature", "1613",
                               Patterns::Double (),
                               "Temperature at the front boundary (at minimal y-value). Units: K.");
            prm.declare_entry ("Back temperature", "1613",
                               Patterns::Double (),
                               "Temperature at the back boundary (at maximal y-value). Units: K.");
          }
        prm.declare_entry ("Use lithosphere thickness files", "false",
                           Patterns::Bool(),
                           "Whether to read in time-dependent AsciiData files that contain the "
                           "lithosphere thickness in km at the side boundaries.");

        prm.enter_subsection("Function");
        {
          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();

      prm.enter_subsection("Boundary temperature model");
      {
      // member variable lithosphere_thickness_file must be declared before it can be parsed (but does not belong to the Plume class!)
      Utilities::AsciiDataBoundary<2>::declare_parameters(prm,"$ASPECT_SOURCE_DIR/data/boundary-temperature/ascii-data/test/","lith_thick_%s_%d");
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Plume<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Plume");
      {
        // Get the path to the data files. If it contains a reference
        // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
        // as a #define
        data_directory        = prm.get ("Data directory");
        {
          const std::string      subst_text = "$ASPECT_SOURCE_DIR";
          std::string::size_type position;
          while (position = data_directory.find (subst_text),  position!=std::string::npos)
            data_directory.replace (data_directory.begin()+position,
                                    data_directory.begin()+position+subst_text.size(),
                                    ASPECT_SOURCE_DIR);
        }

        plume_file_name    = prm.get ("Plume position file name");
        tail_velocity = prm.get_double ("Inflow velocity");
        tail_amplitude = prm.get_double ("Amplitude");
        tail_radius = prm.get_double ("Radius");

        head_amplitude = prm.get_double("Head amplitude");
        head_radius = prm.get_double("Head radius");
        head_velocity = prm.get_double("Head velocity");
        model_time_to_start_plume_tail = prm.get_double ("Model time to start plume tail");

        use_lithosphere_thickness_files = prm.get_bool ("Use lithosphere thickness files");

        if (prm.get ("Side boundary type") == "initial")
          side_boundary_type = initial;
        else if (prm.get ("Side boundary type") == "adiabatic")
          side_boundary_type = adiabatic;
        else if (prm.get ("Side boundary type") == "constant")
          side_boundary_type = constant;
        else
          AssertThrow (false, ExcNotImplemented());

        age_top_boundary_layer = prm.get_double ("Age top boundary layer");
        age_bottom_boundary_layer = prm.get_double ("Age bottom boundary layer");
        // convert input ages to seconds
        if (this->convert_output_to_years())
          {
            tail_velocity /= year_in_seconds;
            head_velocity /= year_in_seconds;
            age_top_boundary_layer *= year_in_seconds;
            age_bottom_boundary_layer *= year_in_seconds;
            model_time_to_start_plume_tail *= year_in_seconds;
          }

        subadiabaticity = prm.get_double ("Subadiabaticity");
        if (this->n_compositional_fields() > 0)
          {
            prm.enter_subsection("Function");
            try
              {
                function.reset (new Functions::ParsedFunction<1>(this->n_compositional_fields()));
                function->parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t<Initial conditions/Adiabatic/Function>\n"
                          << "with expression\n"
                          << "\t<" << prm.get("Function expression") << ">";
                throw;
              }

            prm.leave_subsection();
          }

        switch (dim)
          {
            case 2:
              temperature_[0] = prm.get_double ("Left temperature");
              temperature_[1] = prm.get_double ("Right temperature");
              temperature_[2] = prm.get_double ("Bottom temperature");
              temperature_[3] = prm.get_double ("Top temperature");
              break;

            case 3:
              temperature_[0] = prm.get_double ("Left temperature");
              temperature_[1] = prm.get_double ("Right temperature");
              temperature_[2] = prm.get_double ("Front temperature");
              temperature_[3] = prm.get_double ("Back temperature");
              temperature_[4] = prm.get_double ("Bottom temperature");
              temperature_[5] = prm.get_double ("Top temperature");
              break;

            default:
              Assert (false, ExcNotImplemented());
          }
      }
      prm.leave_subsection ();

      if (use_lithosphere_thickness_files)
        {
          lithosphere_thickness_file.reset(new Utilities::AsciiDataBoundary<dim> ());
          // manually initialize the simulator here, since this does not happen automatically!
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(lithosphere_thickness_file.get()))
            sim->initialize (this->get_simulator());
          prm.enter_subsection("Boundary temperature model");
          {
          lithosphere_thickness_file->parse_parameters(prm);
          }
          prm.leave_subsection ();
        }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Plume,
                                               "plume",
                                               "A model in which the temperature is chosen constant on "
                                               "all the sides of a box.")
  }
}

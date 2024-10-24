#!/bin/bash

# Instructions for how to use this script are provided in the README.

for stokes_degree in 2; do #3
  for cfl in 1 0.5 0.25; do
    for refinement in 5; do
      for particles_per_direction in 10; do # 6 7 10 15 20
        echo "set CFL number = $cfl" > current.prm
        echo "subsection Discretization" >> current.prm
        echo "  set Stokes velocity polynomial degree = $stokes_degree" >> current.prm
        echo "end" >> current.prm

        echo "subsection Postprocess" >> current.prm
        echo "  subsection Particles" >> current.prm
        echo "    subsection Generator" >> current.prm
        echo "      subsection Reference cell" >> current.prm
        echo "        set Number of particles per cell per direction = $particles_per_direction" >> current.prm
        echo "      end" >> current.prm
        echo "    end" >> current.prm
        echo "  end" >> current.prm
        echo "end" >> current.prm

        echo "subsection Mesh refinement" >> current.prm
        echo "  set Initial global refinement = $refinement" >> current.prm
        echo "end" >> current.prm


        echo "set Output directory = Q${stokes_degree}_CFL${cfl}_refinement${refinement}_${particles_per_direction}_first_order_time" >> current.prm
        echo "Starting Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}_${particles_per_direction}"
        cat time_dependent_annulus.prm current.prm | mpirun -np 6 ../plugin/aspect --
      done
    done
  done
done

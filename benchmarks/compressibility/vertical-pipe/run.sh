#!/bin/bash

# Run all combinations of models in this folder

for resolution in 1; do #0 1 2 3
  echo "subsection Mesh refinement" > resolution.prm
  echo "set Initial global refinement = $resolution" >> resolution.prm
  echo "end" >> resolution.prm

  for temperature in adiabatic; do #sub-adiabatic adiabatic super-adiabatic
    for formulation in projected-density; do #ala isothermal hydrostatic projected-density
      output_folder=output-compressibility-resolution-${resolution}-${formulation}-${temperature}
      echo "set Output directory = ${output_folder}" > output.prm
      cat vertical-pipe.prm ${temperature}.prm ${formulation}.prm resolution.prm output.prm | mpirun -np 6 aspect --
      rm output.prm
    done
  done
done

bash make_statistics.sh

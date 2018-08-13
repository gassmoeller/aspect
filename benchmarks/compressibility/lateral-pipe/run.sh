#!/bin/bash

# Run all combinations of models in this folder

for resolution in 0 1 2 3; do
echo "subsection Mesh refinement" > resolution.prm
echo "set Initial global refinement = $resolution" >> resolution.prm
echo "end" >> resolution.prm

for temperature in sub-adiabatic adiabatic super-adiabatic; do #sub-adiabatic adiabatic super-adiabatic
  for formulation in ala isothermal hydrostatic projected-density; do #ala isothermal hydrostatic projected-density
    output_folder=output-lateral-pipe-resolution-${resolution}-${formulation}-${temperature}
    echo "set Output directory = ${output_folder}" > output.prm
    cat lateral-pipe.prm ${temperature}.prm ${formulation}.prm resolution.prm output.prm | mpirun -np 6 aspect --
    rm output.prm
  done
done
done

bash make_statistics.sh

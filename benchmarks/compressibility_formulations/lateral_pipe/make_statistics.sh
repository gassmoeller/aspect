#!/bin/bash

for temperature in sub_adiabatic adiabatic; do
  outputfile=mass_flux_error_${temperature}
  echo -n "# Resolution ala isentropic hydrostatic projected-density" > $outputfile
  for repetitions in 2 4 8 16 32 64 128; do
    echo '' >> $outputfile
    echo -n $repetitions >> $outputfile
    for formulation in ala isentropic hydrostatic projected_density; do
      data_folder=output-lateral_pipe_repetitions_${repetitions}_${formulation}_${temperature}
      left_flux=`grep "Outward mass flux through boundary with indicator 0" $data_folder/statistics | gawk '{print $2}' | sed s/.$//`
      echo $left_flux
      right_flux=`grep "Outward mass flux through boundary with indicator 1" $data_folder/statistics | gawk '{print $2}' | sed s/.$//`
      cat $data_folder/statistics | tail -n 2 | head -n 1 | gawk "{printf \" %g\",sqrt((\$${left_flux}+\$${right_flux})*(\$${left_flux}+\$${right_flux}))/\$${right_flux}}" >> $outputfile
    done
  done
done

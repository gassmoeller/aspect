#!/bin/bash

for temperature in sub_adiabatic adiabatic; do
  outputfile=mass_flux_error_${temperature}
  echo -n "# Resolution ala isentropic hydrostatic projected_density" > $outputfile
  for repetitions in 2 4 8 16 32 64 128; do
    echo '' >> $outputfile
    echo -n $repetitions >> $outputfile
    for formulation in ala isentropic hydrostatic projected_density; do
      data_folder=output-vertical_pipe_repetitions_${repetitions}_${formulation}_${temperature}
      bottom_flux=`grep "Outward mass flux through boundary with indicator 2" $data_folder/statistics | gawk '{print $2}' | sed s/.$//`
      top_flux=`grep "Outward mass flux through boundary with indicator 3" $data_folder/statistics | gawk '{print $2}' | sed s/.$//`
      cat $data_folder/statistics | tail -n 1 | gawk "{printf \" %g\",sqrt((\$${bottom_flux}+\$${top_flux})*(\$${bottom_flux}+\$${top_flux}))/sqrt(\$${top_flux}*\$${top_flux})}" >> $outputfile
    done
  done
  echo '' >> $outputfile
done



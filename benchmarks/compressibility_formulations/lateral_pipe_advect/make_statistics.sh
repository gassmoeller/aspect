#!/bin/bash

outputfile=mass_flux_error
echo -n "# Resolution ala isentropic hydrostatic projected-density" > $outputfile
for repetitions in 2 4 8 16 32 64 128; do
  echo '' >> $outputfile
  echo -n $repetitions >> $outputfile
  for formulation in ala isentropic hydrostatic projected_density; do
    data_folder=output-lateral_pipe_advect_repetitions_${repetitions}_${formulation}
    echo $data_folder
    right_flux=`grep "Outward mass flux through boundary with indicator 1" $data_folder/statistics | gawk '{print $2}' | sed s/.$//`
    cat $data_folder/statistics | tail -n 1 | gawk "{printf \" %g\",sqrt((\$${right_flux}-71410462.65337557)*(\$${right_flux}-71410462.65337557))/(71410462.65337557)}" >> $outputfile
  done
done

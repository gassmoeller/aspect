#!/bin/bash

for temperature in sub-adiabatic adiabatic super-adiabatic; do
	outputfile=mass_flux_error_${temperature}
	echo -n "# Resolution ala isothermal hydrostatic projected-density" > $outputfile
	  for resolution in 0 1 2 3; do
		  nodes=`echo "scale=10; 10*2^$resolution" | bc -l`
		  echo '' >> $outputfile
		  echo -n $nodes >> $outputfile
  for formulation in ala isothermal hydrostatic projected-density; do #projected-density
    data_folder=output-lateral-pipe-resolution-${resolution}-${formulation}-${temperature}
    cat $data_folder/statistics | tail -n 1 | gawk '{printf " %g",sqrt(($21+$22)*($21+$22))/$22}' >> $outputfile
    done
  done
done

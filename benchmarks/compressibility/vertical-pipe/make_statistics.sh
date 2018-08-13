#!/bin/bash

for temperature in sub-adiabatic adiabatic super-adiabatic; do
	outputfile=mass_flux_error_${temperature}
	echo -n "Resolution ala isothermal hydrostatic projected-density" > $outputfile
	  for resolution in 0 1 2 3; do
		  echo '' >> $outputfile
		  echo -n $resolution >> $outputfile
  for formulation in ala isothermal hydrostatic projected-density; do #projected-density
    data_folder=output-compressibility-resolution-${resolution}-${formulation}-${temperature}
    cat $data_folder/statistics | tail -n 1 | gawk '{printf " %g",sqrt(($24+$25)*($24+$25))}' >> $outputfile
    done
  done
done

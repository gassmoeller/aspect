#!/bin/bash

# Instructions for how to use this scipt are provided in the README.

processes=8
ASPECT_EXEC="../plugin/aspect"

for refinement in 6 7; do #0 1 2 3 4 5
  for k in 0 1 2 4 8; do
    echo "subsection Mesh refinement" > current.prm
    echo "  set Initial global refinement = $refinement" >> current.prm
    echo "end" >> current.prm

    echo "subsection Annulus benchmark" >> current.prm
    echo "  set k = $k" >> current.prm
    echo "end" >> current.prm

    echo "set Output directory = output-refinement_${refinement}_k_${k}" >> current.prm

    cat annulus.prm current.prm | mpirun -np $processes $ASPECT_EXEC --
  done
done

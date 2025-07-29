#!/bin/bash

for averaging in "harmonic average"; do # arithmetic/geometric/harmonic average
  for nsinkers in 4 16; do
    for viscosity in 1e6; do
      for refinement in 3 4; do # 3
        for solver in "nsinker_bfbt" "nsinker_gmg_gc" ; do # "nsinker" "nsinker_gmg" "nsinker_direct"
          echo "subsection Material model" > current.prm
          echo "  set Material averaging = $averaging" >> current.prm
          echo "  subsection NSinker" >> current.prm
          echo "    set Number of sinkers = $nsinkers" >> current.prm
          echo "    set Dynamic viscosity ratio = $viscosity" >> current.prm
          echo "  end" >> current.prm
          echo "end" >> current.prm

          echo "subsection Mesh refinement" >> current.prm
          echo "  set Initial global refinement = $refinement" >> current.prm
          echo "  set Initial adaptive refinement = 0" >> current.prm
          echo "end" >> current.prm

          current_model="${solver}_averaging${averaging// /_}_nsinkers${nsinkers}_viscosity${viscosity}_refinement${refinement}"
          echo "set Output directory = output-${current_model}" >> current.prm
          echo "Starting ${current_model}"
          cat ${solver}.prm current.prm | mpirun -np 8 ./aspect-release --
        done
      done
    done
  done
done

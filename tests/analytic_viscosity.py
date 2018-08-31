#!/bin/python

import numpy as np

dislocation_activation_energy = 5.3e5
dislocation_activation_volume = 1.4e-5
dislocation_creep_prefactor = 1.1e-16
adiabatic_pressure = 0.0

dislocation_creep_exponent = 3.5
gas_constant = 8.3144621
temperature = 1600.0


second_strain_rate_invariant = 6.33775e-14

second_strain_rate_invariant = 6.2839e-14

energy_term = np.exp((dislocation_activation_energy + dislocation_activation_volume * adiabatic_pressure) / (dislocation_creep_exponent * gas_constant * temperature))

strain_rate_dependence = (1.0 - dislocation_creep_exponent) / dislocation_creep_exponent;

dislocation_viscosity = dislocation_creep_prefactor**(-1.0/dislocation_creep_exponent) * second_strain_rate_invariant**strain_rate_dependence * energy_term

print (dislocation_viscosity)

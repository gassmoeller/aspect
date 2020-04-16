#!/bin/perl
#
# This script uses the new needs_viscosity_evaluation function
# of the material model input object

while (<>)
{
  s/in.strain_rate.size\(\) > 0/in.needs_viscosity_evaluation\(\)/g;
  s/in.strain_rate.size\(\)/in.needs_viscosity_evaluation\(\)/g;
  print;
}

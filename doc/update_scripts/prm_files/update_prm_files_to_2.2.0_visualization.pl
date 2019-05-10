#!/bin/perl
#
# This script unifies old visualization postprocessor names
# into the newer and more efficient 'material properties'
# postprocessor.

while (<>)
{
  if (/List of output variables/)
  {
  if (/, material properties/ == 0)
  {
  s/density/material properties/g;
  }
  else
  {
    s/density/material properties/g;
  if (/material properties/ == 0)
  {
  s/thermal expansivity/material properties/g;
  }
  s/thermal conductivity/material properties/g;
  s/thermal diffusivity/material properties/g;
  s/specific heat/material properties/g;
  s/viscosity/material properties/g;
  }
  print;
}

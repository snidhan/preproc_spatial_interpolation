#!/bin/bash
ulimit -s unlimited
#ulimit -s 999999999
rm resinterp.x #*.vtk *.dat

ifort -O3 -mcmodel=large -debug -traceback -o resinterp.x resinterp_xyzz_2.f90
#ifort -O3 -mcmodel=large -o resinterp.x resinterp_xyzz_2.f90
./resinterp.x

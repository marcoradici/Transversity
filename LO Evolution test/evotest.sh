#!/bin/bash

FF=gfortran
#CXX=clang++
CXX=g++

#$FF  -O3 -c DriverF.f90
$FF  -O3 -c evotest.f90
$FF  -O3 -c u_valence.f90
$FF  -O3 -c d_valence.f90

$CXX -O3 -fPIC -std=c++11 `apfelxx-config --cppflags` -c Evolution.cc

#$CXX -O3 -c DriverC.cc
#$CXX -o DriverC d_valence.o u_valence.o Evolution.o DriverC.o `apfelxx-config --ldflags`

#$FF  -o DriverF d_valence.o u_valence.o Evolution.o DriverF.o -lc++ -lstdc++ `apfelxx-config --ldflags`
$FF  -o evotest.exe d_valence.o u_valence.o Evolution.o evotest.o -lc++ -lstdc++ `apfelxx-config --ldflags`


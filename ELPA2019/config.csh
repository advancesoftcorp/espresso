#!/bin/csh -f

set FC = mpiifort
set CC = mpiicc
set CXX = mpiicpc
set MKL_HOME = /opt/app/intel2019/mkl
set SCALAPACK = "-L${MKL_HOME}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -Wl,-rpath,${MKL_HOME}/lib/intel64"

./configure --prefix=${PWD} FC=${FC} CC=${CC} CXX=${CXX} SCALAPACK_LDFLAGS="${SCALAPACK}" SCALAPACK_FCFLAGS="${SCALAPACK}"


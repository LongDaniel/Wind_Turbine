#!/bin/bash
#PBS -N 339
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=32:mpiprocs=32:host=n1
#PBS -j oe
#PBS -m bae
#PBS -S /bin/bash

cd $PBS_O_WORKDIR

#mpirun -np 32 ./Channel_DNS.e > result
cp fort.11.second fort.11
mpirun -np 32 ./fst > err1

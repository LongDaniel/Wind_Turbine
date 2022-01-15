#!/bin/bash
### Job name
#PBS -N clipper_tsr5 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=3:ppn=12
#PBS -l walltime=47:59:00
#PBS -j oe

cd $PBS_O_WORKDIR

cp fort.11.second fort.11

mpirun -np 32 ./fst > err1 2>&1

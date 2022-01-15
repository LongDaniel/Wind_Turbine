#!/bin/bash
### Job name
#PBS -N SAFL_mini 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=2:ppn=16
#PBS -l walltime=119:59:00
#PBS -j oe

cd $PBS_O_WORKDIR

cp fort.11.second fort.11

mpirun -np 32 ./fst > err1 2>&1

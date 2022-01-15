#!/bin/bash
### Job name
#PBS -N inflow06_1_pre 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=4:ppn=16
#PBS -l walltime=0:30:00
#PBS -j oe

cd /$PBS_O_WORKDIR
mpirun ./fst > err1 2>&1

#!/bin/bash
### Job name
#PBS -N move05_pre 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=16:ppn=16
#PBS -l walltime=4:59:00
#PBS -j oe

cd /$PBS_O_WORKDIR
mpirun ./main_hos > err1_hos 2>&1
mpirun ./fst > err1 2>&1

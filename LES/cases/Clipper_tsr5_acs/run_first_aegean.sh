#!/bin/bash
### Job name
#PBS -N admr06_pre 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=8:ppn=16
#PBS -l walltime=0:30:00
#PBS -j oe

cd /$PBS_O_WORKDIR
mpirun ./fst > err1 2>&1

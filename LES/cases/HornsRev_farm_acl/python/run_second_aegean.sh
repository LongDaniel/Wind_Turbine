#!/bin/bash
### Job name
#PBS -N move05_1 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=16:ppn=16
#PBS -l walltime=71:59:00
#PBS -j oe

cd /$PBS_O_WORKDIR

cp fort.11.second fort.11

mpirun ./fst > err1 2>&1

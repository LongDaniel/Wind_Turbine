#!/bin/bash
### Job name
#PBS -N inflow05_post 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=4:ppn=10
#PBS -l walltime=47:59:00
#PBS -j oe

cd /$PBS_O_WORKDIR



mpirun python post_inlet_6.py > log.post_inlet 2>&1 

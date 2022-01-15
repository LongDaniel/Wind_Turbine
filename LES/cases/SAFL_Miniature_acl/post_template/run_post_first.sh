#!/bin/bash
### Job name
#PBS -N val_post 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=3:ppn=12
#PBS -l walltime=8:59:00
#PBS -j oe

cd $PBS_O_WORKDIR

cp fort.11.second fort.11

mpirun -np 32 ./fst > log.post.first 2>&1

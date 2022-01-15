#!/bin/bash
### Job name
#PBS -N motion05_post 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=4:ppn=16
#PBS -l walltime=3:59:00
#PBS -j oe

cd /$PBS_O_WORKDIR

cp fort.11.second fort.11

mpirun ./fst > log.post.first 2>&1

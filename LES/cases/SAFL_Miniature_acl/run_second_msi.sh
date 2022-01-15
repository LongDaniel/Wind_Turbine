#!/bin/bash -l
#PBS -N SAFL_mini 
#PBS -m abe
#PBS -l walltime=23:59:00,nodes=2:ppn=24,pmem=1000mb
#PBS -M plyu@umn.edu

cd $PBS_O_WORKDIR

#source ~/env_hosles.source

cp fort.11.second fort.11

mpirun -np 32 ./fst > err1 2>&1
#/safl/software/aegean/openmpi/1.5.5/gcc/4.7.0-fix/bin/mpirun --bind-to-core ./vwis > err1

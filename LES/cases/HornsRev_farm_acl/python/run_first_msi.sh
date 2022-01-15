#!/bin/bash -l
#PBS -N hosles01_1_pre 
#PBS -m abe
#PBS -l walltime=0:30:00,nodes=13:ppn=20,pmem=1000mb
#PBS -M plyu@umn.edu

cd /$PBS_O_WORKDIR

#source ~/env_hosles.source

cp fort.11.first fort.11
cp fort.12.first fort.12

mpirun -np 256 ./main_hos > err1_hos 2>&1
mpirun -np 256 ./fst > err1 2>&1
#/safl/software/aegean/openmpi/1.5.5/gcc/4.7.0-fix/bin/mpirun --bind-to-core ./vwis > err1

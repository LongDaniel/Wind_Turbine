#!/bin/bash
### Job name
#PBS -N post_python1 
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=5:ppn=10
#PBS -l walltime=3:59:00
#PBS -j oe

#cd /$PBS_O_WORKDIR

for i in "0001" "0005" "0008"
do
  echo POST_U_1D1_$i
  cd POST_U_1D1_$i
  cp ../post_u_all.py ./
  mpirun python post_u_all.py > log.post_python1 2>&1 
  cd ..
done

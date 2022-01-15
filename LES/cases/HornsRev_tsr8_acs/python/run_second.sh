#!/bin/bash
### Job name
#PS -N Turbine_HOS_LES
### Mail to user
#PBS -m ae
#PBS -k o
#PBS -l nodes=32:ppn=16
#PBS -l walltime=71:59:00
#PBS -j oe

##module load acml/4.4.0  2015/petsc/gcc-4.9.2/openmpi-1.5.5/hypre-2.8.0b/3.1-rebuild
##module load acml/4.4.0 2015/petsc/gcc-4.9.2/openmpi-1.8.5/3.4.5 
#module load 2015/petsc/gcc-4.9.2/openmpi-1.8.5/hypre-2.10.0b/3.5.3

#source /home/plyu/env_hosnew.source
#echo $LD_LIBRARY_PATH > env_debuginfo 2>&1
cd /$PBS_O_WORKDIR

cp fort.11.second fort.11

#pwd >dir_debuginfo 2>&1
#cp /home/plyu/FotisCode/testt ./
mpirun ./fst > err1 2>&1
#/safl/software/aegean/openmpi/1.5.5/gcc/4.7.0-fix/bin/mpirun --bind-to-core ./vwis > err1

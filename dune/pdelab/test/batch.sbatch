#!/bin/bash
#SBATCH -N 1
#SBATCH -p haswell10c
#SBATCH -n 20

ml llvm/6.0.0
ml openblas
ml openmpi
ml arpack-ng
ml cmake
ml superlu
ml suite-sparse
ml python/3.6.5
ml eigen


#SRUNOPT="--cpu_bind=verbose,core"


export=OMP_NUM_THREADS=1
#--export=OMP_NUM_THREADS=1

srun $SRUNOPT ./testgeneo $1

#!/bin/bash -l
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=hpib
#SBATCH --account=slchen
export OMP_NUM_THREADS=16 

/project/chenyongtin/tools/deepmd-kit/bin/lmp -i input.lammps >lmp.out


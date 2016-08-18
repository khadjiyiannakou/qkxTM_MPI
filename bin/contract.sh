#!/bin/bash -l
#SBATCH --job-name="16x32"
#SBATCH --nodes=1

#SBATCH --time=00:10:00
#SBATCH --output=/users/krikitos/qkxTM_MPI/bin/out.out
#SBATCH --error=/users/krikitos/qkxTM_MPI/bin/out.err

#======START=============================== 
echo "On which nodes it executes"
echo $SLURM_JOB_NODELIST
echo "Now run the MPI tasks..."

date
module switch PrgEnv-cray PrgEnv-intel
export OMP_NUM_THREADS=8

aprun -n 1 -d 8 /users/krikitos/qkxTM_MPI/bin/twop-hadrons 16 16 16 32 1 1 1 1 /users/krikitos/scratch/16x32/confs/conf_sm.1505 /users/krikitos/scratch/16x32/props/prop_1505_SL_up /users/krikitos/scratch/16x32/props/prop_1505_SL_down /users/krikitos/scratch/16x32/props/twop_final 0 0 0 0 50 4.0 0 

date


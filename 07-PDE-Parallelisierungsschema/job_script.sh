#!/bin/bash -l

#SBATCH --nodes=4
#SBATCH --tasks=4
#SBATCH --partition=west
#SBATCH --job-name="MPI-T"
#SBATCH --open-mode=append
#SBATCH --output=test.out

echo -n '' >> test.out
mpiexec ./circle 4
echo '--------------------------' >> test.out
#echo 'fertig' > job_script.out

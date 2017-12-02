#!/bin/bash -l

#SBATCH --nodes=5
#SBATCH --tasks=5
#SBATCH --partition=west
#SBATCH --job-name="MPI-T"
#SBATCH --open-mode=append
#SBATCH --output=test.out

echo -n '' >> test.out
mpiexec ./timempi
echo '...........................' >> test.out
mpiexec ./timempi2
echo '--------------------------' >> test.out
#echo 'fertig' > job_script.out

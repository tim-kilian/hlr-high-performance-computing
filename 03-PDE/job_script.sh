#!/bin/bash -l

#SBATCH --nodes=4
#SBATCH --partition=west
#SBATCH --ntasks=16
#SBATCH --job-name="time-t"
#SBATCH --open-mode=append
#SBATCH --output=timescript.out

echo -n '' >> timescript.out
srun -N 4 -n 16 -p west -J time-t --open-mode=append -o timescript.out ./timescript.sh
echo '--------------------------' >> timescript.out
echo 'fertig' > job_script.out

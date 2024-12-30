#!/usr/bin/bash
#SBATCH -A bramsona
#SBATCH -t 7-00:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -n 105
#SBATCH --mail-type=END
#SBATCH --mail-user=klaferri@purdue.edu

module load matlab
unset DISPLAY

for ia in $(seq 0 20)
do
 srun --exclusive -c5 -n1 matlab -nodisplay -nosplash -r "LagThicknessInput($ia)" &

done

wait

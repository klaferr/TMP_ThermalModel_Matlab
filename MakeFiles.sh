#!/usr/bin/bash
#SBATCH -A bramsona
#SBATCH -t 0-2:00:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH -n 20
#SBATCH --mail-type=END
#SBATCH --mail-user=klaferri@purdue.edu

module load matlab
unset DISPLAY

for i in $(seq 0 20)
do 
srun --exclusive -c1 -n1  matlab -nodisplay -nosplash -r "MakeParamFiles(2, 6, $i)" &
done
wait


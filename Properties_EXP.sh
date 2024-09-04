#!/bin/bash -l

#SBATCH --job-name fugacity_calcs
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu 2000MB
#SBATCH --constraint OS8
#SBATCH --mail-user gustavo.rodriguez.manotas@tuhh.de
#SBATCH --mail-type=FAIL

conda activate /fibus/fs0/02/cgr7735/envs/Fugacity_calcs #activate env


. /etc/profile.d/module.sh
module load anaconda/2023.07-1

echo -e "\nStart time:" ; date

echo "run python"
python Liq_den.py "COMPONENT_NAME" 
echo "finish python"

echo
printenv | grep SLURM
echo
scontrol show job $SLURM_JOBID

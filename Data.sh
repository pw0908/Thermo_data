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

. /etc/profile.d/module.sh
module load anaconda/2023.07-1
module load julia/1.9.2


echo -e "\nStart time:" ; date

conda activate /fibus/fs0/02/cgr7735/envs/Fugacity_calcs #activate env

echo "run julia"
julia EOS_evaluator.jl "COMPONENT_NAME" "EOS"
echo "finish julia"

echo
printenv | grep SLURM
echo
scontrol show job $SLURM_JOBID

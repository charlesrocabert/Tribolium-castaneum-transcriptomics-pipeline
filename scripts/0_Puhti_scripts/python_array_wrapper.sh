#!/bin/bash
# coding: utf-8

#SBATCH --job-name=tribolium_python_array_wrapper
#SBATCH --account=Project_2003847
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=128G
#SBATCH --time=10:00:00
#SBATCH --gres=nvme:100
#SBATCH --output=./jobs/output_%j.txt
#SBATCH --error=./jobs/errors_%j.txt

### Load modules ###
module load gatk
module load biokit
module load allas
module load r-env-singularity

### Reinitialize Allas session ###
source /appl/opt/allas-cli-utils/allas_conf -f -k -u rocabert -p project_2003847

### Export local scratch path to DATADIR ###
export DATADIR=$LOCAL_SCRATCH

### Wrap the python script (launch main task) ###
python3 $* -sindex $SLURM_ARRAY_TASK_ID

### Save the efficiency report ###
seff $SLURM_JOBID > ./jobs/efficiency_report_$SLURM_JOBID.txt



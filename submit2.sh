#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=48
#SBATCH -p sequana_cpu_bigmem
#SBATCH --exclusive
#SBATCH --j HP2NET
#SBATCH --time=00:45:00
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


module load anaconda3/2024.02_sequana
eval "$(conda shell.bash hook)"
CONDA_ENV="/scratch/cenapadrjsd/reiglan.lourenco/RNA-seq/conda_envs/parsl_env/"
CDIR="/scratch/cenapadrjsd/reiglan.loureco/HP2NET"
export CONDA_ENV
conda activate $CONDA_ENV
python parsl_workflow.py -s default.ini -w work.config
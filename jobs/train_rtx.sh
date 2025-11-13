#!/bin/bash
#SBATCH --job-name=t1rtx
#SBATCH --error=jobs/%j.txt
#SBATCH --output=jobs/%j.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=23:59:00
#SBATCH --gres=gpu:rtx:1
cd $HOME/repos/liuqe-ml
echo "running job $SLURM_JOB_ID"
git log -1 --pretty=format:"[%h] %s"
echo
echo "Running training script RTX..."
srun jupyter nbconvert train.ipynb --to python && python train.py && rm -rf train.py
echo "Training script completed."
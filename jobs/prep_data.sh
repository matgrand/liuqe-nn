#!/bin/bash
#SBATCH --job-name=prep 
#SBATCH --error=jobs/%j.txt
#SBATCH --output=jobs/%j.txt
#SBATCH --partition=allgroups
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=4:00:00
cd $HOME/repos/liuqe-ml

echo "Latest commit:"
git log -1 --pretty=format:"[%h] %s"
echo
echo "Running dataset preparation script..."

srun jupyter nbconvert prepare_dataset.ipynb --to python && python prepare_dataset.py && rm -rf prepare_dataset.py
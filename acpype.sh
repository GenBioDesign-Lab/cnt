#!/bin/bash
#SBATCH --job-name=NAMD_CPU
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --partition=cpus
##SBATCH --gres=gpu:a6000:1
#SBATCH --array=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

ml load anaconda/2024.10
conda activate mdpipe
cd /data01/genbiolab/mdanh/data/cnt-gaff/working

# Loop through all .mol2 files in the current directory
for mol2_file in *.mol2; do
    if [ -f "$mol2_file" ]; then
        python /data01/genbiolab/mdanh/data/acpype/run_acpype.py -i "$mol2_file" -c user
    fi
done
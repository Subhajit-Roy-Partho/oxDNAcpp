#!/bin/bash
#SBATCH -N 1                     # number of Node  
#SBATCH -n 20                    # number of tasks
#SBATCH -t 0-04:00               # wall time (D-HH:MM)
#SBATCH -o job.out              # STDOUT (%j = JobId)
#SBATCH -e job.err              # STDERR (%j = JobId)
#SBATCH --job-name="Assembly"
#SBATCH --gres=gpu:1
#SBATCH -p sulcgpu1
#SBATCH -q sulcgpu1

module load gcc/11.2.0
/scratch/sroy85/Software/oxDNA/build/bin/oxDNA_debug input_BIG

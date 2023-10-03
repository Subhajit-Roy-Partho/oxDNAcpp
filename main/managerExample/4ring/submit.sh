#!/bin/bash
#SBATCH -N 1                     # number of Nodes  
#SBATCH -n 1                     # number of tasks
#SBATCH -t 7-00:00               # wall time (D-HH:MM)
#SBATCH -o job.out              # STDOUT (%j = JobId)
#SBATCH -e job.err              # STDERR (%j = JobId)
#SBATCH --job-name="Seeding Blank"

module load gcc/11.2.0
/scratch/sroy85/Software/oxDNA/build/bin/oxDNA_debug input
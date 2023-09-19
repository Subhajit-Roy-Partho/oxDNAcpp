#!/bin/bash
#SBATCH -N 1                     # number of Nodes  
#SBATCH -n 1                     # number of tasks
#SBATCH -t 7-00:00               # wall time (D-HH:MM)
#SBATCH -o job.out              # STDOUT (%j = JobId)
#SBATCH -e job.err              # STDERR (%j = JobId)
#SBATCH -p sulcgpu3
#SBATCH -q sulcgpu1
#SBATCH --job-name="Seeding Blank"

oxDNA2 input
#!/bin/bash
#SBATCH --job-name=annotate_network_thyroid    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu     # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=150gb                     # Job memory request
#SBATCH --time=10:05:00               # Time limit hrs:min:sec
#SBATCH --output=annotate_network_urinary%j.log   # Standard output and error log

ml R

Rscript annotate_network.R thyroid_gland

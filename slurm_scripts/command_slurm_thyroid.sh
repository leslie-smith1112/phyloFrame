#!/bin/sh
#SBATCH --job-name=phyloFrame     # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Run a single task		
#SBATCH --cpus-per-task=2            # Number of CPU cores per task
#SBATCH --mem=50gb                    # Job memory request
#SBATCH --time=15:00:00              # Time limit hrs:min:sec
#SBATCH --output=phyloFrame%j.log   
pwd; hostname; date
path=$1
echo $path
penalty=$2
echo $penalty
disease=$3
echo $disease
version=$4
echo $version
bench_penalty=$5
echo $bench_penalty
var_genes=$6
echo $var_genes
sh /home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/create_dir_structure.sh $path $disease

ml R

Rscript --max-ppsize=500000 /home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/phyloFrame/phyloFrame_run.R -d $disease -m $penalty -p $path -V $version -b $bench_penalty -g $var_genes 

date

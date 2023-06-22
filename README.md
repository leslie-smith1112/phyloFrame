# phyloFrame
Workflow for ancestry unbaised genetic signatures

-- GnomAD VCF parsing done in leslie-smith1112/gnomad_vcf_parsing.  Kept as seperate directory for nor for simplicity in dependencies.

-- Prereqs for the pipeline include the annotated the parsed VCF files from Gnomad and the HumanBase network. This can be done using the script annotate_network.R. All scripts are run using slurm, example of all scripts are provided in the ./slurm_script directory. 

-- Scripts for figures used in the paper are in the ./figure_scripts directory. 

-- In addition to the GnomAD ancestry information and HumanBase tissue specific network, PhyloFrame needs the disease gene expression matrix and subtype information that should be used in the binary logistic regression prediction. Currently read ins of these files are all hardcoded. This will be updated in a future release, if you would like to run phyloFrame on a disease other than Breast, Thyroid, or Uterine from the TCGA data, I would be happy to help edit code to do that. Code for each disease that assigns subtype information to samples is in ./diseases directory. 

-- 



#!/bin/bash
#SBATCH --job-name LDprune
#SBATCH -n 1 --mem-per-cpu 5g 
#SBATCH -t 1:00:00
#SBATCH -p day,pi_zhao,scavenge
#SBATCH --mail-type=ALL --requeue

cd /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/ld_prune

geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/hapmap3/1000G_mac5eur_hapmap3_chr22_forUKB'
window_size=50  # Window size in SNPs
step_size=5  # Step size (number of SNPs to shift the window each step)
r2_threshold=0.2  # LD threshold (SNPs with r^2 > this value will be removed)

module load PLINK/1.9b_6.21-x86_64
plink --bfile ${geno} \
	--indep-pairwise ${window_size} ${step_size} ${r2_threshold} \
	--out 1000G_mac5eur_hapmap3_chr22


# MAF for indepdendent SNPs
plink --bfile ${geno} \
	--extract 1000G_mac5eur_hapmap3_chr22.prune.in \
	--freq \
	--out 1000G_mac5eur_hapmap3_chr22_pruned

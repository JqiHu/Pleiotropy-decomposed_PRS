#!/bin/bash
#SBATCH --job-name SIMU_PRS
#SBATCH -c 8 --mem-per-cpu 5g
#SBATCH -t 1:00:00
#SBATCH -p day,pi_zhao,scavenge
#SBATCH --mail-type=all --requeue

cd /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_pheno/phe1

# Function to get score
getPRS(){
  num_col=$1
#  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/hapmap3/1000G_mac5eur_hapmap3_chr22'
  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/UKB/ukb_unrelated_wb_hapmap3_chr22'
  beta='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simulated_beta_by_region.txt'
  out='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_pheno/phe1/'

  name_col=`head -n 1 ${beta} | tr '\t' '\n' | sed -n "${num_col}p"`

  module load PLINK/1.9b_6.21-x86_64 
  plink --bfile ${geno} \
	  --score ${beta} 2 5 ${num_col} sum \
	  --out ${out}simu_UKB_sum_PRS_${name_col} 
  
}
export -f getPRS

#module load parallel
#parallel getPRS {1} ::: {12..19}
for i in {15..22};do
getPRS ${i}
done


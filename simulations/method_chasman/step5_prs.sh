#!/bin/bash
#SBATCH --job-name SIMU_PRS
#SBATCH -c 4 --mem-per-cpu 5g
#SBATCH -t 1:00:00
#SBATCH -p day,pi_zhao,scavenge
#SBATCH --mail-type=all --requeue

cd /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/prs
# Function to get score
getPRS(){
  num_col=$1 # 7-11 

#  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/hapmap3/1000G_mac5eur_hapmap3_chr22'
  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/UKB/ukb_unrelated_wb_hapmap3_chr22'
  beta='../beta_component/weights_for_PRS.txt'

  name_col=`head -n 1 ${beta} | tr '\t' '\n' | sed -n "${num_col}p"`

  module load PLINK/1.9b_6.21-x86_64 
  plink --bfile ${geno} \
	  --score ${beta} 2 3 ${num_col} sum\
	  --out simu_PRS_${name_col} 
  
}
export -f getPRS

for i in 7 8 9 10 11;do
getPRS ${i}

done

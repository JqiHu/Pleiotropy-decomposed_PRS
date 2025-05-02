#!/bin/bash
#SBATCH --job-name SIMU_PRS
#SBATCH -c 8 --mem-per-cpu 5g
#SBATCH -t 1:00:00
#SBATCH -p day,pi_zhao,scavenge
#SBATCH --mail-type=all --requeue

cd /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/pd_prs
# Function to get score
getPRS(){
  num_col=$1 
  group=$2

#  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/hapmap3/1000G_mac5eur_hapmap3_chr22'
  geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/UKB/ukb_unrelated_wb_hapmap3_chr22'
  beta='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/snp_subset/simulated_beta_group'${group}'.txt'
  out='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/pd_prs/'

  name_col=`head -n 1 ${beta} | tr '\t' '\n' | sed -n "${num_col}p"`

  module load PLINK/1.9b_6.21-x86_64 
  plink --bfile ${geno} \
	  --score ${beta} 2 5 ${num_col} \
	  --out ${out}simu_PRS_${name_col}_${group} 
  
}
export -f getPRS

#module load parallel
#parallel getPRS {1} {2} ::: 12 16 ::: 1 2 3 4
for i in 15 19;do
  for j in 1 2 3 4;do
    getPRS ${i} ${j}
  done
done

#### Calculate an overall PRS as well
geno='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/UKB/ukb_unrelated_wb_hapmap3_chr22'
module load PLINK/1.9b_6.21-x86_64
plink --bfile ${geno} \
	--score /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simulated_beta_by_region.txt 2 5 15 \
	--out /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/pd_prs/simu_PRS_Beta_Trait1_all

plink --bfile ${geno} \
        --score /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simulated_beta_by_region.txt 2 5 19 \
        --out /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/pd_prs/simu_PRS_Observed_Beta_Trait1_all

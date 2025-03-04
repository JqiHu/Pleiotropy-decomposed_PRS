#!/bin/bash 
#SBATCH -c 6 
#SBATCH --job-name=ps_prs
#SBATCH --partition=scavenge,day,pi_zhao
#SBATCH --mail-type=all
#SBATCH --mem-per-cpu=30g -t 1-00:00:00
#SBATCH --requeue


module load PLINK/1.9b_6.21-x86_64
module load parallel

getScore(){
  k=$1

  cd /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data/cluster_${k}/prs_calculation
  
  parallel plink --bfile /gpfs/gibbs/pi/zhao/yy496/ukb_imp/ukb_imp_qc1_set --score ../prs_beta_selection/Cluster{1}.txt 3 4 7 --out Cluster{1} ::: $(seq 1 $k)
}
export -f getScore

parallel getScore {1} ::: 5 7 9 15


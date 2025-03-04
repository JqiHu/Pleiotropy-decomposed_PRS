# GNOVA for correlations between 43 CAD-correlated traits
# generate job file and submit it

cd /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/gnova_43 # folder for script

module load R/4.2.0-foss-2020b
Rscript gnova_43.R # generate job file 

module load dSQ
dsq --job-file gnova_43.job --job-name gnova_43 --time 12:00:00 --mem-per-cpu 10g -p general,pi_zhao,scavenge --mail-type=ALL --batch-file gnova_43.pbs --status-dir /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/dsq/



# SUPERGNOVA for CAD and 43 traits
# generate job file and submit it

cd /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/supergnova # folder for script

module load R/4.2.0-foss-2020b
Rscript job_gene.R # generate job file 

module load dSQ
dsq --job-file su43.job --job-name supergnova_43 --time 12:00:00 --mem-per-cpu 100g -p general,pi_zhao,scavenge --mail-type=ALL --batch-file su43.pbs --status-dir /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/dsq/

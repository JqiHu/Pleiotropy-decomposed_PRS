# GNOVA for correlations between CAD and 2,128 traits
# generate job file and submit it

cd /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/gnova_2128 # folder for script

module load R/4.2.0-foss-2020b
Rscript all_cad_gnova.R # generate job file 

module load dSQ
dsq --job-file all_cad_gnova.job --job-name gnova_cad_2128 --time 12:00:00 --mem-per-cpu 10g -p general,pi_zhao,scavenge --mail-type=ALL --batch-file all_cad_gnova.pbs --status-dir /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/dsq/

sbatch all_cad_gnova.pbs

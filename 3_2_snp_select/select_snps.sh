source ~/.bashrc;
module load GCC;
partition=scavenge,day,pi_zhao;
module load dSQ
dSQ --jobfile select_snps.job  -p ${partition} -n 1 --mem-per-cpu=12g -t 5:00:00 --mail-type=ALL --batch-file select_snps.pbs
sbatch select_snps.pbs

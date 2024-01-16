source ~/.bashrc;
module load GCC;
partition=pi_zhao,general,scavenge;
module load dSQ
dSQ --jobfile beta.job  -p ${partition} -n 1 -C avx2 --mem-per-cpu=20g -t 2-00:00:00 --mail-type=ALL --batch-file beta.pbs

sbatch beta.pbs

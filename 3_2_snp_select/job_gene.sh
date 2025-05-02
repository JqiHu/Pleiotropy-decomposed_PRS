# generate job file to select snps

files=`ls /gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/supergnova/res/*` # list supergnova output files

for file in $files
do
file=`echo $file`
trait=${file##*/}
trait=${trait%_*} # extract trait name

echo 'module load R/4.2; Rscript select_snps.R' $file $trait >> select_snps.job
done

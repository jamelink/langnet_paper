#!/bin/sh
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash 
#$ -e /data/clusterfs/lag/users/jitame/logs/
#$ -o /data/clusterfs/lag/users/jitame/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -N reg_1_sent_exome
#$ -m beas

# Store current date and time in variable and start time of the script
now=$( date )
start=`date +%s`
echo "Start at ${start}" 


#define paths, edit this
base=/data/clusterfs/lag/users/jitame/SENT_CORE
base_reg=${base}/geno/regenie/
in_path=$base_reg/st1_in_plink
out_path=$base_reg/exome/st1_out
pheno_file=$base/pheno/sent_edges_exome.tsv

mkdir -p $out_path/temp
cd $out_path/temp

echo "Run regenie for $pheno_file" 

#load + run regenie
module load regenie/3.2.1
regenie \
--step 1 \
--bed $in_path/sent_all_merged_regenie_step1 \
--covarFile $base/covars/covars_pc10_exome.tsv \
--catCovarList "sex,geno_array_dummy,site_dummy_11025,site_dummy_11026,site_dummy_11027,exome_batch" \
--phenoFile $pheno_file \
--bsize 1000 \
--threads 4 \
--apply-rint \
--lowmem \
--lowmem-prefix temp \
--out $out_path/reg_st1_exome

### Creating tar

out_file=$base/step_1_sent_exome.tar.gz
cd $out_path
tar cfz $out_file *.loco


# Store current date and time in variable and calculate the runtime
now=$( date )
checkpoint=`date +%s`
runtime=$(((checkpoint-start)/60))
printf "\n Elapsed time is "${runtime}" minutes.\n\n"
#!/bin/sh
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash 
#$ -e /data/clusterfs/lag/users/jitame/logs/
#$ -o /data/clusterfs/lag/users/jitame/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -m beas

usage()
{
   echo -e "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

   echo "Script for runninng regenie step 2 on local GWAS data"
   echo " "
   echo " "
   echo "Usage: bash regenie_step2_sent.sh -c <chrom_number>"
   echo -e "\t-c chromosome number - REQUIRED "
   echo " "
   echo "-- Cluster use (single15.q) --"
   echo "When running the script on the cluster (through gridmaster), you might want to provide specific qsub arguments, such as the location where a standard"
   echo "log file is stored, or the email address to which a message should be sent when the script is finished."
   echo "In this case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -N regenie_step2_sent regenie_step2_sent.sh -c <chrom_number>"
   echo -e "------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
   exit 1 # Exit script after printing usage
}

### EVALUATE SOFTWARE AND SET INPUT ###

while getopts "c:" opt
do
   case "$opt" in
	 c ) chr="$OPTARG" ;;
	 ? ) usage ;; # Print usage in case parameter is non-existent
   esac
done

# Store current date and time in variable and start time of the script
now=$( date )
start=`date +%s`
echo "Start at ${start}" 

#define paths
base=/data/clusterfs/lag/users/jitame/SENT_CORE
pred_path=$base/geno/regenie/gwas/st1_out
in_path=$base/geno/regenie/gwas/st2_in
out_path=$base/geno/regenie/gwas/st2_out/c${chr}

#Set up path and move to path
mkdir -p $out_path 
cd $out_path

#create list
awk -F " " '{print $4 }' $in_path/filtered/filter_var/subsetting_reg_st2_GWAS_chr${chr}.snpstats_mfi_hrc.snps2keep > $in_path/filtered/filter_var/subsetting_reg_st2_GWAS_chr${chr}.snplist

#load + run regenie
module load regenie/3.2.1
regenie \
--step 2 \
--bgen $in_path/imagingT1_chr${chr}.bgen \
--sample $in_path/imagingT1_chr${chr}.sample \
--covarFile $base/covars/covars_pc10_gwas.tsv \
--phenoFile $base/pheno/sent_edges_exome.tsv \
--pred $pred_path/reg_st1_gwas__pred.list \
--keep $base/exome_subs_plink.txt \
--extract $in_path/filtered/filter_var/subsetting_reg_st2_GWAS_chr${chr}.snplist \
--catCovarList "sex,geno_array_dummy,site_dummy_11025,site_dummy_11026,site_dummy_11027" \
--minMAC 297 \
--bsize 500 \
--threads 4 \
--apply-rint \
--ref-first \
--verbose \
--lowmem \
--lowmem-prefix temp \
--out $out_path/reg2_gwas_c${chr}

# Store current date and time in variable and calculate the runtime
now=$( date )
checkpoint=`date +%s`
runtime=$(((checkpoint-start)/60))
printf "\n Elapsed time is "${runtime}" minutes.\n\n"
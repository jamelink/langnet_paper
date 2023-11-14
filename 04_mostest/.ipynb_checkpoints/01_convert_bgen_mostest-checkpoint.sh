#!/bin/sh
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash 
#$ -e /data/clusterfs/lag/users/jitame/logs/
#$ -o /data/clusterfs/lag/users/jitame/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -m beas

### USAGE ###

usage()
{
   echo -e "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
   echo "MOSTEST preparataion"
   echo " "
   echo "The script takes as input a (chromosome-specific) bgen binary fileset(s), a snp-list and subject list and outputs PLINK-format (.bed, .bim .fam) and applies filtering steps to select variants to use as input in REGENIE step 1."
   echo "Filtering includes , MAF < 0.01, HWE p-value < 1e-7 and genotype missingness > 0.05 and imputation quality < 0.7."
   echo " "
   echo "Usage: bash regenie_prepare_step1_select_variants.sh -i input_file -o outdir"
   echo -e "\t-c chromosome number"
   echo " "
   echo "-- Cluster use (single15.q) --"
   echo "When running the script on the cluster (through gridmaster), you might want to provide specific qsub arguments, such as the location where a standard"
   echo "log file is stored, or the email address to which a message should be sent when the script is finished."
   echo "In this case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -e /dir/where/err/is/saved -o /dir/where/log/is/saved convert_bgen_plink_mostest.sh -i input_fileset -o outdir -n runname -s subject_list"
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


###  SPECIFY PATHS ###

Software_plink=/home/jitame/bin/software/plink2_221024
base_path=/data/clusterfs/lag/users/jitame/SENT_CORE
outdir=${base_path}/geno/mostest/in

mkdir -p ${outdir}
cd ${outdir}

### Run plink2

${Software_plink}/plink2 \
--bgen ${base_path}/geno/regenie/gwas/st2_in/imagingT1_chr${chr}.bgen 'ref-first' \
--sample ${base_path}/geno/regenie/gwas/st2_in/imagingT1_chr${chr}.sample \
--extract ${base_path}/geno/regenie/gwas/st2_in/filtered/filter_var/subsetting_reg_st2_GWAS_chr${chr}.snpstats_mfi_hrc.snps2keep \
--keep ${base_path}/exome_subs_plink.txt \
--make-bed \
--out ${outdir}/mostest_geno_in_c${chr}


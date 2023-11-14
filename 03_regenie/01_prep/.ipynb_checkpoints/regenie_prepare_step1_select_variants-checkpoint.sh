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
   echo "REGENIE STEP 1 PREPARATION - SELECT VARIANTS FROM ARRAY DATA"
   echo " "
   echo "The script takes as input a (chromosome-specific) plink binary fileset(s) (.bed/.bim/.fam) and applies filtering steps to select variants to use as input in REGENIE step 1."
   echo "Filtering includes removal of high-LD variants (> 0.9), MAF < 0.01, HWE p-value < 1e-15 and genotype missingness > 0.01."
   echo " "
   echo "Usage: bash regenie_prepare_step1_select_variants.sh -i input_file -o outdir"
   echo -e "\t-i input_bed file - REQUIRED - Full path to bed file."
   echo -e "\t-b input_bim file - REQUIRED - Full path to bim file."
   echo -e "\t-f input_fam file - REQUIRED - Full path to fam file."
   echo -e "\t-o outdir - REQUIRED - Output directory."
   echo -e "\t-n runname - REQUIRED - Run name that is appended to the output files."
   echo -e "\t-s subject_list - REQUIRED - List of subjects to include in analysis."
   echo -e "\t-h - Show help."
   echo " "
   echo "-- Cluster use (single15.q) --"
   echo "When running the script on the cluster (through gridmaster), you might want to provide specific qsub arguments, such as the location where a standard"
   echo "log file is stored, or the email address to which a message should be sent when the script is finished."
   echo "In this case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -e /dir/where/err/is/saved -o /dir/where/log/is/saved regenie_prepare_step1_select_variants.sh -i input_fileset -o outdir -n runname -s subject_list"
   echo -e "------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
   exit 1 # Exit script after printing usage
}

### EVALUATE SOFTWARE AND SET INPUT ###

while getopts "i:b:f:o:n:e:s:" opt
do
   case "$opt" in
     i ) input_bed="$OPTARG" ;;
     b ) input_bim="$OPTARG" ;;
     f ) input_fam="$OPTARG" ;;
     o ) outdir="$OPTARG" ;;
	  n ) runname="$OPTARG" ;;
     e ) exclude_variants="$OPTARG" ;;
	  s ) subject_list="$OPTARG" ;;
	  ? ) usage ;; # Print usage in case parameter is non-existent
   esac
done

Software_plink=/home/jitame/bin/software/plink19_220402
#module load plink/1.9b6

basename=$( echo "${input_bed}.bed" | awk -F"/" '{print $NF}' | awk -F"." '{print $1}' )

### SCRIPT
if [ -z "$exclude_variants" ]; then

   if [[ ${basename} =~ "_cX_" ]]; then

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --indep-pairwise 1000 100 0.9 \
      --keep ${subject_list} \
      --snps-only 'just-acgt' \
      --maf 0.01 \
      --geno 0.01 \
      --out ${outdir}/${runname}_${basename}

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --extract ${outdir}/${runname}_${basename}.prune.in \
      --keep ${subject_list} \
      --make-bed \
      --out ${outdir}/${runname}_${basename}

   else

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --indep-pairwise 1000 100 0.9 \
      --keep ${subject_list} \
      --snps-only 'just-acgt' \
      --maf 0.01 \
      --hwe 1e-15 \
      --geno 0.01 \
      --out ${outdir}/${runname}_${basename}

      ${Software_plink}/plink2 \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --extract ${outdir}/${runname}_${basename}.prune.in \
      --keep ${subject_list} \
      --make-bed \
      --out ${outdir}/${runname}_${basename}

   fi

else

   if [[ ${basename} =~ "_cX_" ]]; then

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \ 
      --indep-pairwise 1000 100 0.9 \
      --keep ${subject_list} \
      --exclude ${exclude_variants} \
      --snps-only 'just-acgt' \
      --maf 0.01 \
      --geno 0.01 \
      --out ${outdir}/${runname}_${basename}

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --extract ${outdir}/${runname}_${basename}.prune.in \
      --keep ${subject_list} \
      --make-bed \
      --out ${outdir}/${runname}_${basename}

   else

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --indep-pairwise 1000 100 0.9 \
      --keep ${subject_list} \
      --exclude ${exclude_variants} \
      --snps-only 'just-acgt' \
      --maf 0.01 \
      --hwe 1e-15 \
      --geno 0.01 \
      --out ${outdir}/${runname}_${basename}

      ${Software_plink}/plink \
      --bed ${input_bed} \
      --bim ${input_bim} \
      --fam ${input_fam} \
      --extract ${outdir}/${runname}_${basename}.prune.in \
      --keep ${subject_list} \
      --make-bed \
      --out ${outdir}/${runname}_${basename}

   fi

fi

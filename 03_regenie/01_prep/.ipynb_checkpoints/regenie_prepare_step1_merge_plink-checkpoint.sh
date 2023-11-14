#!/bin/sh
#$ -cwd
#$ -N regenie_step1_prepare_merge_plink
#$ -q single15.q
#$ -S /bin/bash

### USAGE ###

usage()
{
   echo -e "\n------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
   echo "REGENIE STEP 1 PREPARATION - MERGE PLINK BINARY FILESETS"
   echo " "
   echo "The script takes as input per-chromosome plink binary fileset(s) (.bed/.bim/.fam) stored in one directory and merges these in preparation of REGENIE step 1."
   echo " "
   echo "Usage: bash regenie_prepare_step1_merge_plink.sh -i input_file -o outdir"
   echo -e "\t-i workdir - REQUIRED - Full path to folder where per-chromosome plink binary filesets (.bed/.bim/.fam) are stored."
   echo -e "\t-n runname - REQUIRED - Run name that is appended to the output files."
   echo -e "\t-h - Show help."
   echo " "
   echo "-- Cluster use (single15.q) --"
   echo "When running the script on the cluster (through gridmaster), you might want to provide specific qsub arguments, such as the location where a standard"
   echo "log file is stored, or the email address to which a message should be sent when the script is finished."
   echo "In this case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -e /dir/where/err/is/saved -o /dir/where/log/is/saved regenie_prepare_step1_merge_plink.sh -i workdir"
   echo -e "------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
   exit 1 # Exit script after printing usage
}

### EVALUATE SOFTWARE AND SET INPUT ###

while getopts "i:n:" opt
do
   case "$opt" in
      i ) workdir="$OPTARG" ;;
      n ) runname="$OPTARG" ;;
    esac
done

Software_plink=/home/jitame/bin/software/plink19_220402

ls -l ${workdir}/*.bed | awk '{print $9}' | sort -V | awk -F".bed" '{print $1".bed\t"$1".bim\t"$1".fam"}' > ${workdir}/plink_binary_merge_list.txt

first_file=$( head -1 ${workdir}/plink_binary_merge_list.txt | awk -F".bed" '{print $1}' )

sed '1d' ${workdir}/plink_binary_merge_list.txt > ${workdir}/tmp && mv ${workdir}/tmp ${workdir}/plink_binary_merge_list.txt

# Merge all the separate chromosome files in the submitted workdir
${Software_plink}/plink \
--bed ${first_file}.bed \
--bim ${first_file}.bim \
--fam ${first_file}.fam \
--merge-list ${workdir}/plink_binary_merge_list.txt \
--make-bed \
--out ${workdir}/${runname}_merged_regenie_step1_tmp

# Apply --merge-x to the final file to also merge chr X and XY
${Software_plink}/plink \
--bfile ${workdir}/${runname}_merged_regenie_step1_tmp \
--merge-x \
--make-bed \
--out ${workdir}/${runname}_merged_regenie_step1

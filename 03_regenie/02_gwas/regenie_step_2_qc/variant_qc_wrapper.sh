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

   echo "Wrapper for variant_qc"
   echo " "
   echo "This is a wrapper for an R-pipeline for variant QC from Dick Schijven and Amaia Castillona"
   echo "This is both possible for GWAS and exome analysis."
   echo " "
   echo "Usage: bash mrp_regenie_wrapper.sh -i input_file -o outdir"
   echo -e "\t-c chromosome number - REQUIRED "
   echo " "
   echo "-- Cluster use (single15.q) --"
   echo "When running the script on the cluster (through gridmaster), you might want to provide specific qsub arguments, such as the location where a standard"
   echo "log file is stored, or the email address to which a message should be sent when the script is finished."
   echo "In this case, provide the qsub arguments before the script name, and the arguments specific to the script after the script name:"
   echo "qsub -N variant_qc variant_qc_wrapper.sh -c <chrom_number>"
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

#set environment path
module unload R/R-4.0.3
module load R/R-4.1.2

code_path=/home/jitame/bin/code/AICHA/genetics/regenie/gwas

echo "Running variant qc for chromosome ${chr}"

Rscript ${code_path}/imaging40k_variant_qc.R ${code_path}/imaging40k_variant_qc_config.R ${chr}


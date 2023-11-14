#!/bin/sh
#$ -q single15.q
#$ -S /bin/bash
#$ -e /home/jitame/bin/logs/
#$ -o /home/jitame/bin/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -m beas
#written by Jitse S. Amelink
#last update 20220210

#set up path for GCTA toolbox
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/gcta_1.93.0beta/:$PATH

#set base path
base_dir=/data/clusterfs/lag/users/jitame/SENT_CORE

while getopts "htc:f:n:c"  opt
do
   case "$opt" in
	  f ) pheno_file="$OPTARG" ;;
      n ) i="$OPTARG" ;;
      c ) output_name="$OPTARG" ;;
   esac
done

# Store current date and time in variable and start time of the script
now=$( date )
start=`date +%s`
echo "Start at ${start}" 

echo "Entered file is: $pheno_file"
echo "Number of phenotypes is: $i"
echo "Output name is: $output_name "

# GREML analysis (GCTA toolbox)

gcta64 \
--grm /data/workspaces/lag/workspaces/lg-ukbiobank/derived_data/genetic_data/snp/subset_imagingT1_40k/v1_white_british_ancestry/with_rel_filter/grm/imagingT1_allchr \
--pheno $pheno_file \
--mpheno ${i} \
--reml \
--out $output_name \
--thread-num 6



# Store current date and time in variable and calculate the runtime
now=$( date )
checkpoint=`date +%s`
runtime=$(((checkpoint-start)/60))
printf "\n Elapsed time is "${runtime}" minutes.\n\n"

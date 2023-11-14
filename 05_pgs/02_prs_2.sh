#!/bin/sh
#$ -q single15.q
#$ -S /bin/bash
#$ -e /home/jitame/bin/logs/
#$ -o /home/jitame/bin/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -m beas
#written by Jitse S. Amelink
#last update 20220207

#get options
while getopts "htc:f:n:c"  opt
do
   case "$opt" in
      f ) input_efs="$OPTARG" ;;
      n)  outp="$OPTARG" ;;
	  c ) chr="$OPTARG" ;;
   esac
done

# Store current date and time in variable and start time of the script
now=$( date )
start=`date +%s`
echo "Start at ${start}" 
echo "Calculating polygenic scores using plink2 on chromosome $chr" 

base_dir=/data/clusterfs/lag/users/jitame/SENT_CORE

module add plink/1.9b6

plink \
--bfile $base_dir/geno/mostest/in/mostest_geno_in_c$chr \
--score $input_efs 2 4 6 sum --allow-no-sex --out $outp


# Store current date and time in variable and calculate the runtime
now=$( date )
checkpoint=`date +%s`
runtime=$(((checkpoint-start)/60))
printf "\n Elapsed time is "${runtime}" minutes.\n\n"
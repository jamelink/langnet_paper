#!/bin/sh
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash 
#$ -e /data/clusterfs/lag/users/jitame/logs/
#$ -o /data/clusterfs/lag/users/jitame/logs/
#$ -M Jitse.Amelink@mpi.nl
#$ -N cca_pgs
#$ -m beas


python_path="/home/jitame/bin/anaconda3/envs/results_env/bin/python"
file_name="/home/jitame/bin/code/AICHA/06_visualization/03_cca.py"
output_path="/data/clusterfs/lag/users/jitame/SENT_CORE/geno"
cd $output_path

$python_path $file_name


#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J test_speed.graphkir.par
#SBATCH -p ngs92G
#SBATCH -n 14
#SBATCH --mem=92G
#SBATCH -o test_speed.graphkir.par.std.log
#SBATCH -e test_speed.graphkir.par.err.log

cd /home/linnil1tw/kir/
module load pkg/Anaconda3 libs/singularity
conda activate .conda/
date
ls data_tmp/linnil1_syn_s2022.*1.fq | awk -F. '{print $1 "." $2 "." $3}' | parallel -j 14 kirpipe {} --tools graphkir --thread 1 --engine singularity_linnil1
date

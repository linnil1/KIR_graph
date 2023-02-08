#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J test_speed.graphkir
#SBATCH -p ngs92G
#SBATCH -n 14
#SBATCH --mem=92G
#SBATCH -o test_speed.graphkir.std.log
#SBATCH -e test_speed.graphkir.err.log

cd /home/linnil1tw/kir/
module load pkg/Anaconda3 libs/singularity
conda activate .conda/
date
kirpipe data_tmp/linnil1_syn_s2022.{}.30x_s1031 --tools graphkir --thread 14
date

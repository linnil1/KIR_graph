#!/usr/bin/bash
#SBATCH -A MST109178
#SBATCH -J test_speed.ping
#SBATCH -p ngs92G
#SBATCH -n 14
#SBATCH --mem=92G
#SBATCH -o test_speed.ping.std.log
#SBATCH -e test_speed.ping.err.log

cd /home/linnil1tw/kir/
module load pkg/Anaconda3 libs/singularity
conda activate .conda/
date
# Use other_kir instead of kirpipe
# is becuase ping require manual CN thresholding
# So use auto thresholding in research/other_kir.py
# make sure
#  comment out other tools
#  remove all data (otherwise it'll skip)
python research/other_kir.py
date

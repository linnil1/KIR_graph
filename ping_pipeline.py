# Generate
https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
wget https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier2016.06.05linux64.tgz
tar xf artbinmountrainier2016.06.05linux64.tgz
export PATH=$PWD/art_bin_MountRainier/:$PATH
Rscript -e 'install.packages(c("data.table", "stringr"))'
git clone https://github.com/wesleymarin/ping_paper_scripts.git
Rscript syntheticSequenceGen.R

git clone https://github.com/wesleymarin/PING.git
cd PING
docker build . -t ping

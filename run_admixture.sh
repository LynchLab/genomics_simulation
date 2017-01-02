#!/bin/bash
NUMBER=1000
TIME=$((NUMBER))
SITES=640000
#./population_simulation/pedigree_sim $NUMBER $TIME $SITES 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/states.txt.gz 2> ./sequences/var
cat ./sequences/states.txt.gz | gunzip - | python ./analysis_pipelines/correl_stream.py > ./analysis_files/admixture.csv
cd ./analysis_pipelines/
Rscript make_figure_3,4,5.rscript

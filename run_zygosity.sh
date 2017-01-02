#!/bin/bash
NUMBER=100
TIME=$((NUMBER))
SITES=80000
#./population_simulation/pedigree_sim $NUMBER $TIME $SITES 1 0 | gzip - > ./sequences/states.txt.gz 2> ./sequences/var
cat ./sequences/states.txt.gz | gunzip - | python-2.7.9 ./analysis_pipelines/correl.py > ./analysis_files/admixture.csv
cd ./analysis_pipelines/
Rscript make_figure_3,4,5.rscript

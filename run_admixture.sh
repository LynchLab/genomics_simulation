#!/bin/bash
NUMBER=100
TIME=$((8*NUMBER))
SITES=204800
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/states.txt.gz 2> ./sequences/var
cat ./sequences/states.txt.gz | gunzip - | python-2.7.9 ./analysis_pipelines/correl_stream.py > ./analysis_files/admixture.csv
cd ./analysis_pipelines/
Rscript make_figure_3,4,5.rscript
cat ../sequences/states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py > con_long32_6.txt; Rscript make_figure_new.rscript

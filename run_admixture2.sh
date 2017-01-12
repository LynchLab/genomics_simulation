#!/bin/bash
NUMBER=100
TIME=$((6*NUMBER))
SITES=204800

TIME1=$((4*NUMBER))
TIME2=$((2*NUMBER))

./population_simulation/pedigree_sim_ran $NUMBER $TIME $SITES 0.01 0.01 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/ran_states.txt.gz 2> ./sequences/var
./population_simulation/pedigree_sim_ran_co $NUMBER $TIME $SITES 0.01 0.01 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/ran_co_states.txt.gz 2> ./sequences/var
./population_simulation/pedigree_sim_sib $NUMBER $TIME $SITES 0.01 0.01 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/sib_states.txt.gz 2> ./sequences/var
./population_simulation/pedigree_sim_con $NUMBER $TIME $SITES 0.01 0.01 1 0 | cut -d '	' -f 1-100 | gzip - > ./sequences/con_states.txt.gz 2> ./sequences/var

./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 1 0 0 $TIME2 | cut -d '	' -f 1-100 | gzip - > ./sequences/mixi1_states.txt.gz 2> ./sequences/var
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 1 0 $TIME2 $TIME1 | cut -d '	' -f 1-100 | gzip - > ./sequences/mix12_states.txt.gz 2> ./sequences/var
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 1 0 $TIME1 $TIME | cut -d '	' -f 1-100 | gzip - > ./sequences/mix2p_states.txt.gz 2> ./sequences/var

cd ./analysis_pipelines/

cat ../sequences/ran_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/ran.txt.gz; 
cat ../sequences/ran_co_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/ran_co.txt.gz; 
cat ../sequences/sib_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/sib.txt.gz; 
cat ../sequences/con_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/con.txt.gz; 

cat ../sequences/mixi1_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/mixi1.txt.gz; 
cat ../sequences/mix12_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/mix12.txt.gz; 
cat ../sequences/mix2p_states.txt.gz | gunzip - | python-2.7.9 correl_stream2.py 1000 | gzip - > ../analysis_files/mix2p.txt.gz; 

Rscript make_figure_new.rscript

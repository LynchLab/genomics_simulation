#!/bin/bash
NUMBER=100
SITES=20480
DISTANCE=10000

TIME=$(echo "scale=4; 3.0*$NUMBER" | bc) 
TIME1=$(echo "scale=4; 0.4*$NUMBER" | bc) 
TIME2=$(echo "scale=4; 0.8*$NUMBER" | bc) 
TIME3=$(echo "scale=4; 1.2*$NUMBER" | bc) 
TIME4=$(echo "scale=4; 2.6*$NUMBER" | bc) 

echo $TIME1

./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 r 1 0 2> ./sequences/var | gzip - > ./sequences/ran_states.txt.gz 
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 r 1 0 2> ./sequences/var | gzip - > ./sequences/ran_co_states.txt.gz 
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 s 1 0 2> ./sequences/var | gzip - > ./sequences/sib_states.txt.gz 
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 0 2> ./sequences/var | gzip - > ./sequences/con_states.txt.gz 

./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 1 0 $TIME4 2> ./sequences/var | gzip - > ./sequences/mixi1_states.txt.gz
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 1 $TIME4 $TIME3 2> ./sequences/var | gzip - > ./sequences/mix12_states.txt.gz
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 1 $TIME3 $TIME2 2> ./sequences/var | gzip - > ./sequences/mix23_states.txt.gz 
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 1 $TIME2 $TIME1 2> ./sequences/var | gzip - > ./sequences/mix34_states.txt.gz 
./population_simulation/pedigree_sim $NUMBER $TIME $SITES 0.01 0.01 i 1 1 $TIME1 $TIME 2> ./sequences/var | gzip - > ./sequences/mix4p_states.txt.gz 

cd ./analysis_pipelines/

cat ../sequences/ran_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/ran.txt.gz; 
cat ../sequences/ran_co_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/ran_co.txt.gz; 
cat ../sequences/sib_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/sib.txt.gz; 
cat ../sequences/con_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/con.txt.gz; 

cat ../sequences/mixi1_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/mixi1.txt.gz; 
cat ../sequences/mix12_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/mix12.txt.gz; 
cat ../sequences/mix23_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/mix23.txt.gz; 
cat ../sequences/mix34_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/mix34.txt.gz; 
cat ../sequences/mix4p_states.txt.gz | gunzip - | ../fast_correl/fast_correl $DISTANCE $NUMBER 1 | gzip - > ../analysis_files/mix4p.txt.gz; 

Rscript make_figure_new.rscript

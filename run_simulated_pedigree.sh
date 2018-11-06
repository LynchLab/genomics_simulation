bash make_ms.sh
python simulate_on_pedigree.py ms > ms_out
#python ./ms_to_state.py ms > ms_state
#python analysis_pipelines/python_utilities/ms_to_cm.py ./ms_out 100 > cm.dat
#cat ms_state | ./gsl_correl/call_relatedness_lambda -r 50 -g cm.dat > gsl_out

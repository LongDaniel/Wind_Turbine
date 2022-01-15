cp fort.11.first fort.11
JOB_ID_1=$(qsub -V run_first_aegean.sh)
JOB_ID_2=$(qsub -V -W depend=afterok:$JOB_ID_1 run_second_aegean.sh)
echo 'Two tasks submited'
echo $JOB_ID_1
echo $JOB_ID_2

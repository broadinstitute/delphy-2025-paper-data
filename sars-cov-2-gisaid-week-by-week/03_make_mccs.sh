#!/usr/bin/bash

for d in outputs_by_{submission,collection}_date_{a,b}/to_epi_week*
do
    run_name=$(echo "$d" | grep -o 'to_epi_week_.*')
    #../treeannotator2 -heights ca -burnin 30 ${d}/${run_name}.trees ${d}/${run_name}.mcc
    ../delphy_mcc ${d}/${run_name}.trees ${d}/${run_name}.mcc
done

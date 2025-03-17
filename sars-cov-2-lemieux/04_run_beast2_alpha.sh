#!/bin/bash
time ../beast2 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run_alpha/ma_sars_cov_2_alpha.xml

../treeannotator2 -heights ca -burnin 30 beast2_run_alpha/output.trees beast2_run_alpha/ma_sars_cov_2_beast2_alpha.mcc

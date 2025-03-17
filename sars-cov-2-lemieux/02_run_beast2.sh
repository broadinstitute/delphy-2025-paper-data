#!/bin/bash
time ../beast2 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run/ma_sars_cov_2.xml

../treeannotator2 -heights ca -burnin 30 beast2_run/output.trees beast2_run/ma_sars_cov_2_beast2.mcc

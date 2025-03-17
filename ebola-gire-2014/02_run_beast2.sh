#!/bin/bash
time ../beast2 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run/ebola.xml

../treeannotator2 -heights ca -burnin 30 beast2_run/output.trees beast2_run/ebola_beast2.mcc

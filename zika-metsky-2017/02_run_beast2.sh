#!/bin/bash
time ../beast2 \
     -threads -1 \
     -beagle \
     -working \
     beast2_run/zika.xml

../treeannotator2 -heights ca -burnin 30 beast2_run/output.trees beast2_run/zika_beast2.mcc

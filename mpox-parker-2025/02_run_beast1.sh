#!/bin/bash
time ../beast1 \
     -threads 8 \
     -beagle \
     -working \
     beast_run/mpox-parker-2025-beast.xml

../treeannotator2 -heights ca -burnin 30 beast_run/Mpox_2poch_combined.trees beast_run/mpox-parker-2025-beast.mcc

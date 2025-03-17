#!/bin/bash
mkdir -p delphy_outputs_b
mkdir -p beast2_run

time ../delphy \
   --v0-in-fasta delphy_inputs/ebola.fasta \
   --v0-steps 1000000000 \
   --v0-out-log-file delphy_outputs_b/ebola_delphy.log \
   --v0-log-every 200000 \
   --v0-out-trees-file delphy_outputs_b/ebola_delphy.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_b/ebola_delphy.dphy \
   --v0-delphy-snapshot-every 1000000 \
   --v0-out-beast-xml beast2_run/ebola.xml

../treeannotator2 -heights ca -burnin 30 delphy_outputs_b/ebola_delphy.trees delphy_outputs_b/ebola_delphy.mcc

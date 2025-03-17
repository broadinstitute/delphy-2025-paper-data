#!/bin/bash
mkdir -p delphy_outputs_a

time ../delphy \
   --v0-in-fasta delphy_inputs/mpox-parker-2025.fasta \
   --v0-steps 1000000000 \
   --v0-out-log-file delphy_outputs_a/mpox-parker-2025.log \
   --v0-log-every 100000 \
   --v0-out-trees-file delphy_outputs_a/mpox-parker-2025.trees \
   --v0-tree-every 1000000 \
   --v0-out-delphy-file delphy_outputs_a/mpox-parker-2025.dphy \
   --v0-delphy-snapshot-every 1000000 \
   --v0-mpox-hack

../treeannotator2 -heights ca -burnin 30 delphy_outputs_a/mpox-parker-2025.trees delphy_outputs_a/mpox-parker-2025.mcc

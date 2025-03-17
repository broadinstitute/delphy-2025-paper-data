#!/bin/bash
mkdir -p delphy_outputs_b

time ../delphy \
   --v0-in-fasta delphy_inputs/mpox-otoole-2023.fasta \
   --v0-steps 200000000 \
   --v0-out-log-file delphy_outputs_b/mpox-otoole-2023.log \
   --v0-log-every 20000 \
   --v0-out-trees-file delphy_outputs_b/mpox-otoole-2023.trees \
   --v0-tree-every 200000 \
   --v0-out-delphy-file delphy_outputs_b/mpox-otoole-2023.dphy \
   --v0-delphy-snapshot-every 200000 \
   --v0-mpox-hack

../treeannotator2 -heights ca -burnin 30 delphy_outputs_b/mpox-otoole-2023.trees delphy_outputs_b/mpox-otoole-2023.mcc

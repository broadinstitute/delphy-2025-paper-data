#!/bin/bash
mkdir -p delphy_outputs_b

# 1889 samples => ~10,000,000,000 steps
time ../delphy \
   --v0-in-fasta delphy_inputs/h5n1-andersen-ebfdf65-ALL_full_dates_only.fasta \
   --v0-steps   10000000000 \
   --v0-out-log-file delphy_outputs_b/h5n1-andersen-ebfdf65-ALL_full_dates_only.log \
   --v0-log-every   1000000 \
   --v0-out-trees-file delphy_outputs_b/h5n1-andersen-ebfdf65-ALL_full_dates_only.trees \
   --v0-tree-every 10000000 \
   --v0-out-delphy-file delphy_outputs_b/h5n1-andersen-ebfdf65-ALL_full_dates_only.dphy \
   --v0-delphy-snapshot-every 10000000 \
   --v0-site-rate-heterogeneity

../treeannotator2 -heights ca -burnin 30 delphy_outputs_b/h5n1-andersen-ebfdf65-ALL_full_dates_only.trees delphy_outputs_b/h5n1-andersen-ebfdf65-ALL_full_dates_only.mcc

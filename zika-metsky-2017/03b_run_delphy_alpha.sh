#!/bin/bash
mkdir -p delphy_outputs_alpha_b
mkdir -p beast2_run_alpha

time ../delphy \
    --v0-in-fasta delphy_inputs/zika.fasta \
    --v0-steps 2000000000 \
    --v0-out-log-file delphy_outputs_alpha_b/zika_delphy_alpha.log \
    --v0-log-every 400000 \
    --v0-out-trees-file delphy_outputs_alpha_b/zika_delphy_alpha.trees \
    --v0-tree-every 2000000 \
    --v0-out-delphy-file delphy_outputs_alpha_b/zika_delphy_alpha.dphy \
    --v0-delphy-snapshot-every 2000000 \
    --v0-site-rate-heterogeneity \
    --v0-out-beast-xml beast2_run_alpha/zika_alpha.xml

../treeannotator2 -heights ca -burnin 30 delphy_outputs_alpha_b/zika_delphy_alpha.trees delphy_outputs_alpha_b/zika_delphy_alpha.mcc

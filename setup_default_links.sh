#!/usr/bin/bash

[ -f mafft ] || ln -s "`which mafft`"
[ -f delphy ] || ln -s "${HOME}/now/delphy/build/release/delphy"
[ -f delphy_mcc ] || ln -s "${HOME}/now/delphy/build/release/delphy_mcc"
[ -f treeannotator2 ] || ln -s "`which treeannotator2`"  # BEAST2's treeannotator
[ -f loganalyser2 ] || ln -s "`which loganalyser2`"  # BEAST2's loganalyser
[ -f sapling ] || ln -s "${HOME}/now/sapling/build/release/sapling"
[ -f iqtree2 ] || ln -s "${HOME}/tools/iqtree-2.3.6-Linux-intel/bin/iqtree2"
[ -f beast1 ] || ln -s "${HOME}/tools/BEASTv10.5.0/bin/beast" beast1
[ -f beast2 ] || ln -s "${HOME}/github/CompEvol/beast2.6.2/bin/beast" beast2

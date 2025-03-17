#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from pathlib import Path

def plot_tree(mcc_filename, out_pdf_filename):
    mcc_tree = bt.loadNexus(mcc_filename)
    mcc_tree.traverse_tree()
    mcc_tree.setAbsoluteTime(bt.decimalDate("2014-06-16"))  # Hard-coded, but ok

    fig,ax = plt.subplots(figsize=(5,7),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black'
    s_func=lambda k: 1
    
    mcc_tree.sortBranches(descending=True)
    for k in mcc_tree.getInternal():
        k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,x_attr=x_attr,colour='#5B5B5C',width=0.5)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)
    
    ax.set_xlim(2014.15, 2014.45);
    ax.set_ylim(-5, mcc_tree.ySpan+5);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2014-03-01"),
        bt.decimalDate("2014-04-01"),
        bt.decimalDate("2014-05-01"),
        bt.decimalDate("2014-06-01"),
       ])
    ax.set_xticklabels([
        "1 Mar\n2014",
        "1 Apr\n2014",
        "1 May\n2014",
        "1 Jun\n2014",
       ]);
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

plot_tree('./delphy_outputs_a/ebola_delphy.mcc', 'plots/DelphyAMcc.pdf')
plot_tree('./delphy_outputs_alpha_a/ebola_delphy_alpha.mcc', 'plots/DelphyAMccAlpha.pdf')
plot_tree('./delphy_outputs_b/ebola_delphy.mcc', 'plots/DelphyBMcc.pdf')
plot_tree('./delphy_outputs_alpha_b/ebola_delphy_alpha.mcc', 'plots/DelphyBMccAlpha.pdf')
plot_tree('./beast2_run/ebola_beast2.mcc', 'plots/Beast2Mcc.pdf')
plot_tree('./beast2_run_alpha/ebola_beast2_alpha.mcc', 'plots/Beast2MccAlpha.pdf')
plot_tree('./ml/tt/timetree.nexus', 'plots/ML.pdf')
plot_tree('./ml/tt_alpha/timetree.nexus', 'plots/MLAlpha.pdf')

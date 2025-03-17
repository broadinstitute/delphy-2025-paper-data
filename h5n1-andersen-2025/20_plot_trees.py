#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import csv
from scipy import stats
import pandas as pd
import math
import numpy as np

# Read in metadata
# ===========================
print("\nReading in metadata...")
sample_id_2_state = {}
with open('delphy_inputs/h5n1-andersen-ebfdf65_metadata.csv', 'r') as f:
    r = csv.DictReader(f)
    for record in r:
        sample_id_2_state[record['id']] = record['geo']

# State color map
# $ cut -d',' -f2 delphy_inputs/h5n1-andersen-ebfdf65_metadata.csv | sort | uniq -c | sort -nr
#     962 -          # Almost all (all?) of these are SRRs that have no GenBank counterpart
#     808 CA
#     271 CO
#     240 TX
#     142 ID
#     108 MI
#      94 MN
#      56 IA
#      45 WY
#      42 SD
#      40 NM
#      20 OH
#      11 KS
#       4 OK
#       4 NC
#       3 UT
#       1 geo

state_2_color = {
    'CA': 'blue',
    'CO': 'orange',
    'TX': 'red',
    'ID': 'green',    
    'MI': 'yellow',
    'MN': 'cyan',
    'IA': 'darkgrey',
    'WY': 'grey',
    'SD': 'magenta',
    'NM': 'teal',
    'OH': 'purple',
    'KS': 'lightgreen',
    'OK': 'pink',
    'NC': 'brown',
    'UT': 'black',
    '-':  'lightgrey',
}

def state_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    return sample_id_2_state.get(sample_id, "-")

def color_of(leaf_name):
    return state_2_color[state_of(leaf_name)]

def simple_parsimony(tree, tip_state_assigner):
    # Post-order traversal implementing Fitch algorithm
    index_2_fs = {}
    def go(node):
        if not node.is_leaf():
            for child in node.children:
                go(child)
        
        if node.is_leaf():
            fs = tip_state_assigner(node)
            if isinstance(fs, str):
                fs = frozenset([fs])  # Singleton set
            fs = frozenset(fs)  # Make lists, sets and frozensets into just frozensets
        else:
            child_fss = [index_2_fs[child.index] for child in node.children]

            # Fitch rule: intersection if nonempty, otherwise union
            fs = frozenset.intersection(*child_fss)
            if len(fs) == 0:
                fs = frozenset.union(*child_fss)

        index_2_fs[node.index] = fs
    go(tree.root)

    # Pre-order traversal to assign states (pick "lowest" state on ties according to natural sort order of keys)
    index_2_state = {}
    def gogo(node, parent_state):
        node_fs = index_2_fs[node.index]
        
        if parent_state in node_fs:
            node_state = parent_state
        else:
            node_state = sorted(list(node_fs))[0]
            if len(node_fs) > 1:
                print(f'WARNING: Arbitrarily breaking tie between {node_fs} in favor of {node_state}')
        index_2_state[node.index] = node_state

        if not node.is_leaf():
            for child in node.children:
                gogo(child, node_state)
    gogo(tree.root, None)

    return index_2_state


def extract_tMRCA_distr(log_filename, burnin, t0):
    # Extremely brittle way of extracting a tMRCA distribution from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    return [t0 - h for h in data['TreeHeight']]


def plot_tree(mcc_filename, out_pdf_filename, tMRCAs):
    mcc_tree = bt.loadNexus(mcc_filename)
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
    mcc_tree.setAbsoluteTime(bt.decimalDate(max_date))

    index_2_cat = simple_parsimony(mcc_tree, lambda node: state_of(node.name))
    
    fig,ax = plt.subplots(figsize=(6,6),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 15 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    mcc_tree.drawTree()

    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: state_2_color[index_2_cat[node.index]],
                      width=1)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=True)

    mean_tMRCA = np.mean(tMRCAs)  # Slightly different from root height because tMRCA was sampled much more often than trees
    kde = stats.gaussian_kde(tMRCAs)
    tt = np.linspace(bt.decimalDate('2023-10-01'), bt.decimalDate('2024-03-01'), 200)
    ax.plot(tt, -20+60*kde(tt), color='#888888', lw=0.5)
    ax.fill_between(tt, -20+60*kde(tt), -20, color='#000000', alpha=0.2)

    ax.set_xlim(bt.decimalDate("2023-09-01"), bt.decimalDate("2025-01-01"))
    
    ax.set_ylim(-20, mcc_tree.ySpan+20);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    # ax.get_yaxis().set_ticks([
    #     -5+10*0.00,
    #     -5+10*0.50,
    #     -5+10*1.00,
    #     -5+10*1.50,
    # ]);
    # ax.get_yaxis().set_ticklabels([
    #     '0',
    #     '50',
    #     '100',
    #     '150',
    # ]);
    # ax.set_ylabel('Density', loc='bottom')
    # ax.yaxis.set_label_coords(-0.12, 0.12)  # Trial and error

    ax.get_xaxis().set_ticks([
        bt.decimalDate("2023-10-01"),
        bt.decimalDate("2024-01-01"),
        bt.decimalDate("2024-04-01"),
        bt.decimalDate("2024-07-01"),
        bt.decimalDate("2024-10-01"),
        bt.decimalDate("2025-01-01"),
    ])
    ax.set_xticklabels([
        "1 Oct\n2023",
        "1 Jan\n2024",
        "1 Apr\n2024",
        "1 Jul\n2024",
        "1 Oct\n2024",
        "1 Jan\n2025",
    ]);
    ax.set_xlabel('Time (years)')

    plt.axvspan(bt.decimalDate("2023-10-01"), bt.decimalDate("2023-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2023-12-01"), bt.decimalDate("2024-01-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-02-01"), bt.decimalDate("2024-03-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-04-01"), bt.decimalDate("2024-05-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-06-01"), bt.decimalDate("2024-07-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-08-01"), bt.decimalDate("2024-09-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-10-01"), bt.decimalDate("2024-11-01"), color='black', alpha=0.05)
    plt.axvspan(bt.decimalDate("2024-12-01"), bt.decimalDate("2025-01-01"), color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    custom_lines = []
    custom_labels = []
    for cat, color in state_2_color.items():
        if cat == '-':
            continue    # Suppress extremely rare but confusing category in legend
        custom_lines.append(Line2D([0], [0], color=color, lw=3))
        custom_labels.append(cat)
    ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.2,0.2), framealpha=0.5)
    
    def annotate_node(node, text, diamond=False, **kwargs):
        x = kwargs.get('x', x_attr(node))
        y = node.y

        # Allow kwargs to override defaults
        allkwargs = dict(horizontalalignment='right', verticalalignment='bottom') | kwargs
        if 'x' in allkwargs:
            del allkwargs['x']
        plt.text(x-0.03, y, text, **allkwargs)
        if diamond:
            ax.plot(x, y, marker='D', color='black')

    #annotate_node(mrca_node_of(mcc_tree, ['ON563414']), 'B.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'ON563414']), 'A.1.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON676708', 'OP612681']), 'A.1')
    #annotate_node(mrca_node_of(mcc_tree, ['ON674051', 'ON676707']), 'A.2')

    annotate_node(mcc_tree.root, 'tMRCA', diamond=True, x=mean_tMRCA, rotation='vertical', verticalalignment='top')
    
    ax.plot([mean_tMRCA, mean_tMRCA], [mcc_tree.root.y, -20], color='black', linestyle='--', lw=0.75)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

for run in ('a', 'b'):
    mcc_filename = f'./delphy_outputs_{run}/h5n1-andersen-ebfdf65-ALL_full_dates_only.mcc'
    
    mcc_tree = bt.loadNexus(mcc_filename)
    leaves = mcc_tree.traverse_tree()
    max_date = max(leaf.name.split('|')[-1] for leaf in leaves)  # e.g., 2020-03-01
    
    tMRCAs = extract_tMRCA_distr(f'./delphy_outputs_{run}/h5n1-andersen-ebfdf65-ALL_full_dates_only.log', burnin=0.30, t0=bt.decimalDate(max_date))
    plot_tree(mcc_filename, f'plots/{run.upper()}DelphyMcc.pdf', tMRCAs)

#!/usr/bin/env python

import baltic as bt
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import csv

# Read in metadata
# ===========================
print("\nReading in metadata...")
sample_id_2_geo = {}
with open('delphy_inputs/zika_metadata.csv', 'r') as f:
    r = csv.DictReader(f)
    for record in r:
        sample_id_2_geo[record['id']] = record['Geo']

# Simplified geo labels
geo_2_simple_geo = {
    'Suriname':              'Other',  # unlabeled in paper
    'Brazil':                'BRA',
    'Dominican_Republic':    'DOM',
    'Haiti':                 'HTI',
    'USA':                   'USA',
    'Jamaica':               'JAM',
    'Honduras':              'HND',
    
    'El_Salvador/Guatemala': 'UNK',
    'Guatemala':             'UNK',

    'Puerto_Rico':           'PRI',
    'Colombia':              'COL',
    'Martinique':            'MTQ',
}

# Extracted from Metsky et al, 2017
simple_geo_colors = {
    'BRA': '#0f7232',
    'DOM': '#97a5d1',
    'HTI': '#e5007e',
    'USA': '#d83433',
    'JAM': '#009ee2',
    'HND': '#f29100',
    'UNK': '#ffde00',
    'UNK': '#ffde00',
    'PRI': '#63247d',
    'COL': '#b27e4b',
    'MTQ': '#a38dc2',
    'Other': '#000000',  # unlabeled in paper
}

def simple_geo_of(leaf_name):
    sample_id, *rest = leaf_name.split('|')
    geo = sample_id_2_geo.get(sample_id, "?")
    return geo_2_simple_geo.get(geo, 'Other')

def color_of(leaf_name):
    return simple_geo_colors[simple_geo_of(leaf_name)]

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

def plot_tree(mcc_filename, out_pdf_filename):
    mcc_tree = bt.loadNexus(mcc_filename)
    mcc_tree.traverse_tree()
    mcc_tree.setAbsoluteTime(bt.decimalDate("2016-10-10"))  # Hard-coded, but ok

    index_2_simple_geo = simple_parsimony(mcc_tree, lambda node: simple_geo_of(node.name))
    
    fig,ax = plt.subplots(figsize=(5,7),facecolor='w')

    x_attr=lambda k: k.absoluteTime
    c_func=lambda k: 'black' if k.branchType != 'leaf' else color_of(k.name)
    s_func=lambda k: 20 if (k.branchType == 'leaf') else 1
    
    mcc_tree.sortBranches(descending=True)
    for k in mcc_tree.getInternal():
        k.children.reverse()
    mcc_tree.drawTree()
    
    mcc_tree.plotTree(ax,
                      x_attr=x_attr,
                      colour=lambda node: simple_geo_colors[index_2_simple_geo[node.index]],
                      width=0.5)
    mcc_tree.plotPoints(
        ax,
        x_attr=x_attr,
        size=s_func,
        colour=c_func,
        zorder=100,
        outline=False)

    ax.set_xlim(2013.5, 2017.0)
    ax.set_ylim(-5, mcc_tree.ySpan+5);
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.get_yaxis().set_ticks([]);
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2013-01-01"),
        bt.decimalDate("2014-01-01"),
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2016-01-01"),
        bt.decimalDate("2017-01-01"),
       ])
    ax.set_xticklabels([
        "1 Jan\n2013",
        "1 Jan\n2014",
        "1 Jan\n2015",
        "1 Jan\n2016",
        "1 Jan\n2017",
       ]);
    
    plt.axvspan(2014.0, 2015.0, color='black', alpha=0.05)
    plt.axvspan(2016.0, 2017.0, color='black', alpha=0.05)

    # Legend code adapted from https://matplotlib.org/stable/gallery/text_labels_and_annotations/custom_legends.html
    custom_lines = []
    custom_labels = []
    for simple_geo, color in simple_geo_colors.items():
        custom_lines.append(Line2D([0], [0], color=color, lw=3))
        custom_labels.append(simple_geo)
    ax.legend(custom_lines, custom_labels, bbox_to_anchor=(0.15, 0.95))
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

plot_tree('./delphy_outputs_a/zika_delphy.mcc', 'plots/DelphyAMcc.pdf')
plot_tree('./delphy_outputs_b/zika_delphy.mcc', 'plots/DelphyBMcc.pdf')
plot_tree('./delphy_outputs_alpha_a/zika_delphy_alpha.mcc', 'plots/DelphyAMccAlpha.pdf')
plot_tree('./delphy_outputs_alpha_b/zika_delphy_alpha.mcc', 'plots/DelphyBMccAlpha.pdf')
plot_tree('./beast2_run/zika_beast2.mcc', 'plots/Beast2Mcc.pdf')
plot_tree('./beast2_run_alpha/zika_beast2_alpha.mcc', 'plots/Beast2MccAlpha.pdf')
plot_tree('./ml/tt/timetree.nexus', 'plots/ML.pdf')
plot_tree('./ml/tt_alpha/timetree.nexus', 'plots/MLAlpha.pdf')

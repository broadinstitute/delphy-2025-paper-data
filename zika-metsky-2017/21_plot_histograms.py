#!/usr/bin/env python

import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path

def extract_ml_results(iqtree_log_filename, tt_folder_name, look_for_alpha=False):
    """Extremely brittle IQ-Tree log file parser, enough for these plots though"""

    results = {}
    
    with open(iqtree_log_filename, 'r') as f:
        lines = f.readlines()

    # Extract HKY kappa
    # -----------------
    # Looking for a section like this:
    #
    # Rate parameter R:
    # 
    #   A-C: 1.0000
    #   A-G: 9.0197
    #   A-T: 1.0000
    #   C-G: 1.0000
    #   C-T: 9.0197
    #   G-T: 1.0000
    #
    rate_line_idx = [i
                     for i in range(len(lines))
                     if lines[i] == 'Rate parameter R:\n'][0]
    results['kappa'] = float(lines[rate_line_idx+3].strip().split(':')[1])

    # Extract HKY stationary state frequencies
    # ----------------------------------------
    # Looking for a section like this:
    #
    # State frequencies: (estimated with maximum likelihood)
    # 
    #   pi(A) = 0.3192
    #   pi(C) = 0.214
    #   pi(G) = 0.198
    #   pi(T) = 0.2688
    #
    state_freqs_line_idx = [i
                            for i in range(len(lines))
                            if lines[i].startswith('State frequencies')][0]
    results['freqParameter.1'] = float(lines[state_freqs_line_idx+2].split('=')[1])
    results['freqParameter.2'] = float(lines[state_freqs_line_idx+3].split('=')[1])
    results['freqParameter.3'] = float(lines[state_freqs_line_idx+4].split('=')[1])
    results['freqParameter.4'] = float(lines[state_freqs_line_idx+5].split('=')[1])

    # Extract site-rate heterogeneity parameter alpha
    # -----------------------------------------------
    # Looking for a section like this:
    #
    # Model of rate heterogeneity: Gamma with 4 categories
    # Gamma shape alpha: 998.4
    #
    if look_for_alpha:
        alpha_line_idx = [i
                          for i in range(len(lines))
                          if lines[i].startswith('Gamma shape alpha')][0]
        results['gammaShape'] = float(lines[alpha_line_idx].split(':')[1])

    # Extract mutation rate
    # ---------------------
    # Looking in tt/molecular_clock.txt for a line like this:
    #
    # --rate:	1.785e-03
    with (Path(tt_folder_name) / 'molecular_clock.txt').open('r') as f:
        for line in f:
            if 'rate' in line:
                results['clockRate'] = float(line.split(':')[1])

    # Extract tree height
    # -------------------
    # Read all node dates from tt/dates.tsv and extract rate
    with (Path(tt_folder_name) / 'dates.tsv').open('r') as f:
        dates = []
        for line in f:
            if line.startswith('#'):
                continue
            dates.append(float(line.split()[-1]))

        results['TreeHeight'] = max(dates) - min(dates)

    # Extract coalescent prior parameters
    # -----------------------------------
    # Assume tt/skyline.tsv has two useful lines, like this:
    #
    # #Skyline assuming 50.0 gen/year and approximate confidence bounds (+/- 2.000000 standard deviations of the LH)
    # #date 	N_e 	lower 	upper
    # 2014.175	6.230e+00	2.804e+00	1.384e+01
    # 2014.462	1.254e+01	9.726e+00	1.617e+01
    with (Path(tt_folder_name) / 'skyline.tsv').open('r') as f:
        lines = f.readlines()
        assert len(lines) >= 4
        min_date, min_Ne, *_ = map(float, lines[2].split())
        max_date, max_Ne, *_ = map(float, lines[3].split())

        rho = 1./50.  # Default TreeTime assumption = 50 generations per year; rho = generation time in years
        min_Ne_rho = min_Ne*rho
        max_Ne_rho = max_Ne*rho

        results['ePopSize'] = max_Ne_rho
        results['growthRate'] = math.log(max_Ne/min_Ne) / (max_date - min_date)
                
    return results

def compare_logs(delphy_log_filename, beast2_log_filename, out_pdf_filename, burnin=0.10, show_alpha=False, ml_results={}):
    beast2_raw_data = pd.read_table(beast2_log_filename, comment='#')
    num_pts = len(beast2_raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    beast2_data = beast2_raw_data[burnin_pts:]
    print(f'Read in BEAST2 log in {beast2_log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(beast2_data)} usable')

    delphy_raw_data = pd.read_table(delphy_log_filename, comment='#')
    num_pts = len(delphy_raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    delphy_data = delphy_raw_data[burnin_pts:]
    print(f'Read in Delphy log in {delphy_log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(delphy_data)} usable')

    fig,ax = plt.subplots(3, 4, figsize=(11,8.5), facecolor='w')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, hspace=0.3)

    def make_subplot(i, colname, title, adjuster=lambda x: x, colname_beast2=None):
        if colname_beast2 is None:
            colnames_beast2 = [x for x in beast2_data.columns if x.startswith(colname)]
            if len(colnames_beast2) == 0:
                raise ValueError(f'Column {colname} not found in BEAST2 log')
            if len(colnames_beast2) > 1:
                raise ValueError(f'Column {colname} ambiguous in BEAST2 log (matches: {colnames_beast2})')
            colname_beast2 = colnames_beast2[0]
        
        plt.subplot(3,4,i)
        plt.title(title)
        plt.xlabel(None)
        
        this_ax = sb.kdeplot(adjuster(delphy_data[colname]), fill=True, color='green');
        this_ax = sb.kdeplot(adjuster(beast2_data[colname_beast2]), fill=True, color='blue');
        this_ax.set(xlabel=None)
        this_ax.get_yaxis().set_visible(False)

        if colname in ml_results:
            (minx, maxx) = this_ax.get_xaxis().get_data_interval()
            (miny, maxy) = this_ax.get_yaxis().get_data_interval()
            value = adjuster(ml_results[colname])
            if 0.5*minx <= value <= 2*maxx:
                this_ax.axvline(value, color='red')
            else:
                ypos = 0.5*(miny+maxy)
                if value < minx:
                    this_ax.annotate("",
                                     xy=(minx, ypos),
                                     xytext=(minx+0.15*(maxx-minx), ypos),
                                     arrowprops=dict(arrowstyle="->", color='red'))
                else:
                    this_ax.annotate("",
                                     xy=(maxx, ypos),
                                     xytext=(maxx-0.15*(maxx-minx), ypos),
                                     arrowprops=dict(arrowstyle="->", color='red'))
            
        return this_ax

    make_subplot( 1, 'clockRate', r'Mutation rate (x$10^{-3}$/site/year)', lambda mu: mu*1000)
    make_subplot( 2, 'TreeHeight', r'Tree Height (years)')
    this_ax = make_subplot( 3, 'kappa', r'HKY $\kappa$ parameter')
    fig.legend(this_ax.get_legend_handles_labels()[0], ['BEAST2', 'Delphy'], loc='lower right')
    if show_alpha:
        make_subplot( 4, 'gammaShape', r'Heterogeneity Spread $\alpha$')
    else:
        fig.delaxes(ax[0][3])
    
    make_subplot( 5, 'freqParameter.1', r'HKY $\pi_A$', colname_beast2='freqParameter.s:input_alignment1')
    make_subplot( 6, 'freqParameter.2', r'HKY $\pi_C$', colname_beast2='freqParameter.s:input_alignment2')
    make_subplot( 7, 'freqParameter.3', r'HKY $\pi_G$', colname_beast2='freqParameter.s:input_alignment3')
    make_subplot( 8, 'freqParameter.4', r'HKY $\pi_T$', colname_beast2='freqParameter.s:input_alignment4')

    make_subplot( 9, 'CoalescentExponential', r'Coalescent Prior')
    make_subplot(10, 'ePopSize', r'Final Pop Size (years)')
    make_subplot(11, 'growthRate', r'Pop Growth Rate (1/year)')
    fig.delaxes(ax[2][3])
    
    fig.legend(['Delphy', 'BEAST2', 'TreeTime'], loc='lower right', 
           bbox_to_anchor=(0.95, 0.05))
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

Path('plots').mkdir(parents=True, exist_ok=True)

ml_results = extract_ml_results('ml/zika.iqtree', 'ml/tt', False)
compare_logs(
    './delphy_outputs_a/zika_delphy.log',
    'beast2_run/output.log',
    'plots/DelphyAVsBeast2.pdf',
    burnin=0.30,
    ml_results=ml_results
)
compare_logs(
    './delphy_outputs_b/zika_delphy.log',
    'beast2_run/output.log',
    'plots/DelphyBVsBeast2.pdf',
    burnin=0.30,
    ml_results=ml_results
)

ml_results_alpha = extract_ml_results('ml/zika_alpha.iqtree', 'ml/tt_alpha', True)
compare_logs(
    './delphy_outputs_alpha_a/zika_delphy_alpha.log',
    'beast2_run_alpha/output.log',
    'plots/DelphyAVsBeast2Alpha.pdf',
    burnin=0.30,
    show_alpha=True,
    ml_results=ml_results_alpha
)
compare_logs(
    './delphy_outputs_alpha_b/zika_delphy_alpha.log',
    'beast2_run_alpha/output.log',
    'plots/DelphyBVsBeast2Alpha.pdf',
    burnin=0.30,
    show_alpha=True,
    ml_results=ml_results_alpha
)

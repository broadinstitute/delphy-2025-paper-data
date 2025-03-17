#!/usr/bin/env python

import pandas as pd
import subprocess
from io import StringIO
import math
from pathlib import Path
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

def analyse_log(log_filename, burnin=0.10):
    result = subprocess.run(
        [
            '../loganalyser2',
            '-burnin', str(math.floor(burnin*100)),
            log_filename,
        ],
        stdout=subprocess.PIPE
    )
    raw_analysis = result.stdout.decode('utf-8')
    return pd.read_fwf(StringIO(raw_analysis), index_col=0)

def compare_logs(delphy_log_filename, wallclock_delphy_min,
                 beast2_log_filename, wallclock_beast2_min,
                 out_pdf_filename,
                 burnin=0.10, show_alpha=False):
    
    analysis_delphy = analyse_log(delphy_log_filename, burnin=burnin)
    print(analysis_delphy)
    ess_delphy = {
        'Posterior': analysis_delphy['ESS']['posterior_for_Delphy'],
        'Likelihood': analysis_delphy['ESS']['likelihood_really_logG'],
        'Prior': analysis_delphy['ESS']['prior'],
        'Tree Height': analysis_delphy['ESS']['TreeHeight'],
        'Mutation Rate': analysis_delphy['ESS']['clockRate'],
        r'HKY $\pi_A$': analysis_delphy['ESS']['freqParameter.1'],
        r'HKY $\pi_C$': analysis_delphy['ESS']['freqParameter.2'],
        r'HKY $\pi_G$': analysis_delphy['ESS']['freqParameter.3'],
        r'HKY $\pi_T$': analysis_delphy['ESS']['freqParameter.4'],
        r'HKY $\kappa$': analysis_delphy['ESS']['kappa'],
        'Coalescent\nExponential': analysis_delphy['ESS']['CoalescentExponential'],
        'Final Pop Size': analysis_delphy['ESS']['ePopSize'],
        'Pop Growth Rate': analysis_delphy['ESS']['growthRate'],
    }
    if show_alpha:
        ess_delphy[r'Site Rate $\alpha$'] = analysis_delphy['ESS']['gammaShape']

    analysis_beast2 = analyse_log(beast2_log_filename, burnin=burnin)
    print(analysis_beast2)
    ess_beast2 = {
        'Posterior': analysis_beast2['ESS']['posterior'],
        'Likelihood': analysis_beast2['ESS']['likelihood'],
        'Prior': analysis_beast2['ESS']['prior'],
        'Tree Height': analysis_beast2['ESS']['TreeHeight.t:input_alignment'],
        'Mutation Rate': analysis_beast2['ESS']['clockRate.c:input_alignment'],
        r'HKY $\pi_A$': analysis_beast2['ESS']['freqParameter.s:input_alignment1'],
        r'HKY $\pi_C$': analysis_beast2['ESS']['freqParameter.s:input_alignment2'],
        r'HKY $\pi_G$': analysis_beast2['ESS']['freqParameter.s:input_alignment3'],
        r'HKY $\pi_T$': analysis_beast2['ESS']['freqParameter.s:input_alignment4'],
        r'HKY $\kappa$': analysis_beast2['ESS']['kappa.s:input_alignment'],
        'Coalescent\nExponential': analysis_beast2['ESS']['CoalescentExponential.t:input_alignment'],
        'Final Pop Size': analysis_beast2['ESS']['ePopSize.t:input_alignment'],
        'Pop Growth Rate': analysis_beast2['ESS']['growthRate.t:input_alignment'],
    }
    if show_alpha:
        ess_beast2[r'Site Rate $\alpha$'] = analysis_beast2['ESS']['gammaShape.s:input_alignment']

    sorted_names = sorted(list(ess_beast2.keys()), key=lambda x: ess_beast2[x])

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    plt.subplots_adjust(bottom=0.3)
    
    df = pd.DataFrame(
        [
            [f'{i+1}. {n}',
             ess_beast2[n]/wallclock_beast2_min,
             ess_delphy[n]/wallclock_delphy_min]
            for (i, n) in enumerate(sorted_names)
        ],
        columns=['Observable', 'BEAST2', 'Delphy'])
    
    # view data 
    df.plot(x='Observable', 
            kind='bar', 
            stacked=False, 
            ax=ax)

    ax.set_xlabel(None)
    ax.semilogy();
    ax.set_ylabel('Effective Samples per minute');
    #ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter(2));
    plt.xticks(rotation=90);
    plt.legend(loc='upper right')
    ax.set_ylim([1e-3, 1e7])
    
    left, bottom, width, height = [0.21, 0.66, 0.25, 0.20]
    ax2 = fig.add_axes([left, bottom, width, height])
    df2 = pd.DataFrame(
        [
            [f'{i+1}',
             (ess_delphy[n]/wallclock_delphy_min)/(ess_beast2[n]/wallclock_beast2_min)]
            for (i, n) in enumerate(sorted_names)
        ],
        columns=['Observable', 'Speedup'])
    
    # view data 
    df2.plot(x='Observable', 
             kind='bar', 
             stacked=False, 
             ax=ax2,
             legend=False)
    ax2.set_xlabel(None)
    ax2.set_ylabel('Speedup')
    ax2.set_ylim([0.1, 50000]);
    ax2.semilogy();
    plt.xticks(rotation=90);
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')


Path('plots').mkdir(parents=True, exist_ok=True)

wallclock_delphy_min = 46 + 7.566/60.0
wallclock_beast2_min = 345 + 47.019/60.0

compare_logs(
    './delphy_outputs_a/zika_delphy.log', wallclock_delphy_min,
    'beast2_run/output.log', wallclock_beast2_min,
    'plots/ESSDelphyAVsBeast2.pdf',
    burnin=0.30
)

wallclock_delphy_min = 46 + 9.582/60.0
wallclock_beast2_min = 345 + 47.019/60.0

compare_logs(
    './delphy_outputs_b/zika_delphy.log', wallclock_delphy_min,
    'beast2_run/output.log', wallclock_beast2_min,
    'plots/ESSDelphyBVsBeast2.pdf',
    burnin=0.30
)


wallclock_delphy_alpha_min = 59 + 57.758/60.0
wallclock_beast2_alpha_min = 984 + 41.206/60.0

compare_logs(
    './delphy_outputs_alpha_a/zika_delphy_alpha.log', wallclock_delphy_alpha_min,
    'beast2_run_alpha/output.log', wallclock_beast2_alpha_min,
    'plots/ESSDelphyAVsBeast2Alpha.pdf',
    burnin=0.30
)

wallclock_delphy_alpha_min = 59 + 10.386/60.0
wallclock_beast2_alpha_min = 984 + 41.206/60.0

compare_logs(
    './delphy_outputs_alpha_b/zika_delphy_alpha.log', wallclock_delphy_alpha_min,
    'beast2_run_alpha/output.log', wallclock_beast2_alpha_min,
    'plots/ESSDelphyBVsBeast2Alpha.pdf',
    burnin=0.30
)

#!/usr/bin/env python

import baltic as bt
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import math
from pathlib import Path
from scipy import stats

t0 = bt.decimalDate("2023-05-31")  # Hard-coded, but ok

def extract_log_params(log_filename, burnin):
    # Extremely brittle way of extracting distributions from the log files
    raw_data = pd.read_table(log_filename, comment='#')
    num_pts = len(raw_data)
    burnin_pts = math.floor(burnin * num_pts)
    data = raw_data[burnin_pts:]
    print(f'Read in log file {log_filename}')
    print(f' - {num_pts} entries, of which {burnin_pts} burn-in and {len(data)} usable')

    result = {}
    for key in [
            'exponential.popSize',
            'exponential.growthRate',
            'apobec3.clock.rate',
            'non_apobec3.clock.rate',
    ]:
        result[key] = list(data[key])

    return result

def plot_growth_curves(delphy_logs, beast_logs, out_pdf_filename):
    fig,ax = plt.subplots(figsize=(5,4),facecolor='w')

    years = np.linspace(2014.0, 2023.0, 1000)

    for (logs, color, linestyle) in [(delphy_logs, 'green', '--'), (beast_logs, 'blue', 'dotted')]:
        pop_min = []
        pop_2p5 = []
        pop_50 = []
        pop_97p5 = []
        pop_max = []
        for year in years:
            pops = sorted([N0 * np.exp(g*(year-t0))
                           for (N0, g) in zip(logs['exponential.popSize'], logs['exponential.growthRate'])])
            pop_min.append(pops[0])
            pop_2p5.append(pops[int(0.025*len(pops))])
            pop_50.append(pops[int(0.50*len(pops))])
            pop_97p5.append(pops[int(0.975*len(pops))])
            pop_max.append(pops[-1])

        #ax.fill_between(years, pop_min, pop_max, color=color, alpha=0.2)
        ax.fill_between(years, pop_2p5, pop_97p5, color=color, alpha=0.2)
        ax.plot(years, pop_50, color=color, linestyle=linestyle)
    
    ax.set_xlim(2014.5, 2022.7)
    ax.set_ylim(5e-5, 3.2e3);
    ax.set_yscale('log')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.get_yaxis().set_ticks([
        1e-4,
        1e-3,
        1e-2,
        1e-1,
        1e0,
        1e1,
        1e2,
        1e3,
    ]);
    ax.get_yaxis().set_ticklabels([
        '10$^{-4}$',
        '10$^{-3}$',
        '10$^{-2}$',
        '10$^{-1}$',
        '10$^{0}$',
        '10$^{1}$',
        '10$^{2}$',
        '10$^{3}$',
    ]);
    ax.set_ylabel('Effective population size (years)')
    
    ax.get_xaxis().set_ticks([
        bt.decimalDate("2015-01-01"),
        bt.decimalDate("2016-01-01"),
        bt.decimalDate("2017-01-01"),
        bt.decimalDate("2018-01-01"),
        bt.decimalDate("2019-01-01"),
        bt.decimalDate("2020-01-01"),
        bt.decimalDate("2021-01-01"),
        bt.decimalDate("2022-01-01"),
       ])
    ax.set_xticklabels([
        "2015",
        "2016",
        "2017",
        "2018",
        "2019",
        "2020",
        "2021",
        "2022",
       ]);
    ax.set_xlabel('Time (years)')

    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')

def plot_clock_rates(delphy_logs, beast_logs, out_pdf_filename):
    fig,ax = plt.subplots(figsize=(6,6),facecolor='w')

    mean_delphy_apobec = np.sum(delphy_logs['apobec3.clock.rate']) / len(delphy_logs['apobec3.clock.rate'])
    mean_beast_apobec = np.sum(beast_logs['apobec3.clock.rate']) / len(beast_logs['apobec3.clock.rate'])
    mean_delphy_non_apobec = np.sum(delphy_logs['non_apobec3.clock.rate']) / len(delphy_logs['non_apobec3.clock.rate'])
    mean_beast_non_apobec = np.sum(beast_logs['non_apobec3.clock.rate']) / len(beast_logs['non_apobec3.clock.rate'])
    
    bplot = ax.boxplot([delphy_logs['apobec3.clock.rate'],
                        beast_logs['apobec3.clock.rate'],
                        delphy_logs['non_apobec3.clock.rate'],
                        beast_logs['non_apobec3.clock.rate']],
                       widths=0.7,
                       patch_artist=True,
                       tick_labels=['APOBEC3\n(Delphy)', 'APOBEC3\n(BEAST)',
                                    'non-APOBEC3\n(Delphy)', 'non-APOBEC3\n(BEAST)'])

    for patch, color in zip(bplot['boxes'], ['#926568', '#926568',
                                             '#f4e7c8', '#f4e7c8']):
        patch.set_facecolor(color)
    for patch, median, color in zip(bplot['boxes'], bplot['medians'], ['#ded2d2', '#ded2d2',
                                                                       '#dab965', '#dab965']):
        patch.set_edgecolor(color)
        median.set_color(color)
    for whisker, color in zip(bplot['whiskers'], ['#ded2d2', '#ded2d2', '#ded2d2', '#ded2d2',
                                                  '#dab965', '#dab965', '#dab965', '#dab965']):
        whisker.set_color(color)
    for cap, color in zip(bplot['caps'], ['#926568', '#926568', '#926568', '#926568',
                                          '#f4e7c8', '#f4e7c8', '#f4e7c8', '#f4e7c8']):
        cap.set_color(color)
    for flier, color in zip(bplot['fliers'], ['#926568', '#926568',
                                              '#f4e7c8', '#f4e7c8']):
        flier.set_color(color)
        flier.set_marker('*')
        flier.set_markersize(1.0)
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    ax.set_ylim(5e-6, 2e-4);
    #ax.set_ylim(0, 17.6e-5);
    ax.set_yscale('log')
    ax.get_yaxis().set_ticks([
        1e-4,
        1e-5,
       #    0e-5,
       #  2.5e-5,
       #  5.0e-5,
       #  7.5e-5,
       # 10.0e-5,
       # 12.5e-5,
       # 15.0e-5,
       # 17.5e-5,
    ]);
    ax.get_yaxis().set_ticklabels([
        '1e-04',
        '1e-05',
       # '0.000000',
       # '0.000025',
       # '0.000050',
       # '0.000075',
       # '0.000100',
       # '0.000125',
       # '0.000150',
       # '0.000175',
    ])
    ax.set_ylabel('Clock rate (mutations / site / year)')
    ax.set_xlabel('Clock')

    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.1)

    plt.text(1, mean_delphy_apobec+0.000040, f'{mean_delphy_apobec*1e4:.1f}' + ' x 10$^{{-4}}$', horizontalalignment='center')
    plt.text(2, mean_beast_apobec+0.000040, f'{mean_beast_apobec*1e4:.1f}' + ' x 10$^{{-4}}$', horizontalalignment='center')
    plt.text(3, mean_delphy_non_apobec+0.000004, f'{mean_delphy_non_apobec*1e6:.1f}' + ' x 10$^{{-6}}$', horizontalalignment='center')
    plt.text(4, mean_beast_non_apobec+0.000004, f'{mean_beast_non_apobec*1e6:.1f}' + ' x 10$^{{-6}}$', horizontalalignment='center')
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')
    
def plot_doubling_time_dist(delphy_logs, beast_logs, out_pdf_filename):
    fig,ax = plt.subplots(figsize=(5,2.5),facecolor='w')

    for (logs, color, linestyle, label, yoffset) in [(delphy_logs, 'green', '-', '(Delphy)', 0.2),
                                                    (beast_logs, 'blue', 'dotted', '(BEAST)', 0.0)]:
        # e^g = 2^(1/tdbl)  => tdbl = ln(2) / g
        doubling_times = np.log(2) / np.array(logs['exponential.growthRate'])
        times = np.linspace(0, 5, 200)
        kde = stats.gaussian_kde(doubling_times)
        ax.fill_between(times, kde(times), 0, color=color, alpha=0.2)
        ax.plot(times, kde(times), color=color, lw=1.5, linestyle=linestyle)
    
        mean_t_dbl = np.mean(doubling_times)
        plt.vlines(mean_t_dbl, 0, 2.0, linestyle=linestyle, color=color, lw=1)
        plt.text(1.3+mean_t_dbl, 1.5+yoffset, f'Mean {label}: {mean_t_dbl:.2f} years', horizontalalignment='center')

        # Hack to get doubling times into Tracer1 for easy 95% HPD
        with open(f'plots/doubling_times_{label[1:-1]}.log', 'w') as f:
            f.write('Step\ttdbl\n')
            for i, tdbl in enumerate(doubling_times):
                f.write(f'{i}\t{tdbl:.3f}\n')
    
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    
    ax.set_ylabel('Density')
    ax.set_ylim(0, 2.1);
    ax.get_yaxis().set_ticks([
        0.0,
        0.5,
        1.0,
        1.5,
        2.0,
    ]);
    ax.get_yaxis().set_ticklabels([
        '0.0',
        '0.5',
        '1.0',
        '1.5',
        '2.0',
    ])
    
    ax.set_xlabel('Time (years)')
    ax.set_xlim(0, 5)
    ax.get_xaxis().set_ticks([
        0,
        1,
        2,
        3,
        4,
        5,
    ]);
    ax.get_xaxis().set_ticklabels([
        '0',
        '1',
        '2',
        '3',
        '4',
        '5',
    ])

    plt.subplots_adjust(left=0.12, right=0.93, top=0.9, bottom=0.2)
    
    plt.savefig(out_pdf_filename)
    print(f'Plot saved to {out_pdf_filename}')
    
Path('plots').mkdir(parents=True, exist_ok=True)

for run in ('a', 'b'):
    logs = extract_log_params(f'delphy_outputs_{run}/mpox-parker-2025.log', burnin=0.3)
    beast_logs = extract_log_params(f'beast_run/Mpox_2poch_combined.log', burnin=0.3)
    #plot_growth_curves(logs, beast_logs, f'plots/{run.upper()}GrowthCurve.pdf')
    #plot_clock_rates(logs, beast_logs, f'plots/{run.upper()}ClockRates.pdf')
    plot_doubling_time_dist(logs, beast_logs, f'plots/{run.upper()}DoublingTimeDistr.pdf')

from fact.io import read_h5py
from fact.analysis.statistics import calc_proton_obstime
import numpy as np
import click
import matplotlib.pyplot as plt
import yaml
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict
import astropy.units as u


def plot_hists(
        dfs,
        key,
        n_bins=100,
        limits=None,
        transform=lambda x: x,
        xlabel=None,
        yscale='linear',
        ax=None,
        ):

    if ax is None:
        ax = plt.gca()

    if limits is None:
        limits = [
            min(transform(df[key]).quantile(0.01) for df in dfs.values()),
            max(transform(df[key]).quantile(0.99) for df in dfs.values()),
        ]

    for label, df in dfs.items():
        ax.hist(
            transform(df[key].values),
            bins=n_bins,
            range=limits,
            weights=df['weight'],
            label=label.replace('_',r'\_'),
            histtype='step',
        )

    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel or key.replace('_',r'\_'))
    ax.legend()


@click.command()
@click.argument('data_path')
@click.argument('proton_path')
@click.argument('output')
@click.option('-c', '--config')
@click.option('-f', '--sample-fraction', default=0.3)
def main(data_path, proton_path, output, config, sample_fraction):

    if config is not None:
        with open(config) as f:
            config = yaml.safe_load(f)
    else:
        config = {}
    print(config)

    print('start reading')
    crab = read_h5py(data_path, key='events')
    crab_runs = read_h5py(data_path, key='runs')
    protons = read_h5py(proton_path, key='events')
    print('done')

    proton_obstime = calc_proton_obstime(n_events=7998*1500*20,
                                        spectral_index=-2.7,
                                        max_impact=400*u.m,
                                        viewcone=5*u.deg,
                                        e_min=100*u.GeV,
                                        e_max=30*u.TeV)
    # print('Proton obstime: ', proton_obstime)

    protons['weight'] = crab_runs.ontime.sum() / proton_obstime / sample_fraction
    # print('Crab ontime: ', crab_runs.ontime.sum()/3600, 'h')
    # print('Proton_weight_mean: ', protons.weight.mean())
    crab['weight'] = 1

    # print(protons.weight.sum() / len(crab))

    dfs = OrderedDict([
        ('Measured Crab Data', crab),
        ('Simulated Proton Data', protons),
    ])

    keys = set(crab.columns).intersection(set(protons.columns))
    wech = ['run', 'event', 'pointing_position_az', 'pointing_position_zd', 'theta_deg_off_1', 'theta_deg_off_2', 'theta_deg_off_3', 'theta_deg_off_4', 'theta_deg_off_5', 'weight']
    for key in wech:
        keys.remove(key)

    fig = plt.figure()
    # ax_ratio = plt.subplot2grid((4, 1), (3, 0))
    # ax_hist = plt.subplot2grid((4, 1), (0, 0), rowspan=3, sharex=ax_hist)
    ax_hist = fig.add_subplot(1, 1, 1)

    with PdfPages(output) as pdf:
        for key in sorted(list(keys), key=str.lower):
            print(key)

            kwargs = config.get(key, {})
            if 'transform' in kwargs:
                kwargs['transform'] = eval(kwargs['transform'])

            ax_hist.cla()
            plot_hists(dfs, key, ax=ax_hist, **kwargs)

            plt.tight_layout()
            pdf.savefig(fig)


if __name__ == '__main__':
    main()

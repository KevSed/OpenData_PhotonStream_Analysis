import photon_stream as ps
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from fact.io import read_h5py
from collections import OrderedDict
import click
import yaml
import numpy as np

'''
Erstellt data monte carlo Vergleichsplots
'''

def plot_histograms(dfs, key, n_bins=100, limits=None, transform=lambda x: x,
        xlabel=None, yscale='linear', ax=None):

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
            label=label,
            histtype='step',
            density=True
        )

    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel or key)
    ax.legend()


@click.command()
@click.argument('path')
@click.option('--cuts', is_flag=True)
def main(path, cuts):
    config = 'configs/data_mc.yaml'

    if cuts:
        precuts = '_precuts'
    else:
        precuts = ''

    print('Read in data...')
    gammas = read_h5py(path + 'gamma_sample{}.hdf5'.format(precuts), key='events')
    protons = read_h5py(path + 'proton_sample{}.hdf5'.format(precuts), key='events')
    data = read_h5py(path + 'crab_data{}.hdf5'.format(precuts), key='events')
    print('Done...')

    if config is not None:
        with open(config) as f:
            config = yaml.safe_load(f)
    else:
        config = {}
    print(config)

    fig = plt.figure()
    ax_hist = fig.add_subplot(1, 1, 1)

    keys = set(data.columns).intersection(set(protons.columns))
    wech = ['run', 'event', 'pointing_position_az', 'pointing_position_zd', 'theta_deg_off_1', 'theta_deg_off_2', 'theta_deg_off_3', 'theta_deg_off_4', 'theta_deg_off_5']
    for key in wech:
        if key in keys and key in set(gammas.columns):
            keys.remove(key)
    dfs = OrderedDict([
        ('Measured Crab Data', data),
        ('Simulated Proton Data', protons),
        ('Simulated Gamma Data', gammas),
    ])

    print('Plotting keys...')
    with PdfPages('{}gpd_mc_comp_cuts.pdf'.format(path)) as pdf:
        for key in sorted(list(keys)):
            print(key)

            kwargs = config.get(key, {})
            if 'transform' in kwargs:
                kwargs['transform'] = eval(kwargs['transform'])

            ax_hist.cla()
            plot_histograms(dfs, key, ax=ax_hist, **kwargs)

            pdf.savefig(fig)


if __name__ == '__main__':
    main()

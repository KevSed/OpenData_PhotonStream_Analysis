import numpy as np
import photon_stream as ps
from scipy import stats
import matplotlib.pyplot as plt
from fact.io import read_data
import click
import yaml
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from fact.io import read_h5py
from IPython import embed




def plot_histograms2d(dfs, key, h, threshold=0, n_bins=100,
        xlabel=None, yscale='linear', ax=None):

    if ax is None:
        ax = plt.gca()


    for label, df in dfs.items():
        if label == 'data':
            d = {}
            d[key] = df[key][df['gamma_prediction'] > threshold].values
        if label == 'proton':
            p = {}
            p[key] = df[key][df['gamma_prediction'] > threshold].values

    ax.hist2d(
       h['data_' + key][0],
       h['proton_' + key][0],
       bins=[h['data_' + key][1], h['proton_' + key][1]],
       range=[limits, limits],
       # label=label,
       normed=True
    )

    ax.set_xscale(yscale)
    ax.set_yscale(yscale)
    if xlabel is not None:
        ax.set_xlabel(xlabel + ' data')
        ax.set_ylabel(xlabel + ' Proton MC')
    else:
        ax.set_xlabel(key + ' data')
        ax.set_ylabel(key + ' Proton MC')

    # ks = stats.entropy(h['data_' + key][0], h['proton_' + key][0])
    ax.set_title('Proton MC vs. Data {}'.format(key))


def plot_histograms(dfs, key, h, threshold, n_bins=100, limits=None, transform=lambda x: x,
        xlabel=None, yscale='linear', ax=None):

    if ax is None:
        ax = plt.gca()

    if limits is None:
        limits = [
                min(transform(df[key]).quantile(0.01) for df in dfs.values()),
                max(transform(df[key]).quantile(0.99) for df in dfs.values()),
        ]

    for label, df in dfs.items():
#        d = {}
#        if label == 'data':
#            d[key] = df[key][df['gamma_prediction'] > threshold].values
#        if label == 'proton':
#            d[key] = df[key][df['gamma_prediction'] > threshold].values
        h[label + '_' + key] = ax.hist(
            transform(df[key][df['gamma_prediction'] > threshold].values),
            bins=n_bins,
            range=limits,
            label=label,
            histtype='step',
            density=True
        )

    ax.set_yscale(yscale)
    ax.set_xlabel(xlabel or key.replace('_', '-'))
    ax.legend()
    plt.tight_layout()
   # bin_width = (limits[1]-limits[0])/100
   # heights_d = h['data_' + key][0]*bin_width*len(df[key])
   # d = []
   # heights_p = h['proton_' + key][0]*bin_width*len(df[key])
   # p = []
   # for i in range(len(heights_d)):
   #     for j in range(int(heights_d[i])):
   #         d.append(limits[0] + i*bin_width)
   # for i in range(len(heights_p)):
   #     for j in range(int(heights_p[i])):
   #         p.append(limits[0] + i*bin_width)
   # ks = stats.ks_2samp(d[key], p[key])
   # ax.set_title('Proton-Data KS-Test p-value:{:.3e}'.format(ks))
   # ks = stats.entropy(h['data_' + key][0], h['proton_' + key][0])
   # ax.set_title('Proton-Data KL entropy:{:.3}'.format(ks))


@click.command()
@click.argument('path')
@click.option('--cuts', is_flag=True)
@click.option('--threshold', default=0.0, type=float)
def main(path, cuts, threshold):
    config = 'configs/data_mc.yaml'

    if cuts:
        precuts = '_precuts'
    else:
        precuts = ''

    print('Read in data...')
    gamma = read_h5py(path + '/gamma{}.hdf5'.format(precuts), key='events')
    proton = read_h5py(path + '/proton{}.hdf5'.format(precuts), key='events')
    data = read_h5py(path + '/crab_data{}.hdf5'.format(precuts), key='events')
    print('Done...')


    if config is not None:
        with open(config) as f:
            config = yaml.safe_load(f)
    else:
        config = {}
    print(config)

    fig = plt.figure()
    ax_hist = fig.add_subplot(1, 1, 1)

    keys = set(data.columns).intersection(set(proton.columns))
    wech = ['run', 'event', 'pointing_position_az', 'pointing_position_zd', 'theta_deg_off_1', 'theta_deg_off_2', 'theta_deg_off_3', 'theta_deg_off_4', 'theta_deg_off_5']
    for key in wech:
        if key in keys and key in set(gamma.columns):
            keys.remove(key)
    dfs = OrderedDict([
        ('data', data),
        ('proton', proton),
        ('gamma', gamma),
    ])

    # Dict for histograms
    h = {}
    print('Plotting keys...')
    with PdfPages('{}/pdf/feature_comp_cuts.pdf'.format(path)) as pdf:
        for key in sorted(list(keys)):
            print(key)

            kwargs = config.get(key, {})
            if 'transform' in kwargs:
                kwargs['transform'] = eval(kwargs['transform'])

            ax_hist.cla()
            plot_histograms(dfs, key, h, threshold, ax=ax_hist, **kwargs)

            pdf.savefig(fig)

#     print('Plotting keys in 2D...')
#     with PdfPages('{}/pdf/data_mc_2d.pdf'.format(path)) as pdf:
#         for key in sorted(list(keys)):
#             print(key)
# 
#             kwargs = config.get(key, {})
#             if 'transform' in kwargs:
#                 del kwargs['transform']
# 
#             ax_hist.cla()
#             plot_histograms2d(dfs, key, h, ax=ax_hist, **kwargs)
# 
#             pdf.savefig(fig)




if __name__ == '__main__':
    main()

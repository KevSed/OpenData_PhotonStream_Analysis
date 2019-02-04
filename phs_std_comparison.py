from matplotlib.colors import LogNorm, Normalize
import numpy as np
import matplotlib.pyplot as plt
from fact.io import read_data
from matplotlib.backends.backend_pdf import PdfPages
import click
from IPython import embed


def plot_comparison(df, feature, log=False, logz=False):
    plt.figure()
    df = df.dropna()

    if log:
        t = np.log10
    else:
        t = lambda x: x

    x = t(df[feature + '_phs'].dropna())
    y = t(df[feature + '_std'].dropna())
    plt.hist2d(
        x, y,
        bins=100,
        norm=LogNorm() if logz else Normalize(),
        range=[np.nanpercentile(x, [0.05, 99]), np.nanpercentile(y, [0.05, 99])],
    )
    plt.colorbar()
    plt.xlabel('PhotonStream')
    plt.ylabel('FACT-Tools standard analysis')
    plt.title('log10(' + feature + ')' if log else feature)
    plt.tight_layout()


def plot_comparison_hist(df, feature, log=False, logz=False):
    plt.figure()
    df = df.dropna()

    if log:
        t = np.log10
    else:
        t = lambda x: x

    x = t(df[feature + '_phs'].dropna())
    y = t(df[feature + '_std'].dropna())
    plt.hist(
        x,
        bins=100,
        range=np.nanpercentile(x, [0.05, 99]),
        label='PhotonStream',
        histtype='step',
    )
    plt.hist(
        y,
        bins=100,
        range=np.nanpercentile(x, [0.05, 99]),
        label='FACT-Tools Standard',
        histtype='step',
    )
    plt.xlabel('log10(' + feature + ')' if log else feature)
    plt.title('log10(' + feature + ')' if log else feature)
    plt.legend(loc='best')
    plt.tight_layout()


def plot_comparison_hist_all(std, phs, feature, log=False, logz=False):
    plt.figure()
    std = std.dropna()
    phs = phs.dropna()

    if log:
        t = np.log10
    else:
        t = lambda x: x

    x = t(phs[feature].dropna())
    y = t(std[feature].dropna())
    plt.hist(
        x,
        bins=100,
        range=np.nanpercentile(x, [0.05, 99]),
        label='PhotonStream',
        histtype='step',
    )
    plt.hist(
        y,
        bins=100,
        range=np.nanpercentile(x, [0.05, 99]),
        label='FACT-Tools Standard',
        histtype='step',
    )
    plt.xlabel('log10(' + feature + ')' if log else feature)
    plt.title('log10(' + feature + ')' if log else feature)
    plt.legend(loc='best')
    plt.tight_layout()


@click.command()
@click.argument('method', required=True)
@click.argument('type', required=True)
def main(method, type):
    std = read_data(
        '/home/ksedlaczek/Packages/open_crab_sample_analysis/dl2/{}.hdf5'.format(type),
        key='events',
        columns=[
            'run_id',
           #  'corsika_event_header_event_number',
           #  'corsika_event_header_num_reuse',
            'width',
            'length',
            'delta',
            'size',
            'cog_x',
            'cog_y',
           # 'pointing_position_az',
            'skewness_long',
            'skewness_trans',
            'night',
            'event_num',
        ]
    )

    std.rename(columns={
        'run_id': 'run',
        'event_num': 'event',
    }, inplace=True)
    std.set_index(['night', 'run', 'event'], inplace=True)

    phs = read_data('/net/big-tank/POOL/projects/fact/photon-stream/features/{}/{}_data.hdf5'.format(method, type), key='events')

    df = phs.join(
        std,
        how='inner',
        rsuffix='_std',
        lsuffix='_phs',
        on=('night', 'run', 'event'),
    )
    df.sort_index(axis=1, inplace=True)




    with PdfPages('/home/ksedlaczek/OUT_{}/pdf/std_phs_comparison_{}_{}.pdf'.format(method, method, type)) as pdf:
        print("Plotting width")
        plot_comparison(df, 'width')
        pdf.savefig()
        plt.close()

        print("Plotting length")
        plot_comparison(df, 'length')
        pdf.savefig()
        plt.close()

#        print("Plotting delta")
#        plot_comparison(df, 'delta')
#        pdf.savefig()
#        plt.close()

#        print("Plotting skewness_long")
#        plot_comparison(df, 'skewness_long')
#        pdf.savefig()
#        plt.close()

 #       print("Plotting skewness_trans")
 #       plot_comparison(df, 'skewness_trans', logz=False)
 #       pdf.savefig()
 #       plt.close()

        print("Plotting size")
        plot_comparison(df, 'size', log=True)
        pdf.savefig()
        plt.close()

#        print("Plotting cog x")
#        plot_comparison(df, 'cog_x', logz=False)
#        pdf.savefig()
#        plt.close()
#
#        print("Plotting cog y")
#        plot_comparison(df, 'cog_y', logz=False)
#        pdf.savefig()
#        plt.close()

#         print("Plotting 2D")
#         plot_comparison2(df, log=True)
#         pdf.savefig()
#         plt.close()
    with PdfPages('/home/ksedlaczek/OUT_{}/pdf/std_phs_comparison_hist_same_{}_{}.pdf'.format(method, method, type)) as pdf:
        print("Plotting width")
        plot_comparison_hist(df, 'width')
        pdf.savefig()
        plt.close()

        print("Plotting length")
        plot_comparison_hist(df, 'length')
        pdf.savefig()
        plt.close()
#
#        print("Plotting delta")
#        plot_comparison_hist(df, 'delta')
#        pdf.savefig()
#        plt.close()

#         print("Plotting skewness_long")
#         plot_comparison_hist(df, 'skewness_long')
#         pdf.savefig()
#         plt.close()

 #       print("Plotting skewness_trans")
 #       plot_comparison(df, 'skewness_trans', logz=False)
 #       pdf.savefig()
 #       plt.close()

        print("Plotting size")
        plot_comparison_hist(df, 'size', log=True)
        pdf.savefig()
        plt.close()

#        print("Plotting cog x")
#        plot_comparison_hist(df, 'cog_x', logz=False)
#        pdf.savefig()
#        plt.close()
#
#        print("Plotting cog y")
#        plot_comparison_hist(df, 'cog_y', logz=False)
#        pdf.savefig()
#        plt.close()

    with PdfPages('/home/ksedlaczek/OUT_{}/pdf/std_phs_comparison_hist_all_{}_{}.pdf'.format(method, method, type)) as pdf:
        print("Plotting width")
        plot_comparison_hist_all(std, phs, 'width')
        pdf.savefig()
        plt.close()

        print("Plotting length")
        plot_comparison_hist_all(std, phs, 'length')
        pdf.savefig()
        plt.close()
#
#        print("Plotting delta")
#        plot_comparison_hist_all(std, phs, 'delta')
#        pdf.savefig()
#        plt.close()
#
#        print("Plotting skewness_long")
#        plot_comparison_hist_all(std, phs, 'skewness_long')
#        pdf.savefig()
#        plt.close()

 #       print("Plotting skewness_trans")
 #       plot_comparison(df, 'skewness_trans', logz=False)
 #       pdf.savefig()
 #       plt.close()

        print("Plotting size")
        plot_comparison_hist_all(std, phs, 'size', log=True)
        pdf.savefig()
        plt.close()

#         print("Plotting cog x")
#         plot_comparison_hist_all(std, phs, 'cog_x', logz=False)
#         pdf.savefig()
#         plt.close()
#
#         print("Plotting cog y")
#         plot_comparison_hist_all(std, phs, 'cog_y', logz=False)
#         pdf.savefig()
#         plt.close()

if __name__ == '__main__':
    main()

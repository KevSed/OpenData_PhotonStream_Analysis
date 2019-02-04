from matplotlib.colors import LogNorm, Normalize
import numpy as np
import matplotlib.pyplot as plt
from fact.io import read_data
from matplotlib.backends.backend_pdf import PdfPages
import click
from IPython import embed


def plot_comparison(df, log=False, logz=False):
    plt.figure()
    df = df.dropna()

    if log:
        t = np.log10
    else:
        t = lambda x: x

    x = t(df['size'].dropna())
    y = t(df['gamma_energy_prediction'].dropna())
    plt.hist2d(
        x, y,
        bins=50,
        norm=LogNorm() if logz else None, #Normalize(),
        range=[np.nanpercentile(x, [1, 90]), np.nanpercentile(y, [1, 89])],
    )
    plt.colorbar()
    plt.xlabel('log10(size)')
    plt.ylabel('log10(Predicted Energy)')
    plt.title('Energy vs. size (data)')
    plt.tight_layout()



@click.command()
@click.argument('method', required=True)
def main(method):

    phs = read_data('/home/ksedlaczek/OUT_DBSCAN/crab_data_precuts.hdf5', key='events')

    df = phs
   #  df = phs.join(
   #      std,
   #      how='inner',
   #      rsuffix='_std',
   #      lsuffix='_phs',
   #      on=('run', 'event', 'reuse'),
   #  )
    df.sort_index(axis=1, inplace=True)

    with PdfPages('Energy_vs_size_{}.pdf'.format(method)) as pdf:
        print("Plotting 2D")
        plot_comparison(df, log=True)
        pdf.savefig()
        plt.close()


if __name__ == '__main__':
    main()

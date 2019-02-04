# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
from fact.coordinates import horizontal_to_camera
from fact.io import read_data
from astropy.coordinates import SkyCoord, AltAz
from fact.instrument.constants import LOCATION
from fact.coordinates.utils import to_astropy_time
import pandas as pd
from argparse import ArgumentParser
from IPython import embed

parser = ArgumentParser()
parser.add_argument('inputfile')
parser.add_argument('-o', '--outputfile')
parser.add_argument('-t', '--threshold', type=float)

args = parser.parse_args()

facttools = read_data('/home/ksedlaczek/Packages/open_crab_sample_analysis/build/crab_precuts.hdf5', key='events')

facttools = facttools.query(f'gamma_prediction >= 0.8').copy()

phs = read_data(args.inputfile, key='events')

facttools.set_index(['run_id', 'event_num', 'night'], inplace=True)
df = phs.join(
    facttools,
    how='inner',
    rsuffix='_std',
    lsuffix='_phs',
    on=('run', 'event', 'night'),
)
df.sort_index(axis=1, inplace=True)

if  args.threshold:
    phs = phs.query(f'gamma_prediction >= {args.threshold}').copy()

if 'source_position_zd' not in df.columns:
    frame = AltAz(location=LOCATION, obstime=to_astropy_time(pd.to_datetime(df['timestamp'])))

    crab = SkyCoord.from_name('Crab Nebula')
    crab_altaz = crab.transform_to(frame)

    df['source_position_zd_phs'] = crab_altaz.zen.deg
    df['source_position_az_phs'] = crab_altaz.az.deg

df['source_x'], df['source_y'] = horizontal_to_camera(
    df['source_position_zd'],
    df['source_position_az'],
    df['pointing_position_zd_phs'],
    df['pointing_position_az_phs'],
)

df['true_delta'] = np.arctan2(df['cog_y_phs'] - df['source_y'], df['cog_x_phs'] - df['source_x'])

df['delta_delta'] = df['true_delta'] - df['delta_phs']
df['delta_delta'][df['delta_delta'] < -np.pi/2] += 2 * np.pi

mask = df['disp_prediction_phs'] < 0
plt.hist(df.loc[mask, 'delta_delta'], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.3)
plt.hist(df.loc[~mask, 'delta_delta'], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.3)
plt.xlabel(r'$\mathrm{\Delta}\delta$ / rad')
plt.title(r'$\mathrm{\Delta}\delta$ between calculated $\delta$ and $\delta_{\mathrm{true}}$')
plt.tight_layout()
if args.outputfile:
    plt.savefig(args.outputfile + '.pdf', dpi=300)
else:
    plt.show()

plt.clf()

plt.hist2d(df['delta_phs'], df['delta_std'], bins=100, cmap='viridis')
plt.xlabel(r'$\delta$ / rad (PhotonStream)')
plt.ylabel(r'$\delta$ / rad (FACTtools)')
plt.title(r'$\delta_{\mathrm{PhotonStream}}$ vs $\delta_{\mathrm{FACTtools}}$')
plt.tight_layout()
plt.colorbar()
plt.savefig(args.outputfile + '_2d.pdf', dpi=300)

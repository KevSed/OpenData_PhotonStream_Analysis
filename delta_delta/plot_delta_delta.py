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

parser = ArgumentParser()
parser.add_argument('inputfile')
parser.add_argument('-o', '--outputfile')
parser.add_argument('-t', '--threshold', type=float)

args = parser.parse_args()

print("Reading in data...")
df = read_data(args.inputfile, key='events')
print("Done!")

if  args.threshold:
    df = df.query(f'gamma_prediction >= {args.threshold}').copy()

if 'source_position_zd' not in df.columns:
    frame = AltAz(location=LOCATION, obstime=to_astropy_time(pd.to_datetime(df['timestamp'])))

    crab = SkyCoord.from_name('Crab')
    crab_altaz = crab.transform_to(frame)

    df['source_position_zd'] = crab_altaz.zen.deg
    df['source_position_az'] = crab_altaz.az.deg

df['source_x'], df['source_y'] = horizontal_to_camera(
    df['source_position_zd'],
    df['source_position_az'],
    df['pointing_position_zd'],
    df['pointing_position_az'],
)

df['true_delta'] = np.arctan2(df['cog_y'] - df['source_y'], df['cog_x'] - df['source_x'])

df['delta_delta'] = df['true_delta'] - df['delta']
df['delta_delta'][df['delta_delta'] < -np.pi/2] += 2 * np.pi

mask = df['disp_prediction'] < 0
plt.hist(df.loc[mask, 'delta_delta'], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.3, label='sign(disp) < 0')
plt.hist(df.loc[~mask, 'delta_delta'], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.3, label='sign(disp) > 0')
plt.xlabel(r'$\mathrm{\Delta}\delta$ / rad')
plt.title(r'$\mathrm{\Delta}\delta$ between calculated $\delta$ and $\delta_{\mathrm{true}}$')
plt.legend()
plt.tight_layout()
if args.outputfile:
    plt.savefig(args.outputfile, dpi=300)
else:
    plt.show()

import numpy as np
import h5py
from tqdm import tqdm
import pandas as pd

from fact.io import read_h5py
from fact.analysis import (
    li_ma_significance,
    split_on_off_source_independent,
)
import click

columns = [
    'gamma_prediction',
    'theta_deg',
    'theta_deg_off_1',
    'theta_deg_off_2',
    'theta_deg_off_3',
    'theta_deg_off_4',
    'theta_deg_off_5',
]


@click.command()
@click.argument('data_path')
@click.option('--key', help='Key for the hdf5 group', default='events')
def main(data_path, key):

    events = read_h5py(data_path, key='events', columns=columns)

    theta2_cuts = np.arange(0.1, 0.0, -0.001)
    prediction_thresholds = np.arange(0.75, 1, 0.001)

    max_significance = 0
    selected = events
    for threshold in tqdm(prediction_thresholds):
        selected = selected.query('gamma_prediction >= {}'.format(threshold))

        theta2_on = selected.theta_deg**2
        theta2_off = pd.concat([
            selected['theta_deg_off_{}'.format(i)]
            for i in range(1, 6)
        ])**2

        for theta2_cut in theta2_cuts:
            theta2_on = theta2_on[theta2_on <= theta2_cut]
            theta2_off = theta2_off[theta2_off <= theta2_cut]

            n_on = len(theta2_on)
            n_off = len(theta2_off)

            sig = li_ma_significance(n_on, n_off, 0.2)
            if sig >= max_significance:
                max_significance = sig
                best_threshold = threshold
                best_theta2_cut = theta2_cut

    print('Threshold:', best_threshold)
    print('θ² cut:   ', best_theta2_cut)
    print('Li&Ma    :', max_significance)


if __name__ == '__main__':
    main()

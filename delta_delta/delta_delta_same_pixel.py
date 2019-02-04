import photon_stream as ps
from fact.plotting import camera, mark_pixel
import numpy as np
from fact.instrument.camera import get_neighbor_matrix, get_border_pixel_mask, get_pixel_coords
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from feature_stream import cleaning, calc_hillas_features_image, phs2image, calc_hillas_features_phs, facttools_cleaning, is_simulation_event
from tqdm import tqdm
import click
from astropy.table import Table
import warnings
from fact.io import read_data
from astropy.coordinates import SkyCoord, AltAz
from fact.instrument.constants import LOCATION
from fact.coordinates.utils import to_astropy_time
import pandas as pd
from fact.coordinates import horizontal_to_camera
from fact.instrument import camera_distance_mm_to_deg



picture_thresh = 5
boundary_thresh = 2
pix_x, pix_y = get_pixel_coords()


def calc_delta(image, mask):
    # covariance and eigenvalues/vectors for later calculations
    cov = np.cov(pix_x[mask], pix_y[mask], fweights=image[mask])
    eig_vals, eig_vecs = np.linalg.eigh(cov)

    # width, length and delta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        delta = np.arctan(eig_vecs[1, 1] / eig_vecs[0, 1])
    return delta


crab = SkyCoord.from_name('Crab Nebula')


def calc_delta_delta(event, mask):

    # if 'source_position_zd' not in df.columns:
    frame = AltAz(location=LOCATION, obstime=to_astropy_time(pd.to_datetime(event.observation_info.time)))

    crab_altaz = crab.transform_to(frame)

    source_position_zd_phs = crab_altaz.zen.deg
    source_position_az_phs = crab_altaz.az.deg

    source_x, source_y = horizontal_to_camera(
       source_position_zd_phs,
       source_position_az_phs,
       event.zd,
       event.az,
    )

    # safe x, y and t components of Photons. shape = (#photons,3)
    xyt = event.photon_stream.point_cloud
    lol = event.photon_stream.list_of_lists
    x, y, t = xyt.T
    x = np.rad2deg(x) / camera_distance_mm_to_deg(1)
    y = np.rad2deg(y) / camera_distance_mm_to_deg(1)

    cleaned_photons = np.zeros(len(x), dtype=bool)
    for i in range(len(lol)):
        if mask[i]:
            for j in range(len(lol[i])):
                cleaned_photons[i] = True

    cog_x = np.mean(x[cleaned_photons])
    cog_y = np.mean(y[cleaned_photons])

    true_delta = np.arctan2(cog_y - source_y, cog_x - source_x)

    delta = calc_delta(phs2image(event.photon_stream.list_of_lists), mask)
    delta_delta = true_delta - delta
    if delta_delta < -np.pi/2:
        delta_delta += 2 * np.pi
#    if delta_delta < -90:
#        delta_delta += 360

    return delta_delta


@click.command()
@click.argument('method', required=True)
@click.argument('path', required=True)
@click.argument('file', required=True)
@click.argument('feat', required=True)
@click.option('-n', '--number', default=300, type=int, help='Number of events to plot')
def main(method, path, file, feat, number):

    border_pix = get_border_pixel_mask()

    print("Reading in facttools dl1 file...")
    t = Table.read('/net/big-tank/POOL/projects/fact/photon-stream/facttools/crab/{}_dl1.fits'.format(file))
    # print("Reading in facttools dl2 file...")
    # dl2 = read_data('/home/ksedlaczek/Packages/open_crab_sample_analysis/dl2/crab.hdf5', key='events')
    print("Reading in PhotonStream data file...")
    reader = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/{}/{}.phs.jsonl.gz'.format(path, file))
    print("Reading in PhotonStream dl2 data...")
    ph_dl2 = read_data('/home/ksedlaczek/OUT_0.10/crab_data_precuts.hdf5', key='events')
    print("Done...")
    event = next(reader)

    d_delta = []
    disp_prediction = []

    for i in tqdm(range(number)):

        if is_simulation_event(event):
            event_num_phs = event.simulation_truth.event
            reuse_phs = event.simulation_truth.reuse
            run_id_phs = event.simulation_truth.run
        else:
            event_num_phs = event.observation_info.event
            reuse_phs = 42
            run_id_phs = event.observation_info.run

        if path != 'crab':
            run_id = file
            event_num = t[i]['MCorsikaEvtHeader.fEvtNumber']
            reuse = t[i]['MCorsikaEvtHeader.fNumReuse']
        else:
            run_id = file
            event_num = t[i]['EventNum']
            reuse = 42

        assert run_id != run_id_phs

        while (event_num_phs != event_num or reuse != reuse_phs):
            event = next(reader)
            if is_simulation_event(event):
                event_num_phs = event.simulation_truth.event
                reuse_phs = event.simulation_truth.reuse
                run_id_phs = event.simulation_truth.run
            else:
                event_num_phs = event.observation_info.event
                reuse_phs = 42
                run_id_phs = event.observation_info.run

        lol = event.photon_stream.list_of_lists
        # cut away unwanted time slices
        # lol = [[t for t in l if ((35 <= t) & (t < 75))] for l in lol]
        image = phs2image(lol)

        disp_prediction.append(ph_dl2.query('night == 20131104 & run == 162 & event == {}'.format(event_num))['disp_prediction'].values)
        cleaned_pix = t[i]['shower']
        d_delta.append(calc_delta_delta(event, cleaned_pix))

    mask = np.array([True if disp_prediction[i] < 0 else False for i in range(len(disp_prediction))])

    plt.figure()
    plt.hist(np.array(d_delta)[mask], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.4, label='disp < 0')
    plt.hist(np.array(d_delta)[~mask], bins=np.linspace(-np.pi/2, 3/2 * np.pi, 101), alpha=0.4, label='disp >= 0')
#    plt.hist(np.array(d_delta)[mask], bins=100, alpha=0.4, label='disp < 0')
#    plt.hist(np.array(d_delta)[~mask], bins=100, alpha=0.4, label='disp >= 0')


    plt.title(r'$\mathrm{\Delta}\delta$ between phs on facttools cleaning to $\delta_{\mathrm{true}}$ per image')
    plt.legend(loc='best')
    plt.xlabel(r'$\mathrm{\Delta}\delta$ / rad')
    plt.savefig('delta_true_diff_hist_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()



if __name__ == '__main__':
    main()

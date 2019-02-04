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
from IPython import embed


picture_thresh = 5
boundary_thresh = 2
x, y = get_pixel_coords()


def calc_delta(image, mask):
    # covariance and eigenvalues/vectors for later calculations
    cov = np.cov(x[mask], y[mask], fweights=image[mask])
    eig_vals, eig_vecs = np.linalg.eigh(cov)

    # width, length and delta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        delta = np.arctan(eig_vecs[1, 1] / eig_vecs[0, 1])
    return np.rad2deg(delta)


def calc_delta_perc(image, mask, perc=1):

    size = image[mask].sum()

    for i in range(1440):
        if mask[i] and (image[i] < size * perc/100):
            mask[i] = False

    # covariance and eigenvalues/vectors for later calculations
    cov = np.cov(x[mask], y[mask], fweights=image[mask])
    eig_vals, eig_vecs = np.linalg.eigh(cov)

    # width, length and delta
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        delta = np.arctan(eig_vecs[1, 1] / eig_vecs[0, 1])
    return np.rad2deg(delta)

crab = SkyCoord.from_name('Crab Nebula')


def calc_delta_delta(event, mask):

    x, y = get_pixel_coords()

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

    cleaned_pix = np.zeros(len(x), dtype=bool)
    for i in range(len(lol)):
        if mask[i]:
            for j in range(len(lol[i])):
                cleaned_pix[i] = True


    cog_x = np.mean(x[cleaned_pix])
    cog_y = np.mean(y[cleaned_pix])

    true_delta = np.arctan2(cog_y - source_y, cog_x - source_x)

    delta = calc_delta(phs2image(event.photon_stream.list_of_lists), mask)
    delta_delta = true_delta - delta
    if delta_delta < -np.pi/2:
        delta_delta += 2 * np.pi

    return delta_delta


@click.command()
@click.argument('method', required=True)
@click.argument('path', required=True)
@click.argument('file', required=True)
@click.argument('feat', required=True)
@click.option('-n', '--number', default=500, type=int, help='Number of events to plot')
def main(method, path, file, feat, number):

    border_pix = get_border_pixel_mask()

    print("Reading in facttools dl1 file...")
    t = Table.read('/net/big-tank/POOL/projects/fact/photon-stream/facttools/{}/{}_dl1.fits'.format(path, file))
    print("Reading in facttools dl2 file...")
    dl2 = read_data('/home/ksedlaczek/Packages/open_crab_sample_analysis/dl2/{}.hdf5'.format(path), key='events')
    print("Reading in PhotonStream data file...")
    reader = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/{}/{}.phs.jsonl.gz'.format(path, file))
    print("Done...")
    event = next(reader)

    all_pe_diff_mean = []
    all_pe_diff = []
    delta_delta = []
    delta_delta_diff = []
    delta_delta_diff_perc = []
    d_delta = []

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

        if method == 'DBSCAN':
            clustering = ps.photon_cluster.PhotonStreamCluster(event.photon_stream)
            if clustering.number < 1:
                continue
            biggest_cluster = np.argmax(np.bincount(clustering.labels[clustering.labels != -1]))
            mask = clustering.labels == biggest_cluster

            cleaned_pix_phs = np.zeros(len(image), dtype=bool)
            k = 0
            cleaned_img = np.zeros(len(image))
            for h in range(len(lol)):
                for j in range(len(lol[h])):
                    k += 1
                    if mask[k-1]:
                        cleaned_pix_phs[h] = True
        else:
            cleaned_pix_phs = facttools_cleaning(image, lol, picture_thresh, boundary_thresh)
            if sum(cleaned_pix_phs) < 1:
                continue

        cleaned_pix = t[i]['shower']

        t[i]['photoncharge'][t[i]['photoncharge'] < 0] = 0.0
        pe_difference = image - t[i]['photoncharge']

        if is_simulation_event(event):
            delta = np.rad2deg(dl2.query('night == 20131104 & run_id == 162 & event_num == {}'.format(t[i]['EventNum']))['delta'].values[0])
        else:
            delta = np.rad2deg(dl2.query('night == 20131104 & run_id == 162 & event_num == {}'.format(t[i]['EventNum']))['delta'].values[0])


        delta_diff_same = calc_delta(image, cleaned_pix) - delta
        delta_diff_whole = calc_delta(image, cleaned_pix_phs) - delta
        delta_diff_whole_perc = calc_delta_perc(image, cleaned_pix_phs) - delta
        for val in pe_difference:
            all_pe_diff.append(val)
        all_pe_diff_mean.append(np.mean(pe_difference))
        delta_delta.append(delta_diff_same)
        delta_delta_diff.append(delta_diff_whole)
        delta_delta_diff_perc.append(delta_diff_whole_perc)
        d_delta.append(calc_delta_delta(event, cleaned_pix))


    plt.figure()
    plt.hist(all_pe_diff_mean, bins=100, histtype='step')
    plt.title('Means of pe differences per image')
    plt.tight_layout()
    plt.savefig('means_hist_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(all_pe_diff, bins=100, histtype='step', density=True)
    plt.title('PE differences per pixel')
    plt.tight_layout()
    plt.savefig('diffs_hist_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(all_pe_diff, bins=100, histtype='step')
    plt.title('PE differences per pixel')
    plt.semilogy()
    plt.ylabel('events')
    plt.xlabel(r'$\symup{\Delta}$PE')
    plt.tight_layout()
    plt.savefig('diffs_hist_{}_{}_{}_logy.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(delta_delta, bins=70, histtype='step', density=True)
    plt.title('$\mathrm{\Delta}\delta$ between phs and facttools')
    plt.tight_layout()
    plt.savefig('delta_hist_same_pixels_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(delta_delta_diff, bins=70, histtype='step', density=True)
    plt.title(r'$\mathup{\Delta}\delta$ between phs and facttools')
    plt.xlabel(r'$\mathrm{\Delta}\delta / \textdegree$')
    plt.tight_layout()
    plt.savefig('delta_diff_hist_different_cleanings_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(delta_delta_diff_perc, bins=70, histtype='step', density=True)
    plt.title(r'$\mathup{\Delta}\delta$ between phs and facttools')
    plt.xlabel(r'$\mathrm{\Delta}\delta / \textdegree$')
    plt.tight_layout()
    plt.savefig('delta_diff_hist_perc_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(delta_delta_diff, bins=70, histtype='step', density=True, label='DBSCAN')
    plt.hist(delta_delta_diff_perc, bins=70, histtype='step', density=True, label='DBSCAN + > 1%')
    plt.title(r'$\mathup{\Delta}\delta$ between phs and facttools')
    plt.xlabel(r'$\mathrm{\Delta}\delta / \textdegree$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('delta_diff_hist_perc_and_normal_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()

    plt.figure()
    plt.hist(d_delta, bins=100, histtype='step', density=True)
    plt.title(r'$\mathrm{\Delta}\delta$ between phs on facttools cleaning to $\delta_{\mathrm{true}}$ per image')
    plt.tight_layout()
    plt.savefig('delta_true_diff_hist_{}_{}_{}.pdf'.format(method, feat, file))
    plt.clf()


if __name__ == '__main__':
    main()

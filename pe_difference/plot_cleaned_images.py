import photon_stream as ps
from fact.plotting import camera, mark_pixel
import numpy as np
from fact.instrument.camera import get_neighbor_matrix, get_border_pixel_mask, get_pixel_coords
import matplotlib as mpl
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
# from IPython import embed


mpl.rcParams['backend'] = 'pgf'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.size'] = 9
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['pgf.rcfonts'] = False
mpl.rcParams['pgf.texsystem'] = 'lualatex'
mpl.rcParams['pgf.preamble'] = '\input{/net/nfshome/home/ksedlaczek/phs_analysis/header-matplotlib.tex}'

picture_thresh = 5
boundary_thresh = 2
x, y = get_pixel_coords()

crab = SkyCoord.from_name('Crab Nebula')

@click.command()
@click.argument('method', required=True)
@click.argument('path', required=True)
@click.argument('file', required=True)
@click.argument('feat', required=True)
@click.option('-n', '--number', default=130, type=int, help='Number of events to plot')
def main(method, path, file, feat, number):

    border_pix = get_border_pixel_mask()

    print("Reading in facttools dl1 file...")
    t = Table.read('/net/big-tank/POOL/projects/fact/photon-stream/facttools/crab/{}_dl1.fits'.format(file))
    print("Reading in facttools dl2 file...")
    dl2 = read_data('/home/ksedlaczek/Packages/open_crab_sample_analysis/dl2/crab.hdf5', key='events')
    print("Reading in PhotonStream data file...")
    reader = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/{}/{}.phs.jsonl.gz'.format(path, file))
    print("Done...")
    event = next(reader)

    fig, axs = plt.subplots(3, 1, figsize=(4, 10), constrained_layout=True)

    plots = [camera(np.zeros(1440), cmap='inferno', ax=ax) for ax in axs]
    cbars = [fig.colorbar(plot, ax=ax) for plot, ax in zip(plots, axs)]
    plots[1].set_cmap('RdBu_r')

    with PdfPages('pe_difference_{}_{}.pdf'.format(feat, file)) as pdf:
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
                    break

            cleaned_pix = t[i]['shower']

            t[i]['photoncharge'][t[i]['photoncharge'] < 0] = 0.0
            pe_difference = image - t[i]['photoncharge']

            max_abs = np.max(np.abs(pe_difference))

            plots[0].set_array(image)
            plots[1].set_array(pe_difference)
            plots[2].set_array(t[i]['photoncharge'])

            plots[0].autoscale()
            plots[1].set_clim(-max_abs, max_abs)
            plots[2].autoscale()
            # mark_pixel(t[i]['shower'], color=(128/255, 186/255, 38/255), linewidth=2.5)
            for ax in axs:
                ax.axis('off')

            for cbar, plot in zip(cbars, plots):
                cbar.update_bruteforce(plot)

            # embed()
            if is_simulation_event(event):
                fig.suptitle('run {} event {} reuse {} mean {:.2f}'.format(run_id, event_num, reuse, np.mean(pe_difference)))
            else:
                # fig.suptitle('{} event {} mean {:.2f}'.format(file, event_num, np.mean(pe_difference)))
                fig.suptitle('{} event {}'.format(file[:8] + ' ' + file[9:], event_num))
            pdf.savefig(fig)


#    if method == "thresholds":
#        with PdfPages('cleaning_thresh_{}_{}.pdf'.format(feat, file)) as pdf:
#            for i in tqdm(range(number)):
#                fig = plt.figure()
#                ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
#                #ax.set_axis_off()
#
#                lol = event.photon_stream.list_of_lists
#                lol = [[t for t in l if ((35 <= t) & (t < 75))] for l in lol]
#                image = phs2image(lol)#, lower=30, upper=70)
#                cleaned_pix = facttools_cleaning(image, lol, 35, 75, picture_thresh, boundary_thresh)
#                with warnings.catch_warnings():
#                    warnings.simplefilter("ignore")
#                    arrival_times = np.array([np.nanmedian(l) for l in lol])
#                # cleaned_pix = cleaning(image, lol, picture_thresh, boundary_thresh)
#                if len(cleaned_pix[cleaned_pix != 0]) > 1:
#                    # border_ph = [(border_pix[i] and cleaned_pix[i]) for i in range(1440)]
#                    # leakage = image[border_ph].sum()/image[cleaned_pix].sum()
#                    df = calc_hillas_features_image(image, cleaned_pix)
#                    # ell = Ellipse(
#                    #     [df['cog_x'], df['cog_y']],
#                    #     df['length']*2,
#                    #     df['width']*2,
#                    #     angle=np.rad2deg(df['delta']),
#                    #     fill=False, linewidth=2, color='b'
#                    # )
#                    # ax.add_patch(ell)
#                    ell = Ellipse(
#                        [df['cog_x'], df['cog_y']],
#                        df['length']*4,
#                        df['width']*4,
#                        angle=np.rad2deg(df['delta']),
#                        fill=False, linewidth=1.5, color='b'
#                    )
#                    # ax.add_patch(ell)
#                    if is_simulation_event(event):
#                        fig.suptitle('run {} event {} reuse {}'.format(event.simulation_truth.run, event.simulation_truth.event, event.simulation_truth.reuse))
#                    else:
#                        fig.suptitle('{} event {} delta {}'.format(file, event.observation_info.event, df['delta']))
#                    if feat == 'arrival_times':
#                        with warnings.catch_warnings():
#                            warnings.simplefilter("ignore")
#                            x = arrival_times-np.nanmean(arrival_times)
#                        x[np.isnan(x)] = 0
#                        c = camera(x, cmap='Spectral', ax=ax)
#                        mark_pixel(cleaned_pix, color='k', linewidth=2.5)
#                    else:
#                        c = camera(image, cmap='viridis', ax=ax)
#                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
#                    ax.axis('off')
#                    fig.colorbar(c)
#                    pdf.savefig(fig)
#                    ax.cla()
#                plt.close(fig)
#
#    if method == "DBSCAN":
#        with PdfPages('cleaning_DBSCAN_biggest_{}_{}.pdf'.format(feat, file)) as pdf:
#            for i in tqdm(range(number)):
#                fig = plt.figure()
#                ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
#                event = next(reader)
#
#                # clustering of events
#                clustering = ps.photon_cluster.PhotonStreamCluster(event.photon_stream, eps=0.08)
#
#                if clustering.number > 0:
#
#                    lol = event.photon_stream.list_of_lists
#                    image = phs2image(lol)
#                    with warnings.catch_warnings():
#                        warnings.simplefilter("ignore")
#                        arrival_times = np.array([np.nanmedian(l) for l in lol])
#                    # biggest cluster:
#                    biggest_cluster = np.argmax(np.bincount(clustering.labels[clustering.labels != -1]))
#
#                    mask = clustering.labels == biggest_cluster
#                    # mask = clustering.labels != -1
#
#                    xyt = event.photon_stream.point_cloud
#                    x, y, t = xyt.T
#                    cleaned_pix = np.zeros(len(image), dtype=bool)
#                    k = 0
#                    cleaned_img = np.zeros(len(image))
#                    for i in range(len(lol)):
#                        for j in range(len(lol[i])):
#                            k += 1
#                            if mask[k-1]:
#                                cleaned_pix[i] = True
#                                cleaned_img[i] += 1
#
#                    df = calc_hillas_features_phs(event.photon_stream, clustering)
#                    # ell = Ellipse(
#                    #     [df['cog_x'], df['cog_y']],
#                    #     df['length']*2,
#                    #     df['width']*2,
#                    #     angle=np.rad2deg(df['delta']),
#                    #     fill=False, linewidth=2, color='b'
#                    # )
#                    # ax.add_patch(ell)
#                    ell = Ellipse(
#                        [df['cog_x'], df['cog_y']],
#                        df['length']*4,
#                        df['width']*4,
#                        angle=np.rad2deg(df['delta']),
#                        fill=False, linewidth=1.5, color='b'
#                    )
#                    # ax.add_patch(ell)
#                    if is_simulation_event(event):
#                        fig.suptitle('run {} event {} reuse {}'.format(event.simulation_truth.run, event.simulation_truth.event, event.simulation_truth.reuse))
#                    else:
#                        fig.suptitle('{} event {} delta {}'.format(file, event.observation_info.event, df['delta']))
#                    if feat == 'arrival_times':
#                        with warnings.catch_warnings():
#                            warnings.simplefilter("ignore")
#                            x = arrival_times-np.nanmean(arrival_times)
#                        c = camera(x, cmap='viridis', ax=ax)
#                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
#                    else:
#                        c = camera(image, cmap='viridis', ax=ax)
#                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
#                    ax.axis('off')
#                    fig.colorbar(c)
#                    pdf.savefig(fig)
#                    ax.cla()
#                plt.close(fig)


if __name__ == '__main__':
    main()

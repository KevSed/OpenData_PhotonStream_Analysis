import photon_stream as ps
from fact.plotting import camera, mark_pixel
import numpy as np
from fact.instrument.camera import get_neighbor_matrix, get_border_pixel_mask
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from feature_stream import cleaning, calc_hillas_features_image, phs2image, calc_hillas_features_phs, facttools_cleaning, is_simulation_event
from tqdm import tqdm
import click
from astropy.table import Table
import warnings
from fact.io import read_data
# from IPython import embed


picture_thresh = 5
boundary_thresh = 2


@click.command()
@click.argument('method', required=True)
@click.argument('path', required=True)
@click.argument('file', required=True)
@click.argument('feat', required=True)
@click.option('-n', '--number', default=100, type=int, help='Number of events to plot')
def main(method, path, file, feat, number):

    border_pix = get_border_pixel_mask()
    if method == "thresholds":
        reader = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/{}/{}.phs.jsonl.gz'.format(path, file))
        with PdfPages('cleaning_thresh_{}_{}.pdf'.format(feat, file)) as pdf:
            for i in tqdm(range(number)):
                fig = plt.figure()
                ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
                #ax.set_axis_off()
                event = next(reader)

                lol = event.photon_stream.list_of_lists
                lol = [[t for t in l if ((35 <= t) & (t < 75))] for l in lol]
                image = phs2image(lol)#, lower=30, upper=70)
                cleaned_pix = facttools_cleaning(image, lol, 35, 75, picture_thresh, boundary_thresh)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    arrival_times = np.array([np.nanmedian(l) for l in lol])
                # cleaned_pix = cleaning(image, lol, picture_thresh, boundary_thresh)
                if len(cleaned_pix[cleaned_pix != 0]) > 1:
                    # border_ph = [(border_pix[i] and cleaned_pix[i]) for i in range(1440)]
                    # leakage = image[border_ph].sum()/image[cleaned_pix].sum()
                    df = calc_hillas_features_image(image, cleaned_pix)
                    # ell = Ellipse(
                    #     [df['cog_x'], df['cog_y']],
                    #     df['length']*2,
                    #     df['width']*2,
                    #     angle=np.rad2deg(df['delta']),
                    #     fill=False, linewidth=2, color='b'
                    # )
                    # ax.add_patch(ell)
                    ell = Ellipse(
                        [df['cog_x'], df['cog_y']],
                        df['length']*4,
                        df['width']*4,
                        angle=np.rad2deg(df['delta']),
                        fill=False, linewidth=1.5, color='b'
                    )
                    # ax.add_patch(ell)
                    if is_simulation_event(event):
                        fig.suptitle('run {} event {} reuse {}'.format(event.simulation_truth.run, event.simulation_truth.event, event.simulation_truth.reuse))
                    else:
                        fig.suptitle('{} event {} delta {}'.format(file, event.observation_info.event, df['delta']))
                    if feat == 'arrival_times':
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            x = arrival_times-np.nanmean(arrival_times)
                        x[np.isnan(x)] = 0
                        c = camera(x, cmap='Spectral', ax=ax)
                        mark_pixel(cleaned_pix, color='k', linewidth=2.5)
                    else:
                        c = camera(image, cmap='viridis', ax=ax)
                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
                    ax.axis('off')
                    fig.colorbar(c)
                    pdf.savefig(fig)
                    ax.cla()
                plt.close(fig)

    if method == "DBSCAN":
        reader = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/{}/{}.phs.jsonl.gz'.format(path, file))
        with PdfPages('cleaning_DBSCAN_biggest_{}_{}.pdf'.format(feat, file)) as pdf:
            for i in tqdm(range(number)):
                fig = plt.figure()
                ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
                event = next(reader)

                # clustering of events
                clustering = ps.photon_cluster.PhotonStreamCluster(event.photon_stream)

                if clustering.number > 0:

                    lol = event.photon_stream.list_of_lists
                    image = phs2image(lol)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        arrival_times = np.array([np.nanmedian(l) for l in lol])
                    # biggest cluster:
                    biggest_cluster = np.argmax(np.bincount(clustering.labels[clustering.labels != -1]))

                    mask = clustering.labels == biggest_cluster
                    # mask = clustering.labels != -1

                    xyt = event.photon_stream.point_cloud
                    x, y, t = xyt.T
                    cleaned_pix = np.zeros(len(image), dtype=bool)
                    k = 0
                    cleaned_img = np.zeros(len(image))
                    for i in range(len(lol)):
                        for j in range(len(lol[i])):
                            k += 1
                            if mask[k-1]:
                                cleaned_pix[i] = True
                                cleaned_img[i] += 1

                    cleaned_pix_perc = np.zeros(1440, dtype=bool)
                    for i in range(1440):
                        if cleaned_pix[i] and (cleaned_img[i] > mask.sum() / 200):
                            cleaned_pix_perc[i] = True
                    df = calc_hillas_features_phs(event.photon_stream, clustering)
                    # ell = Ellipse(
                    #     [df['cog_x'], df['cog_y']],
                    #     df['length']*2,
                    #     df['width']*2,
                    #     angle=np.rad2deg(df['delta']),
                    #     fill=False, linewidth=2, color='b'
                    # )
                    # ax.add_patch(ell)
                    ell = Ellipse(
                        [df['cog_x'], df['cog_y']],
                        df['length']*4,
                        df['width']*4,
                        angle=np.rad2deg(df['delta']),
                        fill=False, linewidth=1.5, color='b'
                    )
                    # ax.add_patch(ell)
                    if is_simulation_event(event):
                        fig.suptitle('run {} event {} reuse {}'.format(event.simulation_truth.run, event.simulation_truth.event, event.simulation_truth.reuse))
                    else:
                        fig.suptitle('{} event {} delta {:.2f}'.format(file, event.observation_info.event, np.rad2deg(df['delta'])))
                    if feat == 'arrival_times':
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            x = arrival_times-np.nanmean(arrival_times)
                        c = camera(x, cmap='viridis', ax=ax)
                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
                    else:
                        c = camera(image, cmap='viridis', ax=ax)
                        mark_pixel(cleaned_pix, color=(128/255, 186/255, 38/255), linewidth=2.5)
                        mark_pixel(cleaned_pix_perc, color='red', linewidth=1.5)
                    ax.axis('off')
                    fig.colorbar(c)
                    pdf.savefig(fig)
                    ax.cla()
                plt.close(fig)


    if method == "facttools":
        print('facttools')
        with PdfPages('cleaning_facttools_{}_{}.pdf'.format(feat, file)) as pdf:
            t = Table.read('/net/big-tank/POOL/projects/fact/photon-stream/facttools/{}/{}_dl1.fits'.format(path, file))
            dl2 = read_data('/home/ksedlaczek/Packages/open_crab_sample_analysis/dl2/crab.hdf5', key='events')

            for i in tqdm(range(number)):
                fig = plt.figure()
                ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
#                if path != 'crab':
#                    fig.suptitle('run {} event {} reuse {}'.format(file, t[i]['MCorsikaEvtHeader.fEvtNumber'], t[i]['MCorsikaEvtHeader.fNumReuse']))
#                else:
#
#                    fig.suptitle('{} event {} delta {:.4f}'.format(file, t[i]['EventNum'], dl2.query('night == 20131104 & run_id == 162 & event_num == {}'.format(t[i]['EventNum']))['delta'].values[0]))

                t[i]['photoncharge'][t[i]['photoncharge'] < 0] = 0.0
                if feat == 'arrival_times':
                    c = camera(t[i]['arrivalTime']-t[i]['arrivalTime'].mean(), cmap='Spectral', ax=ax)
                    # mark_pixel(t[i]['shower'], color='k', linewidth=2.5)
                    ax.axis('off')
                    cb = fig.colorbar(c)
                    cb.set_label(label=r'$t-\bar{t}$ / ns', fontsize=16)
                else:
                    c = camera(t[i]['photoncharge'], cmap='viridis', ax=ax)
                    ax.axis('off')
                    cb = fig.colorbar(c)
                    cb.set_label(label=r'Number of Photons', fontsize=16)

                    #mark_pixel(t[i]['shower'], color=(128/255, 186/255, 38/255), linewidth=2.5)
                # mark_pixel(t[i]['shower'], color=(128/255, 186/255, 38/255), linewidth=2.5)
                pdf.savefig(fig)
                ax.cla()
                plt.close(fig)



if __name__ == '__main__':
    main()

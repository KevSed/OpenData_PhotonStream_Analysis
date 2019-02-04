import photon_stream as ps
import numpy as np
import matplotlib.pyplot as plt
import click


@click.command()
@click.argument('file', required=True, type=click.Path(exists=True))
@click.argument('ev', required=True, type=int)
def main(file, ev):

    reader = ps.EventListReader(file)
    for event in reader:
        if event.observation_info.event == ev:
            break

    pc = event.photon_stream.point_cloud
    pc[:, 0] = np.rad2deg(pc[:, 0])
    pc[:, 1] = np.rad2deg(pc[:, 1])
    pc[:, 2] *= 1e9*0.35

    min_time = pc[:, 2].min()
    max_time = pc[:, 2].max()
    max_space = pc[:, 0].max()

    fig = plt.figure(figsize=[6.4, 8.6])
    ax1 = fig.add_subplot(121)
    ax1.set_xlim([-1.01*max_space, 1.01*max_space])
    ax1.set_ylim([0.99*min_time, 1.01*max_time])
    ax1.set_aspect('equal')
    ax1.scatter(pc[:, 0], pc[:, 2], lw=0, alpha=0.075, s=35., color='b')
    ax1.set_xlabel(r'$c_{x}$ / $\degree$')
    ax1.set_ylabel(r'$c_{t}$ / $\degree$')

    ax2 = fig.add_subplot(122)
    ax2.set_xlim([-1.01*max_space, 1.01*max_space])
    ax2.set_ylim([0.99*min_time, 1.01*max_time])
    ax2.set_aspect('equal')
    ax2.scatter(pc[:, 1], pc[:, 2], lw=0, alpha=0.075, s=35., color='b')
    ax2.set_xlabel(r'$c_{y}$ / $\degree$')
    ax2.set_ylabel(r'$c_{t}$ / $\degree$')

    plt.savefig('/net/nfshome/home/ksedlaczek/phs_analysis/point_cloud.pdf')


if __name__=='__main__':
    main()

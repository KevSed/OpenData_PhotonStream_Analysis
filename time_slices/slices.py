# coding: utf-8
import photon_stream as ps
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import os
from IPython import embed


size_cut = 0
photons_per_pixel = 0
proton_path = '/net/big-tank/POOL/projects/fact/photon-stream/stream_data/proton/uwe/'
gamma_path = '/net/big-tank/POOL/projects/fact/photon-stream/stream_data/gamma/gustav/'

def plot_time(lol2, event):
    lol = event.photon_stream.list_of_lists
    clustering = ps.PhotonStreamCluster(event.photon_stream)
    if clustering.number > 0:
        biggest_cluster = np.argmax(np.bincount(clustering.labels[clustering.labels != -1]))
        mask = clustering.labels == biggest_cluster
        if mask.sum() > size_cut:
            for i in range(len(lol)):
                if len(lol[i]) > photons_per_pixel:
                    for k in range(len(lol[i])):
                        lol2.append(lol[i][k])


def safe_plot(lol_d, lol_g, lol_p):
    plt.figure()
    plt.hist(lol_d, bins=np.arange(0,150,1)-0.5, histtype='step', density=True, label='data')
    plt.hist(lol_g, bins=np.arange(0,150,1)-0.5, histtype='step', density=True, label='gamma')
    plt.hist(lol_p, bins=np.arange(0,150,1)-0.5, histtype='step', density=True, label='proton')
    plt.legend()
    plt.xlabel('time slices')
    plt.tight_layout()
    plt.savefig('all_slices_min_{}_per_pixel.pdf'.format(photons_per_pixel))
    plt.clf()


def main():

    gammas = list(filter(lambda a: a.endswith('phs.jsonl.gz'), os.listdir(gamma_path)))
    protons = list(filter(lambda a: a.endswith('phs.jsonl.gz'), os.listdir(proton_path)))

    data = ps.EventListReader('/net/big-tank/POOL/projects/fact/photon-stream/stream_data/crab/20131101_160.phs.jsonl.gz')

    print('Reading data...')
    lol_d = []
    events = 0
    for event in tqdm(data, total=2000):
        events += 1
        plot_time(lol_d, event)
        if events > 2000:
            break

    print('Reading gamma...')
    lol_g = []
    events = 0
    for file in gammas:
        gamma = ps.EventListReader(gamma_path + file)
        for event in tqdm(gamma, total=2000):
            events += 1
            plot_time(lol_g, event)
            if events > 2000:
                break
        if events > 2000:
            break

    print('Reading proton...')
    lol_p = []
    events = 0
    for file in protons:
        proton = ps.EventListReader(proton_path + file)
        for event in tqdm(proton, total=2000):
            events += 1
            plot_time(lol_p, event)
            if events > 2000:
                break
        if events > 2000:
            break

    safe_plot(lol_d, lol_g, lol_p)


if __name__ == '__main__':
    main()

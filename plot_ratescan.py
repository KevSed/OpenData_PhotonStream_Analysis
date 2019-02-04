from fact.io import read_data
import numpy as np
import matplotlib.pyplot as plt


size_cuts = np.logspace(0, 5, 500)


for eps in [0.03, 0.05, 0.06, 0.07, 0.1, None]:
    if eps:
        d = read_data(f'/net/big-tank/POOL/projects/fact/photon-stream/features/{eps:.2f}/crab_data.hdf5', key='events', columns=['size'])
    else:
        d = read_data(f'/net/big-tank/POOL/projects/fact/data/open/dl2/FACT-Tools/v1.1.2/open_crab_sample_facttools_dl2.hdf5', key='events', columns=['size'])

    n_events = []
    size = d['size'].values
    for cut in size_cuts:
        size = size[size >= cut]
        n_events.append(len(size))

    plt.plot(size_cuts, n_events, label=r'$\varepsilon =$ {:.2}'.format(eps) if eps else 'FACT-Tools')

plt.yscale('log')
plt.xscale('log')
plt.ylabel('events')
plt.xlabel('size-cut')
plt.grid()
plt.title(r'Data rates for different $\varepsilon$')
plt.legend()
plt.tight_layout()
plt.savefig('data_rates_eps.pdf')

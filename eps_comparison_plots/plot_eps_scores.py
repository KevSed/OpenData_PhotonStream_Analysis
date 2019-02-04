import matplotlib.pyplot as plt
import numpy as np


def main():
    df = {}
    df['eps'] = [0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10]
    df['eps_sigma'] = [0.03, 0.05, 0.06, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10]
    df['AUC'] = [0.8928, 0.8937, 0.8931, 0.8908, 0.8842, 0.8821, 0.8757]
    df['sigma'] = [5.6, 16.2, 18.2, 22.8, 20.8, 22.3, 21.9, 22.3, 21.6, 19.8]
    df['R2_disp'] = [0.589, 0.588, 0.597, 0.616, 0.618, 0.615, 0.610]
    df['ACC'] = [0.8446, 0.8483, 0.8519, 0.8488, 0.8422, 0.8390, 0.8319]
    df['R2_E'] = [0.673, 0.674, 0.6781, 0.6817, 0.6834, 0.6847, 0.6849]
    df['theta_cut'] = [0.1, 0.099, 0.1, 0.1, 0.099, 0.096, 0.093]

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    plt.plot(df['eps'], df['AUC'], 'bx', label='AUC')
    plt.plot(df['eps'], df['ACC'], 'kx', label='sgn(disp) ACC')
    plt.plot(df['eps'], df['R2_E'], 'rx', label='Energy R2 score')
    plt.plot(df['eps'], df['R2_disp'], 'gx', label='|disp| R2 score')

    ax.set_xticks(df['eps'])
    ax.set_yticks(np.arange(0, 1, 0.1))
    ax.set_yticks(np.arange(0, 1, 0.05), minor=True)

    # And a corresponding grid
    ax.grid(which='both')

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.5)
    ax.grid(which='major', alpha=0.7)

    # plt.axis([min(df['eps'])-0.002, max(df['eps'])+0.002, 0, 1])
    plt.xlabel(r'$\varepsilon$')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('eps_scores.pdf')
    plt.clf()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
#     ax.set_yticklabels(np.arange(0, 22, 2))
    print(df['eps_sigma'])
    plt.plot(df['eps_sigma'], df['sigma'], 'kx', label='significance Li\&Ma')
    plt.xlim([0.028, 0.102])
    plt.ylim([0.0, 25])
    ax.grid()
    plt.xlabel(r'$\varepsilon$')
    ax.set_ylabel(r'$\sigma$')
    plt.title(r'$\sigma$ Li\&Ma')
    plt.tight_layout()
    plt.savefig('eps_sigma.pdf')
    plt.clf()

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_ylim(0, 0.102)
    # ax.set_yticklabels(np.arange(0.0, 0.11, 0.01))
    plt.plot(df['eps'], df['theta_cut'], 'kx')
    plt.xlim([0.068, 0.102])
    ax.grid()
    plt.xlabel(r'$\varepsilon$')
    ax.set_ylabel(r'$\theta^{2}$ cut')
    plt.title(r'$\theta^{2}$ cut')
    plt.tight_layout()
    plt.savefig('eps_theta_cut.pdf')
    plt.clf()

if __name__=='__main__':
    main()

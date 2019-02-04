from h5py import File
import matplotlib.pyplot as plt
import click


@click.command()
@click.argument('file', required=True, type=click.Path(exists=True))
@click.argument('out', required=True, type=click.Path(exists=False))
@click.argument('data', required=True)
@click.option('-t', '--threshold', default=0.0, type=float, help='threshold for GAMMA prediction')
def main(file, out, data, threshold):
    print("Loading File...")
    f = File(file)
    print("Plotting {} skymap...".format(data))

    plt.figure()
    plt.hist2d(f['events']['source_x_prediction'][f['events']['gamma_prediction'][:]>threshold],
            f['events']['source_y_prediction'][f['events']['gamma_prediction'][:]>threshold],
            bins=100)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('{} source predictions'.format(data))
    plt.colorbar()
    plt.savefig(out)
    plt.cla()
    plt.clf()


if __name__=='__main__':
    main()

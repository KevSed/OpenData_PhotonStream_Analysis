from fact.io import read_data, to_h5py
import click


@click.command()
@click.argument('file', type=click.Path(exists=True))
def main(file):
    run_meta = read_data('~/phs_analysis/open_crab_sample_runs.csv')
    to_h5py(run_meta.iloc[:], file, key='runs', mode='a')


if __name__=='__main__':
    main()

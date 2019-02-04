from features_from_phs import cluster_labels, is_simulation_file
from fact.io import to_h5py
import click
from tqdm import tqdm
from multiprocessing import Pool, cpu_count


@click.command()
@click.argument('output_file', type=click.Path(exists=False))
@click.argument('input_file', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('-n', '--n-jobs', default=-1, type=int, help='Number of cores to use')
def main(output_file, input_file, n_jobs):

    if n_jobs == -1:
        n_jobs = 15

    print('Calculating features using', n_jobs, 'cores')

    if is_simulation_file(input_file[0]):
        print('Received simulation files as input.')
    else:
        print('Received data files as input.')

    with Pool(n_jobs) as pool:
        results = pool.imap_unordered(cluster_labels, input_file)
        for df in tqdm(results, total=len(input_file)):
            to_h5py(df, output_file, key="events", mode='a', index=False)


if __name__ == '__main__':
    main()

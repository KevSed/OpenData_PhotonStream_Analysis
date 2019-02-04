from features_from_phs import gen_features, is_simulation_file
from fact.io import to_h5py
import click
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from time import sleep


def as_completed(futures):
    futures = list(futures)
    while futures:
        for f in futures.copy():
            if f.ready():
                futures.remove(f)
                yield f.get()
        sleep(0.1)


@click.command()
@click.argument('output_file', type=click.Path(exists=False)
@click.argument('input_file', nargs=-1, required=True, type=click.Path(exists=True))
@click.option('-e', '--eps', default=0.1, type=float, help='Eps parameter for DBSCAN')
@click.option('-n', '--n-jobs', default=-1, type=int, help='Number of cores to use')
def main(output_file, input_file, eps, n_jobs):

    if n_jobs == -1:
        n_jobs = 48 # cpu_count()

    print('Calculating features using', n_jobs, 'cores')

    if is_simulation_file(input_file[0]):
        print('Received simulation files as input.')
    else:
        print('Received data files as input.')

    with Pool(n_jobs) as pool:
        results = [pool.apply_async(gen_features, kwds={'data_file': f}) for f in input_file]
        for df in tqdm(as_completed(results), total=len(input_file)):
            to_h5py(df, output_file, key="events", mode='a', index=False)


if __name__ == '__main__':
    main()

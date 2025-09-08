import os

import pandas as pd


def read_pool_file(fname):
    df = pd.read_csv(fname, usecols=['Pool_name', 'Sequence'])
    df = df[~df['Sequence'].str.contains('-')]
    df = df[df['Sequence'].str.len() <= 30]
    df = df.drop_duplicates('Sequence').reset_index(drop=True)
    return df


def read_pool_files(fnames):
    return pd.concat([read_pool_file(f) for f in fnames]).drop_duplicates('Sequence').reset_index(drop=True)


def read_pool_dir(dir_path, ext='.txt'):
    all_files = [os.path.join(dir_path, f) for f in os.listdir(dir_path) if f.endswith(ext)]
    return read_pool_files(all_files)


def write_db(db, out_path):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    db['Sequence'].to_csv(out_path, index=False, header=False)

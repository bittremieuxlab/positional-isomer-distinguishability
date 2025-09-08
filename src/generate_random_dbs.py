import os

import pandas as pd

from src.utils import read_pool_file, read_pool_dir, write_db


def generate_random_db(dir_path="data", ext=".txt", size=1000):
    df = read_pool_dir(dir_path, ext)
    return df.sample(n=int(size))

def generate_random_dbs(dir_path="data", ext=".txt", sizes=None):
    return [generate_random_db(dir_path, ext, s) for s in sizes]


if __name__ == '__main__':
    random_dbs = generate_random_dbs("data", sizes=[1e3, 1e4, 1e5])
    for db in random_dbs:
        write_db(db, out_path=f"random_dbs/{len(db)}.txt")

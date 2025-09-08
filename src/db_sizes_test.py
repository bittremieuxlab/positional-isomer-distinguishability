import time

import pandas as pd

from generate_random_dbs import generate_random_dbs
from msci_runner import find_indistinguishable_peptides
from utils import write_db

if __name__ == "__main__":
    sizes = [1e3, 2e3, 4e3, 6e3, 8e3, 1e4]
    # random_dbs = generate_random_dbs("data", sizes=sizes)
    # for db in random_dbs:
    #     t = time.time()
    #     print(f"Testing a db size: {len(db)}")
    #     pep_file = f"random_dbs/{len(db)}.txt"
    #     write_db(db, out_path=pep_file)
    #     find_indistinguishable_peptides(pep_file)
    #     print(f"Size {len(db)} took {time.time() - t:.2f} seconds")

    for s in sizes:
        t = time.time()
        print(f"Testing a db size: {s}")
        db = pd.read_csv("data/generated_positional_isomers.csv", header=None)
        db = db.sample(n=int(s))
        pep_file = f"random_dbs/pos_iso_{s}.txt"
        write_db(db, out_path=pep_file, col_name=0)
        find_indistinguishable_peptides(pep_file)
        print(f"Size {s} took {time.time() - t:.2f} seconds")

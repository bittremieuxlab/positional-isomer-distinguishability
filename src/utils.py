import pandas as pd
def read_pool_file(fname):
    return pd.read_csv(fname, usecols=['Pool_name', 'Sequence'])

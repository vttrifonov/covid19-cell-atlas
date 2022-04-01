#%%
from pathlib import Path
import xarray as xa
import numpy as np
import pandas as pd

#%%
class _config:
    project='covid19-cell-atlas'
    cache = Path.home()/'.cache'/project
    root = Path(__file__).parent

config = _config()

def xa_dense(x):
    if isinstance(x, xa.DataArray):
        if hasattr(x.data, 'todense'):
            x = x.copy()
            x.data = x.data.todense()
    else:
        for k in x:
            x[k] = xa_dense(x[k])
    return x


xa.DataArray.todense = xa_dense
xa.Dataset.todense = xa_dense

def to_series_sparse(x):
    d = x.data
    if not hasattr(d, 'coords'):
        if hasattr(d, 'tocoo'):
            d = d.tocoo()
        else:
            raise ValueError('not sparse')
    if not hasattr(d, 'data'):
        raise ValueError('not sparse')
    
    s = pd.DataFrame(d.coords.T, columns=x.dims)
    s[x.name] = d.data
    for k in x.dims:
        s[k] = x[k].data[s[k]]
    s = s.set_index(list(x.dims))
    return s
xa.DataArray.to_series_sparse = to_series_sparse

def xa_mmult(x, y, dim_x, dim_y):
    return xa.apply_ufunc(
        np.matmul, x, y,
        input_core_dims=[dim_x, dim_y],
        output_core_dims=[[dim_x[0], dim_y[1]]],
        join='inner'
    )

def round_float(x, n=1):
    e = np.floor(np.log10(x))
    m = np.ceil(x*10**(n-e))/10**n
    return m*10**e

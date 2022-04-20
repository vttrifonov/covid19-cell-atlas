#%%
from pathlib import Path
import xarray as xa
import numpy as np
import pandas as pd
import sparse

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
    s, x = np.sign(x), np.abs(x)
    e = np.floor(np.log10(x))
    m = np.ceil(x*10**(n-e))/10**n
    return s*m*10**e

def dataarray_from_series(x, fill_value=None):
    i = x.index.to_frame(index=False)
    for c in i:
        if not isinstance(i[c], pd.CategoricalDtype):
            i[c] = i[c].astype('category')

    d = sparse.COO(
        [i[c].cat.codes.to_numpy() for c in i],
        x.to_numpy(),
        shape = tuple(len(i[c].cat.categories) for c in i),
        fill_value = fill_value
    )
    d = xa.DataArray(
        d,
        coords=[(c, i[c].cat.categories) for c in i],
        name=x.name
    )
    return d


def loess(x, y, w = None, t = None):
    from rpy2.robjects import r, Formula, pandas2ri, numpy2ri, default_converter
    from rpy2.robjects.conversion import localconverter
    import rpy2.rinterface as ri
    r_loess = r['loess']

    if w is None:
        w = np.ones_like(x)

    if t is None:
        t = x

    with localconverter(
        default_converter+pandas2ri.converter+numpy2ri.converter
    ):
        try:
            f = r_loess(
                Formula('y~x'), 
                pd.DataFrame({'x': x, 'y': y}), w,
            )
            f = r['predict'](f, pd.DataFrame({'x': t}))
        except ri.embedded.RRuntimeError:
            f = np.full(len(t), np.nan)
        return f
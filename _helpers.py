#%%
from pathlib import Path
import xarray as xa

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
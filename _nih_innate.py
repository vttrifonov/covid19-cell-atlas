
#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import scanpy
import xarray as xa
import numpy as np
import sparse
from .helpers import config
from .common.caching import lazy, compose, XArrayCache

#%%
def _data1():
    d = scanpy.read(
        config.cache/'data'/'nih-innate.h5ad',
        backed='r'
    )
    return d

def _data2():
    d = scanpy.read(
        config.cache/'data'/'nih-innate.h5ad'
    )
    return d

class _nih_innate:
    pass

def _():
    _nih_innate.storage = config.cache/'nih-innate'

    @compose(property, lazy)
    def data(self):
        return _data1()
    _nih_innate.data = data

    @compose(property, lazy)
    def X(self):
        import dask.array as da
        import numpy as np
        import sparse

        dtype = self.data.X.dtype
        shape = self.data.shape
        chunks = (shape[0]//31, shape[1]//1)
        def _(x, n):
            s = x//n
            c = (n,)*s
            s = x-s*n
            if s>0:
                c = c+(s,)                
            return c
        chunks = [_(x, n) for x, n in zip(shape, chunks)]

        path = self.storage/'X'
        if not path.exists():
            path.mkdir(parents=True, exist_ok=False)
            X = _data2().X
        def _(block_info=None):
            bi = block_info[None]
            al = bi['array-location']
            file = path/f'{al[0][0]}_{al[1][0]}.npz'
            if file.exists():
                x = sparse.load_npz(file)
                return x
            x = X[slice(*al[0]), slice(*al[1])]
            x = sparse.COO.from_scipy_sparse(x)
            sparse.save_npz(file, x)
            return x
        x = da.map_blocks(_, chunks=chunks, meta=sparse.zeros((0,0), dtype=dtype))

        return x
    _nih_innate.X = X

    @compose(property, lazy, XArrayCache())
    def var_ensembl(self):
        from .ensembl import annotations
        display = annotations.gene_display(db='homo_sapiens_core_104_38') 
        display = display[['stable_id', 'display_label']]

        annot = annotations.gene_entrez(db='homo_sapiens_core_104_38')
        annot = annot[['stable_id', 'dbprimary_acc']]

        annot = annot.merge(display, on='stable_id')
        annot = annot.rename(columns={
            'dbprimary_acc': 'entrez',
            'display_label': 'var'
        })
        annot = annot[annot['var'].isin(self.data.var_names)]
        annot = annot.set_index('var').to_xarray()
        return annot
    _nih_innate.var_ensembl = var_ensembl
_()

nih_innate = _nih_innate()

#%%
if __name__=='__main__':
    self = nih_innate

#%%
    import sparse
    x1 = self.X.map_blocks(
        lambda x: x.sum(keepdims=True), 
        meta=sparse.zeros((0,0), dtype=self.X.dtype)
    )
    x1 = x1.compute()

# %%

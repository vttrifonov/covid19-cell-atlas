
#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import scanpy
import xarray as xa
import numpy as np
import sparse
from .helpers import config
from .common.caching import lazy, compose, PickleCache

#%%
def _data1():
    d = scanpy.read(
        config.cache/'data'/'nih-innate.h5ad',
        backed='r'
    )
    return d

class _dataset1:
    pass

def _():
    _dataset1.storage = config.cache/'nih-innate'

    @compose(property, lazy)
    def data(self):
        return _data1()
    _dataset1.data = data

    @compose(property, lazy)
    def X(self):
        import dask.array as da
        import numpy as np
        import sparse

        dtype = self.data.X.dtype
        shape = self.data.shape
        chunks = (shape[0]//39, shape[1]//1)
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
        def _(block_info=None):
            bi = block_info[None]
            al = bi['array-location']
            file = path/f'{al[0][0]}_{al[1][0]}.npz'
            if file.exists():
                x = sparse.load_npz(file)
                return x
            data = _data1()
            x = data.chunk_X(select=np.arange(al[0][0], al[0][1]))
            x = x[:, slice(*al[1])]
            x = sparse.COO.from_numpy(x, idx_dtype=np.int32)
            sparse.save_npz(file, x)
            return x
        x = da.map_blocks(_, chunks=chunks, meta=sparse.zeros((0,0), dtype=dtype))

        return x
    _dataset1.X = X

    @compose(property, lazy, PickleCache())
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
    _dataset1.var_ensembl = var_ensembl

    @compose(property, lazy, PickleCache())
    def sigs(self):
        from .msigdb import sigs
        e = self.var_ensembl[['entrez']].to_dataframe().reset_index()
        e['e'] = 1
        e = e.drop_duplicates()
        e['entrez'] = e.entrez.astype(int)
        e = e.set_index(['var', 'entrez']).e
        e = xa.DataArray.from_series(e, sparse=True)
        e = e.fillna(0)
        e.data = sparse.GCXS(e.data, compressed_axes=[1])

        s = sigs.db.db.copy().astype(int)
        s.data = sparse.GCXS(s.data, compressed_axes=[1])
        s = xa.apply_ufunc(
            np.matmul, e, s,
            input_core_dims=[['var', 'entrez'], ['entrez', 'sig']],
            output_core_dims=[['var', 'sig']],
            join='inner'
        )
        s.data = sparse.GCXS(s.data, compressed_axes=[0])
        s = (s>0).astype(np.float64)
        return s
    _dataset1.sigs = sigs
_()

dataset1 = _dataset1()

#%%
if __name__=='__main__':
    self = dataset1
    import sparse
    x1 = self.X.map_blocks(
        lambda x: x.sum(keepdims=True), 
        meta=sparse.zeros((0,0), dtype=self.X.dtype)
    )
    x1 = x1.compute()

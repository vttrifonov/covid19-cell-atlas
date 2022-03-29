
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
        config.cache/'data'/'nih-adaptive.h5ad',
        backed='r'
    )
    return d

def _data2():
    d = scanpy.read(
        config.cache/'data'/'nih-adaptive.h5ad'
    )
    return d

class _nih_adaptive:
    pass

def _():
    _nih_adaptive.storage = config.cache/'nih-adaptive'

    @compose(property, lazy)
    def data(self):
        return _data1()
    _nih_adaptive.data = data

    @compose(property, lazy)
    def X1(self):
        storage = self.storage/'X1.npz'
        if storage.exists():
            x = sparse.load_npz(storage)
            return x
                    
        x = _data2().X
        x = sparse.COO.from_scipy_sparse(x)
        storage.parent.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(storage, x)
        return x
    _nih_adaptive.X1 = X1

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
    _nih_adaptive.var_ensembl = var_ensembl
_()

nih_adaptive = _nih_adaptive()

#%%
if __name__=='__main__':
    self = nih_adaptive

#%%
    import sparse
    x1 = self.X.map_blocks(
        lambda x: x.sum(keepdims=True), 
        meta=sparse.zeros((0,0), dtype=self.X.dtype)
    )
    x1 = x1.compute()

# %%

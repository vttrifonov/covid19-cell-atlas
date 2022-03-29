
#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import scanpy
import xarray as xa
import numpy as np
import sparse
from .helpers import config
from .common.caching import lazy, compose, XArrayCache, SparseCache

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

    @compose(property, lazy, SparseCache())
    def X1(self):
        x = _data2().X
        x = sparse.COO.from_scipy_sparse(x)
        return x
    _nih_innate.X1 = X1

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

# %%

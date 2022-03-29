
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

    @compose(property, lazy, SparseCache())
    def X(self):
        x = _data2().X
        x = sparse.COO.from_scipy_sparse(x)
        return x
    _nih_adaptive.X = X

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

    @compose(property, lazy, XArrayCache())
    def obs(self):
        o = self.data.obs.copy()

        o['sample'] = o.donor.astype(str) + '_' + o.timepoint.astype(str)
        
        _ = o.dsm_severity_score.astype(str)
        _ = _.where(_!='', 'nan')
        o['dsm_severity_score'] = _.astype(np.float32)

        samples = o.drop(columns=[
            'cell_type_ontology_term_id', 
            'cell_type'
        ]).drop_duplicates()
        donors = samples.drop(columns=[
            'sample', 'timepoint',
            'days_since_onset', 'days_since_hospitalized'
        ]).drop_duplicates().set_index('donor').to_xarray()
        samples = samples[[
            'sample', 'donor', 'timepoint',
            'days_since_onset', 'days_since_hospitalized'
        ]].rename(columns={'donor': 'DonorID'}).\
            set_index('sample').to_xarray()
        cells = o[['sample']].\
            rename(columns={'sample': 'SampleID'}).\
            to_xarray().rename(index='cell')

        r = xa.merge([samples, donors, cells])

        return r
    _nih_adaptive.obs = obs
_()

nih_adaptive = _nih_adaptive()

#%%
if __name__=='__main__':
    self = nih_adaptive

# %%

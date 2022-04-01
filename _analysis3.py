#%%

if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xa
from pathlib import Path
import sparse

from .common.caching import compose, lazy, XArrayCache
from ._helpers import config, xa_mmult
from .covid19_time_resolved_paper import data as paper
from ._analysis1 import analysis1
from ._analysis2 import analysis2

#%%
class _analysis3:
    storage = Path(config.cache)/'analysis3'

    @compose(property, lazy)
    def cytokines(self):
        x1 = paper.cytokines
        x1 = x1.to_dataframe().reset_index(drop=True)
        x1['sample'] = x1['DonorID'] + '_' + x1['days_from_symptom_onset_to_test'].astype(str)
        x1 = xa.merge([
           x1[['level', 'cytokine', 'sample']].set_index(['cytokine', 'sample']).to_xarray(),
           x1[['sample', 'DonorID', 'days_from_symptom_onset_to_test']].\
               drop_duplicates().\
               set_index('sample').to_xarray(),
        ])

        x2 = self.pseudobulk.drop_dims(['cell_type', 'gene'])
        x2 = x2.drop([
            'tot',
            'days_since_onset', 'days_since_hospitalized',
            'timepoint', 'subset',
            'assay', 'assay_ontology_term_id', 'material_type',
            'assay_ontology_term_id', 'tissue', 
            'tissue_ontology_term_id',
        ])
        x2 = x2.to_dataframe().reset_index(drop=True).drop_duplicates()
        x2 = x2.merge(
            x1[['DonorID']].to_dataframe().reset_index().\
                drop_duplicates(),
            left_on='donor', right_on='DonorID'
        ).drop(columns=['donor']).\
            set_index('sample').to_xarray()
        
        x1 = xa.merge([x1, x2], join='inner')

        x1 = x1.rename(
            days_from_symptom_onset_to_test='days_since_onset',
            DonorID='donor'
        )

        return x1
        
    @compose(property, lazy)
    def pseudobulk(self):
        o1 = analysis1.obs.copy()
        o1['sample'] = o1['sample']+'_innate'
        o1['subset'] = ('sample', ['innate']*o1.sizes['sample'])

        o2 = analysis2.obs.copy()
        o2['sample'] = o2['sample']+'_adaptive'
        o2['subset'] = ('sample', ['adaptive']*o2.sizes['sample'])

        o = xa.merge([o1, o2])

        x2 = o.drop_dims(['sample', 'cell'])
        x2 = x2.rename(donor='DonorID')
        x3 = np.where(
            x2.severity=='', 'Control',
            x2.severity+'-'+x2.outcome
        )
        x2['status'] = ('DonorID', x3)
        x2 = x2.to_dataframe().reset_index()

        x3 = o.drop_dims(['donor', 'cell'])
        _ = x3.days_since_onset
        _ = _.where(_!='', 'nan').astype(np.float32)
        x3['days_since_onset'] = _
        x3 = x3.to_dataframe().reset_index()
        
        x4 = x3.merge(x2, on='DonorID')
        x4 = x4.set_index('sample').to_xarray()
        x4 = x4.rename(DonorID='donor')

        x1 = analysis1.X2.copy()
        x2 = analysis2.X2.copy()

        x1['sample'] = x1.sample + '_innate'
        x2['sample'] = x2.sample + '_adaptive'
        x = xa.concat([x1, x2], dim='sample')

        x['tot'] = x.X.sum(dim=['cell_type', 'gene'])
        x['X'] = 1e6*x.X/x.tot

        r = xa.merge([x, x4], join='inner')

        return r
    
    @compose(property, lazy, XArrayCache())
    def fit1(self):
        from dask import delayed, compute
        from dask.distributed import LocalCluster, Client
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm
        from .sigs.fit import multipletests

        x1 = self.pseudobulk
        x1 = x1[['X', 'days_since_onset', 'dsm_severity_score_group', 'subset']]
        x1 = x1.sum(dim='cell_type')
        x1['X'] = np.log1p(x1.X)/np.log(2)

        def f1(x):
            x3 = smf.ols('X~days_since_onset', data=x).fit()
            x4 = smf.ols('X~days_since_onset * dsm_severity_score_group', data=x).fit()
            x5 = anova_lm(x3, x4).loc[[1],:].reset_index(drop=True)
            x5 = pd.concat([x5, pd.DataFrame(x4.params).T], axis=1)
            return x5
        
        def f2(x):
            x6 = x1.sel(gene=x).to_dataframe().reset_index()
            x6 = x6[x6.dsm_severity_score_group!='']
            x6 = x6.groupby(['subset', 'gene'])
            x6 = x6.apply(f1)
            return x6

        x7 = np.array_split(x1.gene.data, 28)
        with LocalCluster(n_workers=28) as cluster:
            with Client(cluster) as client:
                x8 = compute(delayed(f2)(g) for g in x7)
                x8 = pd.concat(x8[0])
                x8 = x8.reset_index()
                x8 = x8.drop(columns='level_2')

        x8['q'] = x8.groupby('subset')['Pr(>F)'].transform(multipletests, method='fdr_bh')
        x8 = x8.set_index(['subset', 'gene']).to_xarray()
        return x8

    @compose(property, lazy, XArrayCache())
    def fit2(self):
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm
        from .sigs.fit import multipletests

        x1 = self.cytokines
        x1 = x1[['level', 'days_since_onset', 'dsm_severity_score_group']]
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[x1.dsm_severity_score_group!='']

        def f1(x):
            x3 = smf.ols('level~days_since_onset', data=x).fit()
            x4 = smf.ols('level~days_since_onset * dsm_severity_score_group', data=x).fit()
            x5 = anova_lm(x3, x4).loc[[1],:].reset_index(drop=True)
            x5 = pd.concat([x5, pd.DataFrame(x4.params).T], axis=1)
            return x5
        
        x8 = x1.groupby(['cytokine']).apply(f1)
        x8 = x8.reset_index().drop(columns='level_1')
        x8['q'] = multipletests(x8['Pr(>F)'], method='fdr_bh')
        x8 = x8.set_index(['cytokine']).to_xarray()
        return x8

    @compose(property, XArrayCache())
    def symbol_entrez(self):
        from .sigs.entrez import symbol_entrez

        x = self.pseudobulk.gene.data
        x = symbol_entrez(x)
        return x

    @property
    def sigs(self):
        from .sigs._sigs import sigs
        x = sigs.all1
        x = x.sel(sig=x.sig_prefix.isin([
            'REACTOME', 'KEGG1', 
            'BIOCARTA', 'HALLMARK'
        ]))
        return x

analysis3 = _analysis3()

#%%
@compose(property, lazy)
def sigs1(self):
    g = self.symbol_entrez
    g = g.rename(Entrez_Gene_ID='gene')
    g = g/g.sum(dim='symbol').todense()
    g.data = g.data.tocoo()

    t = self.fit1.F.rename(gene='symbol')
    t = np.log(t).fillna(0)
    t.data = sparse.COO.from_numpy(t.data)
    t = xa_mmult(
        t, g,
        dim_x=['subset', 'symbol'], dim_y=['symbol', 'gene']
    )

    s = self.sigs

    r = xa.merge([t.rename('t'), s.rename('s')], join='inner')
    return r
_analysis3.sigs1 = sigs1

@compose(property, lazy, XArrayCache())
def enrich1(self):
    from .sigs.fit import enrichment
    x = self.sigs1
    x = enrichment(x.t, x.s)
    return x
_analysis3.enrich1 = enrich1

@compose(property, lazy, XArrayCache())
def enrich2(self):
    from .sigs.fit import fit_gsea
    x = self.sigs1
    x = fit_gsea(x.t, x.s, 1e5)
    return x
_analysis3.enrich2 = enrich2

#%%
if __name__ == '__main__':
    self = analysis3

# %%

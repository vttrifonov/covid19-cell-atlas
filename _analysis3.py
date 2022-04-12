#%%

if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from contextlib import AbstractContextManager
from pathlib import Path
from re import X
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

    @compose(property, lazy, XArrayCache())
    def enrich1(self):
        from .sigs.fit import enrichment
        x = self.sigs1
        x = enrichment(x.t, x.s)
        return x

    @compose(property, lazy, XArrayCache())
    def enrich2(self):
        from .sigs.fit import fit_gsea
        x = self.sigs1
        x = fit_gsea(x.t, x.s, 1e5)
        return x

analysis3 = _analysis3()

#%%
if __name__ == '__main__':
    self = analysis3

# %%
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    from .sigs.fit import multipletests
    from ._helpers import dataarray_from_series
    from sklearn.decomposition import PCA
    from plotnine import *
    from scipy import stats

    x1 = self.cytokines
    x1 = x1[['level', 'days_since_onset', 'dsm_severity_score', 'dsm_severity_score_group', 'status', 'donor']]
    #x1 = x1.sel(sample=x1.dsm_severity_score_group=='DSM_high')
    x1 = x1.sel(cytokine=~x1.cytokine.isin(['IFN-alpha2a']))
    x1['level'] = np.log1p(x1.level)/np.log(2)
    x1 = x1.to_dataframe().reset_index()
    x1 = x1[x1.dsm_severity_score_group!='']
    x1['time1'] = x1.days_since_onset

#%%
    def f1(x):
        x = x[~x['level'].isna()].copy()
        for d in range(2,4):
            x[f'level_{d}'] = x['level']**d
        x3 = smf.ols('time1~1', data=x).fit()
        x4 = smf.ols('time1~level', data=x).fit()
        x5 = anova_lm(x3, x4).loc[[1],:].reset_index(drop=True)
        x5 = pd.concat([x5, pd.DataFrame(x4.params).T], axis=1)
        x5['rsq'] = x5.ss_diff/(x5.ss_diff+x5.ssr)
        return x5
    
    x8 = x1.groupby(['cytokine']).apply(f1)
    x8 = x8.reset_index().drop(columns='level_1')
    x8['q'] = multipletests(x8['Pr(>F)'], method='fdr_bh')
    x8 = x8.set_index(['cytokine']).to_xarray()

#%%
    x9 = x8.to_dataframe().reset_index()
    x9 = x9.sort_values('Pr(>F)')
    x9 = x9[x9.q<0.05]
    x9

#%%
    for x11 in x9.itertuples():
        x10 = x1[x1.cytokine==x11.cytokine][['time1', 'level']]
        x10 = x10[~x10.level.isna()]
        x12 = pd.DataFrame(dict(
            level=np.linspace(x10.level.min(), x10.level.max(), num=100)
        ))
        x12['time1'] = x11.Intercept+x11.level*x12['level']#+x11.level_2*x12['level']**2+x11.level_3*x12['level']**3
        print(
            ggplot()+aes('level', 'time1')+
                geom_point(data=x10)+
                geom_line(data=x12)+
                labs(title=x11.cytokine)
        )

#%%
    x2 = x1[x1.cytokine.isin(x9.cytokine)].copy()
    x2 = x2.groupby(['donor', 'time1', 'cytokine', 'dsm_severity_score', 'dsm_severity_score_group', 'status'])['level'].mean()
    x2 = x2.reset_index()

    x5 = x2.copy()
    x5 = x5.pivot_table(index=['time1', 'donor', 'dsm_severity_score', 'dsm_severity_score_group', 'status'], columns='cytokine', values='level')
    x6 = "Q('"+x5.columns+"')"
    x6 = '+'.join(x6)
    x10 = x5.to_numpy()
    x18 = x5.columns
    x5 = x5.reset_index()
    x11 = np.isnan(x10).sum(axis=1)
    x5 = x5[x11==0]
    x10 = x10[x11==0,:]
    x10 = np.c_[x10]#, x10**2, x10**3]

    #x7 = smf.ols(f'time1~{x6}', data=x5).fit()
    #print(x7.summary())

    x12 = PCA(n_components=x10.shape[1]).fit(x10)
    x13 = x12.transform(x10)
    x13 = [
        np.c_[[1]*x13.shape[0], x13],
        np.r_[
            np.c_[[1], x12.mean_.reshape(1, -1)],
            np.c_[[0]*x12.n_components, x12.components_]
        ]
    ]
    #print(np.quantile(x13[0] @ x13[1] - sm.add_constant(x10), [0, 0.5, 1]))
    x13[1] = np.linalg.inv(x13[1].T @ x13[1]) @ x13[1].T
    
    x14 = [
        sm.OLS(x5.time1.to_numpy(), x13[0][:,:i]).fit()
        for i in range(2, x13[0].shape[1]+1)
    ]
    #print(x14[-1].summary())
    x14 = sorted(x14, key=lambda x: x.f_pvalue)
    print(x14[0].summary())

    x15 = x13[1][:,:len(x14[0].params)]@x14[0].params    

    print(pd.DataFrame(dict(
        coef=['const'] + x18.to_list(),
        value=x15
    )))

    x5['time2'] = sm.add_constant(x10)@x15
    print(
        ggplot(x5)+aes('time1', 'time2', color='status')+
            geom_point()+
            geom_abline(slope=1, intercept=0)
    )

    print(
        ggplot(x5)+aes('status', 'time2-time1')+
            geom_boxplot()+
            geom_point()
    )

    print(
        ggplot(x5)+aes('dsm_severity_score', 'time2-time1')+
            geom_point()+geom_smooth(method='lm')
    )

    print(np.sqrt(np.mean((x5.time1-x5.time2)**2)))

    x5['delta_time'] = x5.time2-x5.time1
    stats.ttest_ind(
        x5.delta_time[x5.dsm_severity_score_group=='DSM_high'],
        x5.delta_time[x5.dsm_severity_score_group=='DSM_low']
    )

# %%

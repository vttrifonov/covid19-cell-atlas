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
from ._helpers import config, xa_mmult, loess
from .covid19_time_resolved_paper import data as paper
from ._analysis1 import analysis1
from ._analysis2 import analysis2

#%%
def gmm_score1(x1, comp, covariance_type='full', n_perm = 0):
    from sklearn.mixture import GaussianMixture

    x = x1.data

    def gmm(comp):
        resp = np.zeros((x.shape[0], comp.max()+1))
        resp[np.arange(resp.shape[0]), comp] = 1.0
        gmm = GaussianMixture(
            n_components = resp.shape[1], 
            covariance_type = covariance_type
        )
        gmm._initialize(x, resp)
        return gmm

    gmm1 = gmm(comp)
    
    s, p = gmm1.score(x), np.nan
    if n_perm > 0:
        p = np.array([
            gmm(np.random.permutation(comp)).score(x)>s
            for _ in range(n_perm)
        ]).mean()

    x4 = xa.Dataset(coords={'donor': x1.donor.data, 't': x1.t.data})
    x4['means'] = (('clust', 't'), gmm1.means_)
    x4['covs'] = (('clust', 't', 't'), gmm1.covariances_)
    x4['weights'] = (('clust',), gmm1.weights_)
    x4['resp'] = (('donor',), comp)
    x4['score'] = s
    x4['p'] = p

    return x4

def gmm_score2(x2, n):
    from sklearn.mixture import GaussianMixture
    gmm = GaussianMixture(
        n_components = n,
        covariance_type = 'full'
    )
    gmm.fit(x2.data)
    clust = gmm.predict(x2.data)

    x5 = xa.Dataset(coords={'donor': x2.donor.data, 't': x2.t.data})
    x5['means'] = (('clust', 't'), gmm.means_)
    x5['covs'] = (('clust', 't', 't'), gmm.covariances_)
    x5['weights'] = (('clust',), gmm.weights_)
    x5['resp'] = (('donor',), clust)
    x5['score'] = gmm.score(x2.data)

    return x5

def contingency(x6):
    import statsmodels.api as sm
    x6 = sm.stats.Table.from_data(x6)
    return (
        x6.table_orig, np.round(x6.fittedvalues), x6.resid_pearson,
        x6.test_nominal_association().pvalue
    )

def hclust(x, l):
    from scipy.cluster.hierarchy import dendrogram, linkage

    x2 = linkage(x.data, 'average')
    x3 = dendrogram(x2, labels=None, color_threshold=l, no_plot=True)

    x4 = np.array(x3['ivl']).astype(np.int32)

    x5 = pd.DataFrame({
        'donor': x.donor.data[x4],
        'clust': x3['leaves_color_list'], 
        'ord': range(len(x4))
    }).set_index('donor').to_xarray()
    x5['linkage'] = (('node', 'link'), x2)
    x5['threshold'] = l

    return x5


#%%
class _analysis3:
    storage = Path(config.cache)/'analysis3'

    @compose(property, lazy)
    def donor(self):
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
        x2 = x2.set_index('donor').to_xarray()
        return x2

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

        x2 = self.donor.to_dataframe().reset_index()
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

    @compose(property, lazy, XArrayCache())
    def fit_time1(self):
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm
        from .sigs.fit import multipletests

        x1 = self.cytokines
        x1 = x1.sel(cytokine=~x1.cytokine.isin(['IFN-alpha2a']))
        x1 = x1[['level', 'days_since_onset', 'donor']]
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1['level'].isna()]

        def f1(x):
            x3 = smf.ols('days_since_onset~1', data=x).fit()
            x4 = smf.ols('days_since_onset~level', data=x).fit()
            x5 = anova_lm(x3, x4).loc[[1],:].reset_index(drop=True)
            x5 = pd.concat([x5, pd.DataFrame(x4.params).T], axis=1)
            x5['rsq'] = x5.ss_diff/(x5.ss_diff+x5.ssr)
            return x5
        
        x8 = x1.groupby(['cytokine']).apply(f1)
        x8 = x8.reset_index().drop(columns='level_1')
        x8['q'] = multipletests(x8['Pr(>F)'], method='fdr_bh')
        x8 = x8.set_index(['cytokine']).to_xarray()

        return x8

    @compose(property, lazy, XArrayCache())
    def fit_time2(self):
        import statsmodels.api as sm
        from sklearn.decomposition import PCA
        from statsmodels.stats.anova import anova_lm

        x9 = self.fit_time1
        x9 = x9.to_dataframe().reset_index()
        x9 = x9.sort_values('Pr(>F)')
        x9 = x9[x9.q<0.05]

        x1 = self.cytokines
        x1 = x1[['level', 'days_since_onset', 'donor']]
        x1 = x1.sel(cytokine=x9.cytokine.to_list())
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1['level'].isna()]

        x5 = x1.groupby(['donor', 'days_since_onset', 'cytokine'])['level'].mean()
        x5 = x5.reset_index()
        x5 = x5.pivot_table(index=['days_since_onset', 'donor'], columns='cytokine', values='level')
        x10 = x5.to_numpy()
        x18 = x5.columns
        x5 = x5.reset_index()
        x11 = np.isnan(x10).sum(axis=1)
        x5 = x5[x11==0]
        x10 = x10[x11==0,:]

        x12 = PCA(n_components=x10.shape[1]).fit(x10)
        x13 = x12.transform(x10)
        x13 = [
            np.c_[[1]*x13.shape[0], x13],
            np.r_[
                np.c_[[1], x12.mean_.reshape(1, -1)],
                np.c_[[0]*x12.n_components, x12.components_]
            ]
        ]
        x13[1] = np.linalg.pinv(x13[1])
        
        x14 = [
            sm.OLS(x5.days_since_onset.to_numpy(), x13[0][:,:i]).fit()
            for i in range(1, x13[0].shape[1]+1)
        ]
        x14_1 = x14[0]
        x14 = sorted(x14[1:], key=lambda x: x.f_pvalue)[0]

        x15 = x13[1][:,:len(x14.params)]@x14.params
        x15 = pd.Series(x15, index=['const'] + x18.to_list())
        x16 = anova_lm(x14_1, x14).loc[[1],:].reset_index(drop=True)
        x16['rsq'] = x16.ss_diff/(x16.ss_diff+x16.ssr)
        x16['sd'] = np.sqrt(x16.ssr/x16.df_resid)
        x16 = x16.T[0]

        x17 = xa.merge([
            x15.to_xarray().rename(index='coef').rename('coef_value'),
            x16.to_xarray().rename(index='stat').rename('stat_value')
        ])

        return x17

    @compose(property, lazy, XArrayCache()) 
    def fit_level1(self):
        import statsmodels.formula.api as smf
        from statsmodels.stats.anova import anova_lm
        from .sigs.fit import multipletests

        x1 = self.cytokines
        x1 = x1.sel(cytokine=~x1.cytokine.isin(['IFN-alpha2a', 'Ferritin']))
        x1 = x1[[
            'level', 'days_since_onset', 'donor',
            'dsm_severity_score', 'dsm_severity_score_group', 'status'
        ]]
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1['level'].isna()]
        x1['days_since_onset_2'] = x1['days_since_onset']**2

        x2 = x1[['donor', 'cytokine']].value_counts().reset_index()
        x2 = x2[x2[0]>=3].set_index(['donor', 'cytokine'])[0].to_xarray()
        x2 = x2.sel(cytokine=x2.isnull().sum(dim='donor')==0)
        x1 = x1[x1.donor.isin(x2.donor.data)]
        x1 = x1[x1.cytokine.isin(x2.cytokine.data)]

        def f1(x):
            x3 = smf.ols('level~1', data=x).fit()
            x4 = smf.ols('level~days_since_onset', data=x).fit()
            x5 = anova_lm(x3, x4).loc[[1],:].reset_index(drop=True)
            x5 = pd.concat([x5, pd.DataFrame(x4.params).T], axis=1)
            x5['rsq'] = x5.ss_diff/(x5.ss_diff+x5.ssr)
            return x5
        
        x8 = x1.groupby(['cytokine', 'donor']).apply(f1)
        x8 = x8.reset_index().drop(columns='level_2')
        x8['q'] = multipletests(x8['Pr(>F)'], method='fdr_bh')
        x8 = x8.set_index(['donor', 'cytokine']).to_xarray()

    @compose(property, lazy)
    def fit_level2_data(self):
        x1 = self.cytokines
        x1 = x1.sel(cytokine=~x1.cytokine.isin(['IFN-alpha2a', 'Ferritin']))
        x1 = x1[[
            'level', 'days_since_onset', 'donor',
            'dsm_severity_score', 'dsm_severity_score_group', 'status'
        ]]
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1['level'].isna()]
        nd = 1+2
        for d in range(nd):
            x1[f'days_since_onset_{d}'] = x1['days_since_onset']**d 
        f = [f'days_since_onset_{d}' for d in range(nd)]
        f = '+'.join(f)
        f = f'level~0+{f}'

        return x1, nd, f

    @compose(property, lazy, XArrayCache())
    def fit_level2(self):
        import statsmodels.formula.api as smf

        x1, _, f = self.fit_level2_data

        def f2(x):
            def f3(d):
                x2 = x.groupby(x.donor==d).size()
                x.loc[x.donor==d, 'weight'] = 10*1/x2[True]
                x.loc[x.donor!=d, 'weight'] = 1*1/x2[False]
                x5 = smf.wls(f, data=x, weights=x.weight).fit()
                x5 = pd.concat([
                    pd.DataFrame({'p': [x5.f_pvalue]}),
                    pd.DataFrame(x5.params).T
                ], axis=1)
                return x5
            x5 = [f3(d).assign(donor=d) for d in x.donor.drop_duplicates()]
            x5 = pd.concat(x5)
            return x5

        x8 = x1.groupby(['cytokine']).apply(f2)
        x8 = x8.reset_index().drop(columns='level_1')
        x8 = x8.set_index(['donor', 'cytokine']).to_xarray()
        return x8

    def _fit_level3(self, c, w):
        x, _, _ = self.fit_level2_data
        x = x[x.cytokine==c]
        x6 = np.quantile(x.days_since_onset, [0,1])
        x6 = np.linspace(*x6, 100)

        x1 = x.days_since_onset.to_numpy()
        y1 = x.level.to_numpy()
        g = x.donor.to_numpy()

        g3 = np.unique(g)
        x3 = np.zeros((len(g3), len(x6)))
        for i in range(len(g3)):
            x4 = g==g3[i]
            _, x2 = np.unique(x4, return_counts=True)
            x4 = np.where(x4, w*1/x2[1], 1/x2[0])
            x3[i,:] = loess(x1, y1, w=x4, t=x6)
        x3 = xa.DataArray(
            x3,
            coords=[('donor', g3), ('t', x6)]
        ).expand_dims(cytokine=[c])

        return x3

    def _fit_level4(self, c, w, sigma2, nt):
        from scipy.spatial.distance import pdist, squareform

        x, _, _ = self.fit_level2_data
        x = x[x.cytokine==c]
        x = x[x.days_since_onset<=40]
        x6 = np.quantile(x.days_since_onset, [0,1])
        x6 = np.linspace(*x6, 50)

        x1 = x.days_since_onset.to_numpy()
        y1 = x.level.to_numpy()
        g = x.donor.to_numpy()

        g3, g4 = np.unique(g, return_inverse=True)
        x3 = np.zeros((len(g3), len(x6)))
        x5 = np.ones((len(g3), len(g3)))
        for j in range(nt):
            for i in range(len(g3)):
                x7 = x5[i,g4]
                #x7 = x7/x7.sum()
                x4 = g==g3[i]
                x7[x4] = w*x7[x4]/x7[x4].sum()
                x7[~x4] = x7[~x4]/x7[~x4].sum()
                x8 = ~np.isnan(x7)
                x3[i,:] = loess(x1[x8], y1[x8], w=x7[x8], t=x6)
            x5 = squareform(pdist(x3)**2)
            x5 = np.exp(-x5/sigma2)

            #import matplotlib.pyplot as plt
            #p = plt.imshow(x5)
            #plt.show()
            
        x3 = xa.DataArray(
            x3,
            coords=[('donor', g3), ('t', x6)]
        ).expand_dims(cytokine=[c])

        return x3

    @compose(property, lazy, XArrayCache())
    def fit_level4(self):
        x, _, _ = self.fit_level2_data
        x = x.cytokine.drop_duplicates().to_list()
        def _(c):
            print(c)
            return self._fit_level4(c, 3, 100, 10).rename('level').to_dataset()
        x = [_(c) for c in x]
        x = xa.concat(x, dim='cytokine')
        return x

    @compose(property, lazy, XArrayCache())
    def fit_level4_pval1(self):
        def f1(x1):
            def f2():
                x2 = gmm_score2(x1.level, 2)
                x3 = xa.merge([
                    x2.resp, x1.dsm_severity_score_group.drop('cytokine')
                ]).to_dataframe()
                _, _, _, x3 = contingency(x3)
                return x3
            
            print(x1.cytokine.data)
            x3 = [f2() for _ in range(100)]
            x3 = np.exp(np.mean(np.log(x3)))
            return xa.DataArray(x3, name='p')

        x = self.fit_level4
        x = xa.merge([x, self.donor.dsm_severity_score_group], join='inner')
        x = x.sel(donor=x.dsm_severity_score_group!='')
        x = x.groupby('cytokine').apply(f1)
        return x
        

analysis3 = _analysis3()

#%%
if __name__ == '__main__':
    self = analysis3

#%%
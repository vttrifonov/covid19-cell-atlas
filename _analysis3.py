
#%%

if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xa
from pathlib import Path
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests
from dask import delayed, compute
from dask.distributed import LocalCluster, Client


from .common.caching import lazy, XArrayCache, compose, CSVCache
from .helpers import config
from .covid19_time_resolved_paper import data as paper
from ._analysis1 import analysis1
from ._analysis2 import analysis2

#%%

class _analysis3:
    storage = Path(config.cache)/'analysis3'

    @property
    def metadata(self):
        return paper.metadata
    
    @compose(property, lazy)
    def cytokines(self):
        return paper.cytokines
        
    @property
    def obs(self):
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
        return x4

    @property
    def X(self):
        x1 = analysis1.X2.copy()
        x2 = analysis2.X2.copy()

        x1['sample'] = x1.sample + '_innate'
        x2['sample'] = x2.sample + '_adaptive'
        x = xa.concat([x1, x2], dim='sample')

        x['tot'] = x.X.sum(dim=['cell_type', 'gene'])
        x['X'] = 1e6*x.X/x.tot
        return x

    def data2(self, gene):
        x2 = self.obs[['severity', 'outcome', 'dsm_severity_score', 'dsm_severity_score_group']].to_dataframe()
        x3 = np.where(
            x2.severity=='', 'Control',
            x2.severity+'-'+x2.outcome
        )
        x3 = pd.Categorical(
            x3,
            categories=[
                'Control', 'Moderate-alive', 'Severe-alive', 
                'Critical-alive', 'Critical-deceased'
            ]
        )        
        x2['status'] = x3

        x3 = self.obs[['days_since_onset', 'DonorID', 'timepoint', 'subset']]
        x3 = x3.to_dataframe().reset_index()
        x3 = x3.rename(columns={'DonorID': 'donor'})

        _ = x3.days_since_onset
        _ = _.where(_!='', 'nan').astype(np.float32)
        x3['days_since_onset'] = _
        
        x1 = self.X2.X
        x1 = 1e6*x1/x1.sum(dim=['cell_type', 'gene'])
        x1 = x1.sel(gene=gene)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1.X.isna()]
        x1 = x1.merge(x3, on='sample')
        x1 = x1.merge(x2, on='donor')
        return x1

    @property
    def app1_genes(self):
        x1 = self.cytokines
        return x1.cytokine.to_series().drop_duplicates().to_list()

    @property
    def app1_genes1(self):
        return self.X.gene.data

    def app1_plot1_data(self, gene):
        x1 = self.cytokines
        x1 = x1.sel(subject_test_day=x1.cytokine.isin([gene]))
        x1 = x1.to_dataframe().reset_index(drop=True)

        x2 = self.metadata.DSM_group.to_dataframe().reset_index()

        x3 = x1.merge(x2, left_on='DonorID', right_on='donor')
        del x3['DonorID']
        x3 = x3.rename(columns={'days_from_symptom_onset_to_test': 'days_since_onset'})
        x3 = x3[~x3.DSM_group.isna()]
        x3 = x3[x3.days_since_onset<=40]        
        return x3

    def app1_plot2_data(self, gene):
        x1 = xa.merge([self.X.sel(gene=gene), self.obs])
        x1 = x1.sel(sample=x1.dsm_severity_score_group!='')
        x1 = x1.sum(dim='cell_type')
        x1 = x1.to_dataframe().reset_index()
        x3 = pd.Categorical(
            x1.status,
            categories=[
                'Control', 'Moderate-alive', 'Severe-alive', 
                'Critical-alive', 'Critical-deceased'
            ]
        )
        x1['status'] = x3
        return x1


    def app1_plot3_data(self, gene):
        x1 = xa.merge([self.X.sel(gene=gene), self.obs])
        x1 = x1.sel(sample=x1.dsm_severity_score_group!='')
        x1 = x1.to_dataframe().reset_index()
        x3 = pd.Categorical(
            x1.status,
            categories=[
                'Control', 'Moderate-alive', 'Severe-alive', 
                'Critical-alive', 'Critical-deceased'
            ]
        )
        x1['status'] = x3        
        return x1

    @compose(property, lazy, XArrayCache())
    def fit1(self):
        x1 = xa.merge([self.X, self.obs])
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

        x9 = ~x8['Pr(>F)'].isna()
        x8.loc[x9, 'q'] = x8[x9].groupby('subset')['Pr(>F)'].\
            transform(lambda x: multipletests(x, method='fdr_bh')[1])
        x8 = x8.set_index(['subset', 'gene']).to_xarray()
        return x8

analysis3 = _analysis3()

#%%
if __name__ == '__main__':
    from plotnine import *
    self = analysis3

#%%
    x3 = self.app1_plot1_data(['IL-6'])
    print(        
        ggplot(x3)+aes('days_since_onset', 'np.log1p(level)/np.log(10)')+
            geom_point(aes(fill='DSM_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='DSM_group'), alpha=0.5, method='lm')+
            facet_grid('cytokine~.', scales='free')+
            theme(figure_size=(6, 4))+
            labs(y='log2(pg/mL)')
    )

#%%
    x2 = self.app1_plot2_data('IL6')
    print(        
        ggplot(x2)+aes('days_since_onset', 'np.log1p(X)/np.log(2)')+
            geom_point(aes(fill='dsm_severity_score_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='dsm_severity_score_group'), alpha=0.5)+
            facet_grid('subset+gene~.', scales='free')+
            theme(figure_size=(6, 4))+
            labs(y='log2RPM')
    )

#%%
    x2 = self.app1_plot3_data('IL18BP')
    print(        
        ggplot(x2)+aes('days_since_onset', 'np.log1p(X)/np.log(2)')+
            geom_point(aes(fill='dsm_severity_score_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='dsm_severity_score_group'), alpha=0.5)+
            facet_grid('cell_type+gene~.', scales='free_x')+
            theme(figure_size=(4, 30))+
            labs(y='log2RPM')
    )

#%%
    x2 = x1[['donor', 'status', 'dsm_severity_score', 'dsm_severity_score_group']].drop_duplicates()
    print(
        ggplot(x2)+aes('status', 'dsm_severity_score')+
            geom_violin(aes(fill='status'))+
            geom_boxplot(width=0.05)+
            geom_point(aes(color='dsm_severity_score_group'))+
            theme(axis_text_x=element_text(angle=45, ha='right'))
    )

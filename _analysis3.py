
#%%

if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xa
from pathlib import Path

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
        return o

    @property
    def X2(self):
        x1 = analysis1.X2.copy()
        x2 = analysis2.X2.copy()

        x1['sample'] = x1.sample + '_innate'
        x2['sample'] = x2.sample + '_adaptive'
        x = xa.concat([x1, x2], dim='sample')
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
        return self.X2.X.gene.data

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
        x1 = self.data2(gene)
        x2 = x1[[
            'donor', 'days_since_onset', 'status', 
            'dsm_severity_score_group',
            'gene', 'subset',         
            'X'
        ]]
        x2 = x2[x2.dsm_severity_score_group!='']
        x2 = x2.groupby(list(set(x2.columns)-set(['X'])), observed=True)
        x2 = x2.X.sum().reset_index()
        return x2

    def app1_plot3_data(self, gene):
        x1 = self.data2(gene)
        x1 = x1[x1.dsm_severity_score_group!='']
        return x1

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
            geom_line(aes(group='DonorID'), alpha=0.1)+
            geom_smooth(aes(color='DSM_group'), alpha=0.5)+
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

#%%
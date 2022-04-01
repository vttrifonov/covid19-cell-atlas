#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import pandas as pd
import numpy as np

from ._analysis3 import analysis3
from ._helpers import round_float
from .common.caching import compose, lazy

class _app1:
    analysis = analysis3

    @compose(property, lazy)
    def genes(self):
        return self.analysis.pseudobulk.gene.data

    def plot2_data(self, gene):
        x1 = self.analysis.pseudobulk.sel(gene=gene)
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

    def plot3_data(self, gene):
        x1 = self.analysis.pseudobulk.sel(gene=gene)
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

    @compose(property, lazy)
    def genes_table(self):
        x = self.analysis.fit1.to_dataframe().reset_index()
        x = x[~x['Pr(>F)'].isna()]
        x = x.sort_values('q')
        x = x[[
            'subset', 'gene',
            'Pr(>F)', 'q', 'F',
            'Intercept',
            'days_since_onset',
            'dsm_severity_score_group[T.DSM_low]',
            'days_since_onset:dsm_severity_score_group[T.DSM_low]'
        ]]

        for c in ['Pr(>F)', 'q']:            
            x[c] = round_float(x[c], 1)

        for c in x:
            if c in ['q', 'Pr(>F)', 'subset', 'gene']:
                continue
            x[c] = round(x[c], 2)
        return x

    @compose(property, lazy)
    def enrich1_table(self):
        x = self.analysis.enrich1
        x = x.to_dataframe().reset_index()
        x = x[x.coef>0]
        x = x.sort_values('p')

        for c in ['p']:
            x[c] = round_float(x[c], 1)

        for c in ['coef', 'r2', 'se']:
            x[c] = round(x[c], 2)

        return x

app1 = _app1()

#%%
if __name__ == '__main__':
    from plotnine import *
    self = app1

#%%
    x2 = self.plot2_data('IL6')
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
    x2 = self.plot3_data('IL6')
    print(        
        ggplot(x2)+aes('days_since_onset', 'np.log1p(X)/np.log(2)')+
            geom_point(aes(fill='dsm_severity_score_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='dsm_severity_score_group'), alpha=0.5)+
            facet_grid('subset+cell_type+gene~.', scales='free_x')+
            theme(figure_size=(4, 30))+
            labs(y='log2RPM')
    )

#%%
    x2 = self.plot2_data('IL6')
    print(
        ggplot(x2)+aes('status', 'dsm_severity_score')+
            geom_violin(aes(fill='status'))+
            geom_boxplot(width=0.05)+
            geom_point(aes(color='dsm_severity_score_group'))+
            theme(axis_text_x=element_text(angle=45, ha='right'))
    )
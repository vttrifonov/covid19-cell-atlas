#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import pandas as pd
import numpy as np

from ._analysis3 import analysis3
from ._helpers import round_float

class _app2:
    analysis = analysis3

    @property
    def cytokines(self):
        return self.analysis.cytokines.cytokine.data

    def plot1_data(self, cytokine):
        x1 = self.analysis.cytokines.sel(cytokine=cytokine)
        x1 = x1.sel(sample=x1.dsm_severity_score_group!='')
        x1 = x1.to_dataframe().reset_index()
        return x1

    @property
    def cytokines_table(self):
        x = self.analysis.fit2.to_dataframe().reset_index()
        x = x[~x['Pr(>F)'].isna()]
        x = x.sort_values('q')
        x = x[[
            'cytokine',
            'Pr(>F)', 'q', 'F',
            'Intercept',
            'days_since_onset',
            'dsm_severity_score_group[T.DSM_low]',
            'days_since_onset:dsm_severity_score_group[T.DSM_low]'
        ]]
        x['Pr(>F)'] = round_float(x['Pr(>F)'], 1)
        x['q'] = round_float(x.q, 1)
        for c in x:
            if c in ['q', 'Pr(>F)', 'subset', 'cytokine']:
                continue
            x[c] = round(x[c], 2)
        return x

app2 = _app2()

#%%
if __name__ == '__main__':
    from plotnine import *
    self = app2

#%%
    x3 = self.plot1_data(['IL-6'])
    x3 = x3[x3.days_since_onset<=40]
    print(        
        ggplot(x3)+aes('days_since_onset', 'np.log1p(level)/np.log(2)')+
            geom_point(aes(fill='dsm_severity_score_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='dsm_severity_score_group'), alpha=0.5, method='lm')+
            facet_grid('cytokine~.', scales='free')+
            theme(figure_size=(6, 4))+
            labs(y='log2(pg/mL)')
    )

#%%

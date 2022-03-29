
#%%

if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xa
from pathlib import Path

from . import nih_innate
from .common.caching import lazy, XArrayCache, compose, CSVCache
from .helpers import config
from .covid19_time_resolved_paper import data as paper

#%%

class _analysis1:
    pass

def _():
    _analysis1.storage = Path(config.cache)/'analysis1'
    _analysis1.dataset = nih_innate

    @property
    def metadata(self):
        return paper.metadata
    _analysis1.metadata = metadata
    
    @compose(property, lazy)
    def cytokines(self):
        return paper.cytokines
    _analysis1.cytokines = cytokines
        
    @compose(property, lazy, XArrayCache())
    def obs(self):
        o = self.dataset.data.obs
        o['sample'] = o.donor.astype(str) + '_' + o.timepoint.astype(str)
        o = o.to_xarray().rename(index='cell')        
        return o
    _analysis1.obs = obs

    @compose(property, lazy, XArrayCache())
    def X1(self):
        import dask 

        dataset = self.dataset

        cs = dataset.data.var.copy()
        cs['col'] = range(cs.shape[0])
        cs = cs.col

        rs = pd.Series(range(dataset.data.obs.shape[0]), index=dataset.data.obs.index)

        def _(cell_type, x):
            print(cell_type)
            rs1 = rs[x.cell.data]

            with dask.config.set(**{'array.slicing.split_large_chunks': False}):
                mat = dataset.X[rs1.to_list(),:].compute()
            mat = xa.DataArray(mat, [('cell', rs1.index.to_list()), ('gene', cs.index.to_list())])
            mat = xa.merge([mat.rename('X'), x.sample])
            def _(sample, x):
                return xa.merge([
                    x.X.todense().sum(dim='cell').rename('X'),
                    xa.DataArray(x.sizes['cell'], name='n')
                ]).expand_dims({'sample': [sample]})
            _ = [dask.delayed(_)(*x) for x in mat.groupby('sample')]
            _ = dask.compute(*_)
            mat = xa.concat(_, dim='sample')
            return mat.expand_dims({'cell_type': [cell_type]})
        _ = [_(*x) for x in self.obs.groupby('cell_type')]
        mat = xa.concat(_, dim='cell_type')
        return mat
    _analysis1.X1 = X1

    @compose(property, lazy)
    def data1(self):
        x2 = self.obs[['donor', 'timepoint']].to_dataframe()
        x2 = x2.drop_duplicates()
        x2['sample'] = x2.donor + '_' + x2.timepoint
        x2 = x2.set_index('sample').to_xarray()

        x1 = xa.merge([self.X1.drop('n'), x2])
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1.X.isna()].copy()
        x1['X'] = 1e6*x1.X/x1.groupby(['sample']).X.transform('sum')
        return x1
    _analysis1.data1 = data1

    def data2(self, gene):
        x2 = self.obs[['severity', 'outcome', 'sample']].to_dataframe()
        x2 = x2.drop_duplicates().set_index('sample')
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
        x2 = pd.Series(x3, index=x2.index, name='status')
        x2 = x2.reset_index()

        x3 = self.metadata[['DonorID', 'days_from_symptom_onset_to_sample_drawn', 'visit']]
        x3 = x3.rename(visit='Time', days_from_symptom_onset_to_sample_drawn='days')
        x3 = x3.to_dataframe().reset_index(drop=True)

        x4 = self.metadata[['DSM']]
        x4['DSM_group'] = (
            'donor', 
            np.where(
                np.isnan(x4.DSM), 'nan',
                np.where(x4.DSM > np.nanmedian(x4.DSM), 'DSM_high', 'DSM_low')
            )
        )
        x4 = x4.rename(donor='DonorID')
        x4 = x4.to_dataframe().reset_index()
        
        x1 = self.data1
        x1 = x1[x1.gene==gene].copy()
        x1 = x1.merge(x2, on='sample')
        x1 = x1.merge(x3, left_on=['donor', 'timepoint'], right_on=['DonorID', 'Time'], how='left')
        x1 = x1.merge(x4, left_on=['donor'], right_on=['DonorID'], how='left')
        x1['days'] = x1.days.where(x1.timepoint!='HC', 0)
        return x1
    _analysis1.data2 = data2
_()

analysis1 = _analysis1()

#%%
if __name__ == '__main__':
    from plotnine import *
    self = analysis1

#%%
    x1 = self.cytokines
    x1 = x1.sel(subject_test_day=x1.cytokine.isin(['IL-1ra/IL-1F3']))
    x1 = x1.to_dataframe().reset_index(drop=True)

    x2 = self.metadata.DSM_group.to_dataframe().reset_index()

    x3 = x1.merge(x2, left_on='DonorID', right_on='donor')
    x3 = x3.rename(columns={'days_from_symptom_onset_to_test': 'days'})
    #x3['DSM_group'] = x3.DSM_group.where(~x3.DSM_group.isna(), 'nan')
    x3 = x3[~x3.DSM_group.isna()]
    x3 = x3[x3.days<=40]

    print(        
        ggplot(x3)+aes('days', 'np.log1p(level)/np.log(10)')+
            geom_point(aes(fill='DSM_group'), alpha=0.5)+
            geom_line(aes(group='DonorID'), alpha=0.1)+
            geom_smooth(aes(color='DSM_group'), alpha=0.5)+
            facet_grid('cytokine~.', scales='free')+
            theme(figure_size=(6, 4))+
            labs(y='log2(pg/mL)')
    )

#%%
    x1 = self.data2('IL6')  

#%%
    x2 = x1[['donor', 'status', 'DSM', 'DSM_group']].drop_duplicates()
    print(
        ggplot(x2)+aes('status', 'DSM')+
            geom_violin(aes(fill='status'))+
            geom_boxplot(width=0.05)+
            geom_point(aes(color='DSM_group'))+
            theme(axis_text_x=element_text(angle=45, ha='right'))
    )

#%%
    x2 = x1.groupby(['donor', 'gene', 'timepoint', 'days', 'status', 'DSM_group'], observed=True).sum().reset_index()
    x2 = x2[x2.DSM_group!='nan']
    print(        
        ggplot(x2)+aes('days', 'np.log1p(X)/np.log(2)')+
            geom_point(aes(fill='DSM_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='DSM_group'), alpha=0.5)+
            facet_grid('gene~.', scales='free')+
            theme(figure_size=(6, 4))+
            labs(y='log2RPM')
    )

#%%
    x2 = x1[x1.DSM_group!='nan']
    print(        
        ggplot(x2)+aes('days', 'np.log1p(X)/np.log(2)')+
            geom_point(aes(fill='DSM_group'), alpha=0.5)+
            geom_line(aes(group='donor'), alpha=0.1)+
            geom_smooth(aes(color='DSM_group'), alpha=0.5)+
            facet_grid('cell_type+gene~.', scales='free_x')+
            theme(figure_size=(4, 18))+
            labs(y='log2RPM')
    )

#%%
#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

import pandas as pd
import numpy as np
import xarray as xa

from ._analysis3 import analysis3
from .common.caching import compose, lazy
from ._helpers import round_float, loess

class _app3:
    analysis = analysis3

    @compose(property, lazy)
    def fit_level2(self):
        x1 = self.analysis.cytokines[[
            'donor', 'status', 'dsm_severity_score', 'dsm_severity_score_group'
        ]]
        x1 = x1.to_dataframe().drop_duplicates().reset_index(drop=True)

        x2, nd, _ = self.analysis.fit_level2_data
        x2 = x2.groupby(['cytokine', 'donor']).size().rename('n')
        x2 = x2.reset_index()

        x = self.analysis.fit_level2
        x = x.to_dataframe().reset_index()

        x = x2.merge(x, on=['cytokine', 'donor'])
        x = x1.merge(x, on='donor')

        for c in ['dsm_severity_score'] + [f'days_since_onset_{d}' for d in range(nd)]:
            x[c] = np.round(x[c], 2)
        x['p'] = round_float(x['p'], 2)

        return x

    def plot1_data(self, c, d):
        import statsmodels.formula.api as smf

        x, nd, f = self.analysis.fit_level2_data

        x = x[x.cytokine==c].copy()
        x2 = x.groupby(x.donor==d).size()
        x.loc[x.donor==d, 'weight'] = 10*1/x2[True]
        x.loc[x.donor!=d, 'weight'] = 1*1/x2[False]

        x4 = smf.wls(f, data=x, weights=x.weight).fit()
        x5 = smf.ols(f, data=x[x.donor==d]).fit()
        x7 = smf.ols(f, data=x).fit()

        def f1(p, t):
            p = pd.DataFrame(p).T
            p = np.array([p[f'days_since_onset_{d}'] for d in range(nd)]).T
            t = np.vstack([(t**d) for d in range(nd)])
            return p@t

        x6 = np.quantile(x.days_since_onset, [0,1])
        x6 = np.linspace(*x6, 100)
        x6 = pd.DataFrame({
            'days_since_onset': x6,
            'pred1': f1(x4.params, x6)[0,:],
            'pred2': f1(x5.params, x6)[0,:],
            'pred3': f1(x7.params, x6)[0,:]
        })
        x6 = x6.melt(id_vars='days_since_onset')
        x6 = x6.rename(columns={'value': 'level'})

        x7 = np.quantile(x['level'], [0,1])
        x6 = x6[x6.level>=x7[0]]
        x6 = x6[x6.level<=x7[1]]
        return x, x6

    @compose(property, lazy)
    def data2(self):
        _, nd, _  = self.analysis.fit_level2_data
        x8 = self.analysis.fit_level2        
        def _(t):
            x = [(x8[f'days_since_onset_{d}']*(t**d)).expand_dims(d=[d]) for d in range(nd)]
            x = xa.concat(x, dim='d')
            return x.sum(dim='d')
        x9 = [_(t).expand_dims(t=[t]) for t in np.linspace(0, 28, num=29)]
        x9 = xa.concat(x9, dim='t')
        return x9

    @compose(property, lazy)
    def data2_pca(self):    
        from sklearn.decomposition import PCA

        x9 = self.data2
        pca = PCA(n_components=2)
        x10 = xa.apply_ufunc(
            lambda x: pca.fit_transform(x), x9,
            input_core_dims=[['donor', 'cytokine']],
            output_core_dims=[['donor', 'pc']],
            vectorize=True
        )
        x10 = x10.assign_coords(pc=['pc1', 'pc2'])
        x10 = x10.rename('value')

        for i in range(1, x10.sizes['t']):
            for j in range(2):
                x12 = x10[i-1,:,j]
                x13 = x10[i,:,j]
                x14 = [
                    ((x12+x13)**2).sum().data,
                    ((x12-x13)**2).sum().data
                ]
                x14 = 2*np.argmin(x14)-1
                x10[i,:,j] = x14*x13
    
        x1 = self.analysis.cytokines[['donor', 'dsm_severity_score', 'dsm_severity_score_group', 'status']]
        x1 = x1.to_dataframe().drop_duplicates()
        x2 = x10.to_dataframe().reset_index()
        x2 = x2.pivot_table(index=['t', 'donor'], columns='pc', values='value').reset_index()
        x2 = x2.merge(x1, on='donor')
        x2 = x2[x2.dsm_severity_score_group!='']
        return x2

    @compose(property, lazy)
    def plot3_cytokines(self):
        return self.fit_level2.cytokine.drop_duplicates().to_list()


    def plot3_data(self, c):
        x1 = self.data2.rename('level')
        x1 = x1.sel(cytokine=c)
        x1 = x1.to_dataframe().reset_index()
        x2 = self.analysis.cytokines[['donor', 'dsm_severity_score_group']]
        x2 = x2.to_dataframe().drop_duplicates()
        x1 = x1.merge(x2, on='donor')
        x1 = x1[x1.dsm_severity_score_group!='']
        return x1

    def plot4_data(self, c, d, w=1):
        x, _, _ = self.analysis.fit_level2_data
        x = x[x.cytokine==c].copy()
        x2 = x.groupby(x.donor==d).size()
        x.loc[x.donor==d, 'weight'] = w*1/x2[True]
        x.loc[x.donor!=d, 'weight'] = 1*1/x2[False]

        #x['nw'] = np.round(x.weight*(1/x.weight).max())
        #x4 = np.concatenate([np.repeat(*x) for x in enumerate(x.nw)])
        #x4 = x.iloc[x4,:]

        x4 = x.copy()
        x5 = x[x.donor==d].copy().assign(weight=1)
        x7 = x.copy().assign(weight=1)

        #def f1(d, t):
        #    import statsmodels.api as sm        
        #    return sm.nonparametric.lowess(d.level, d.days_since_onset, xvals=t)

        def f1(d, t):
            return loess(
                d.days_since_onset.to_numpy(),
                d.level.to_numpy(),
                d.weight.to_numpy(),
                t
            )

        x6 = np.quantile(x.days_since_onset, [0,1])
        x6 = np.linspace(*x6, 100)
        x6 = pd.DataFrame({
            'days_since_onset': x6,
            'pred1': f1(x4, x6),
            'pred2': f1(x5, x6),
            'pred3': f1(x7, x6)
        })
        x6 = x6.melt(id_vars='days_since_onset')
        x6 = x6.rename(columns={'value': 'level'})

        x7 = np.quantile(x['level'], [0,1])
        x6 = x6[x6.level>=x7[0]]
        x6 = x6[x6.level<=x7[1]]
        return x, x6

    def _fit_level3(self, c, w):
        x, _, _ = self.analysis.fit_level2_data
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

        x, _, _ = self.analysis.fit_level2_data
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

            import matplotlib.pyplot as plt
            p = plt.imshow(x5)
            plt.show()
            
        x3 = xa.DataArray(
            x3,
            coords=[('donor', g3), ('t', x6)]
        ).expand_dims(cytokine=[c])

        return x3
    
    def plot4_data1(self, x3, d):
        def loess2(x1, y1, t):
            x1 = x1.to_numpy()
            y1 = y1.to_numpy()
            return loess(x1, y1, t=t)

        x1, _, _ = self.analysis.fit_level2_data
        x1 = x1[x1.cytokine==x3.cytokine.data[0]].copy()
        x2 = x1[x1.donor==d]
        x3 = x3.sel(donor=d)            

        x6 = x3.t.data
        x6 = pd.DataFrame({
            'days_since_onset': x6,
            'pred1': x3.data[0,:],
            'pred2': loess2(x2.days_since_onset, x2.level, x6),
            'pred3': loess2(x1.days_since_onset, x1.level, x6)
        })
        x6 = x6.melt(id_vars='days_since_onset')
        x6 = x6.rename(columns={'value': 'level'})

        x7 = np.quantile(x1['level'], [0,1])
        x6 = x6[x6.level>=x7[0]]
        x6 = x6[x6.level<=x7[1]]
        return x1, x6


app3 = _app3()

#%%
if __name__ == '__main__':
    from plotnine import *
    from plotnine.animation import PlotnineAnimation
    self = app3

#%%
self.analysis.fit2.to_dataframe().reset_index().sort_values('Pr(>F)').head(100)

#%%
    from scipy.spatial.distance import pdist, squareform
    from scipy.cluster.hierarchy import dendrogram, linkage
    import matplotlib.pyplot as plt

    x1 = self._fit_level4('ST2/IL-33R', 3, 100, 10).rename('level')
    x2 = linkage(x1.data[0], 'average')

    x4 = x1.donor.to_dataframe().reset_index(drop=True).reset_index()
    x3 = self.analysis.cytokines[['donor', 'dsm_severity_score_group']]
    x3 = x3.to_dataframe().drop_duplicates()
    x4 = x4.merge(x3, on='donor')    
    x4['label'] = x4['index'].astype(str)+' '+x4.dsm_severity_score_group

#%%
    import statsmodels.api as sm

    plt.figure(figsize=(10, 4))
    p = dendrogram(x2, labels=None, color_threshold=7.5)

    x4.loc[np.array(p['ivl']).astype(np.int32), 'ord'] = range(x4.shape[0])
    x4['clust'] = pd.Series(p['leaves_color_list'])[x4.ord].to_list()

    x5 = np.array(p['ivl']).astype(np.int32)
    x5 = x1.data[0][x4.sort_values('ord')['index'],:]
    #x5 = np.clip(x5, 0, 10)
    plt.figure(figsize=(10, 15))
    q = plt.imshow(x5)
    plt.show()

    x6 = x4[['clust', 'dsm_severity_score_group']]
    x6 = x6[x6.dsm_severity_score_group!='']
    x6 = sm.stats.Table.from_data(x6)
    
    print(x6.table_orig)
    print(np.round(x6.fittedvalues))
    print(x6.resid_pearson)
    print(x6.chi2_contribs)
    print(x6.test_nominal_association().pvalue)

#%%
    x7 = x1.to_dataframe().reset_index()
    x7 = x7.merge(x4, on='donor')
    x7 = x7[x7.dsm_severity_score_group!='']
    print(
        ggplot(x7)+aes('t', 'level', color='dsm_severity_score_group')+
            geom_line(aes(group='donor'))+
            facet_grid('clust~.', scales='free_y')
    )


#%%   
    for d in x4[x4.clust=='C2'].donor:
        x5, x6 = self.plot4_data1(x1, d)
        print(
            ggplot(x5)+
                aes('days_since_onset', 'level')+
                geom_point(aes(color=f'donor=="{d}"'))+
                geom_line(aes(group='donor', color=f'donor=="{d}"'))+
                geom_line(data=x6[x6.variable=='pred1'], linetype='dotted')+
                geom_line(data=x6[x6.variable=='pred2'], linetype='dashed')+
                geom_line(data=x6[x6.variable=='pred3'], linetype='solid')+
                theme(legend_position='none')+
                scale_color_manual(['lightgray', 'red'])
        )

#%%
    
    
#%%
    self.fit_level2

#%%
    x = self.fit_level2.iloc[0]
    x1, x2 = self.plot1_data(x.cytokine, x.donor)
    print(
        ggplot(x1)+
            aes('days_since_onset', 'level')+
            geom_point()+
            geom_line(aes(group='donor', color=f'donor=="{x.donor}"'))+
            geom_line(data=x2[x2.variable=='pred1'], linetype='dotted')+
            geom_line(data=x2[x2.variable=='pred2'], linetype='dashed')+
            geom_line(data=x2[x2.variable=='pred3'], linetype='solid')+
            theme(legend_position='none')
    )

#%%
    x2 = self.data2_pca
    x4 = np.quantile(x2.pc1, [0,1])
    x5 = np.quantile(x2.pc2, [0,1])
    x3 = [
        ggplot(x)+
            aes(
                'pc1', 'pc2',
                color='dsm_severity_score_group'
            )+
            geom_point()+
            lims(x=x4, y=x5)
        for t, x in x2.groupby('t')
    ]

    #from matplotlib import rc
    #rc('animation', html='html5')
    #%matplotlib ipympl

    p = PlotnineAnimation(
        x3, interval=28e3/len(x3), repeat_delay=500
    )
    p.save('xxx.mp4')
    #p

# %%
    x1 = self.plot3_data('IL-18/IL-1F4')
    print(
        ggplot(x1)+
            aes('t', 'level', color='dsm_severity_score_group')+
            geom_line(aes(group='donor'))
    )
# %%

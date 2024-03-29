#%%
if __name__=='__main__':
    __package__ = 'covid19_cell_atlas'

from re import X
import pandas as pd
import numpy as np
import xarray as xa

from ._analysis3 import analysis3, gmm_score1, gmm_score2, hclust, contingency
from .common.caching import compose, lazy
from ._helpers import round_float, loess


#%%
class _app4:
    analysis = analysis3

    @compose(property, lazy)
    def fit_level4(self):
        x1 = self.analysis.fit_level4
        x1 = xa.merge([x1, self.analysis.donor.dsm_severity_score_group], join='inner')
        x1 = x1.sel(donor=x1.dsm_severity_score_group!='')
        x1 = xa.merge([x1, self.analysis.fit_level4_pval1])
        return x1

    @compose(property, lazy)
    def cytokines(self):
        return self.fit_level4.p.to_series().sort_values().reset_index()

    def plot1_data(self, c):
        x1 = self.fit_level4
        x1 = x1.sel(cytokine=c)
        x1 = x1.squeeze().drop('cytokine')

        x5 = gmm_score2(x1.level, 2)
        x5 = xa.merge([x5, x1[['level', 'dsm_severity_score_group']]])
        x5 = x5.rename(dsm_severity_score_group='group')

        return [
            contingency(x5[['resp', 'group']].to_dataframe()),

            x5[['level', 'group', 'resp']].\
                to_dataframe().reset_index().\
                assign(
                    group = lambda x: 'obs_'+x.group
                ),

            x5[['means']].\
                rename(clust='group', means='level').\
                to_dataframe().reset_index().\
                assign(
                    resp = lambda x: x.group,
                    donor = lambda x: x.group.astype(str),
                    group = lambda x: 'mean_'+x.group.astype(str)
                ),

            c
        ]

    def plot2_data(self, c, d):
        def loess2(x1, y1, t):
            x1 = x1.to_numpy()
            y1 = y1.to_numpy()
            return loess(x1, y1, t=t)

        x1, _, _ = self.analysis.fit_level2_data
        x1 = x1[x1.cytokine==c].copy()
        x2 = x1[x1.donor==d]

        x3 = self.fit_level4
        x3 = x3.sel(cytokine=c)
        x3 = x3.sel(donor=d)

        x6 = x3.t.data
        x6 = pd.DataFrame({
            'days_since_onset': x6,
            'pred1': x3.level.data,
            'pred2': loess2(x2.days_since_onset, x2.level, x6),
            'pred3': loess2(x1.days_since_onset, x1.level, x6)
        })
        x6 = x6.melt(id_vars='days_since_onset')
        x6 = x6.rename(columns={'value': 'level'})

        x7 = np.quantile(x1['level'], [0,1])
        x6 = x6[x6.level>=x7[0]]
        x6 = x6[x6.level<=x7[1]]
        return x1, x6

    @compose(property, lazy)
    def fit_level4_deg(self):
        x = self.analysis.fit_level4_deg
        x = x.to_dataframe().reset_index()
        x = x.sort_values('q')
        return x

    def plot3_data(self, c, t):
        x = self.fit_level4
        x = x.sel(cytokine=c)
        x = x.sel(t=t)
        x = x.to_dataframe().reset_index()
        return x

    def fit_level4_pca(self, n = None):    
        from sklearn.decomposition import PCA
        
        if n is None:
            n = self.cytokines.shape[0]
        
        x1 = self.cytokines.iloc[:n].cytokine.to_list()

        x9 = self.fit_level4.level
        x9 = x9.sel(cytokine=x1)
        pca = PCA(n_components=2)
        def _(x):
            f = pca.fit(x)
            return f.transform(x), f.components_
        x10 = xa.apply_ufunc(
            _, x9,
            input_core_dims=[['donor', 'cytokine']],
            output_core_dims=[['donor', 'pc'], ['pc', 'cytokine']],
            vectorize=True
        )
        x10 = xa.merge([
            x10[0].assign_coords(pc=['pc1', 'pc2']).rename('trans'),
            x10[1].assign_coords(pc=['pc1', 'pc2']).rename('comp')
        ])

        for i in range(1, x10.sizes['t']):
            for j in range(2):
                x12 = x10.trans[i-1,:,j]
                x13 = x10.trans[i,:,j]
                x14 = [
                    ((x12+x13)**2).sum().data,
                    ((x12-x13)**2).sum().data
                ]
                x14 = 2*np.argmin(x14)-1
                x10.trans[i,:,j] = x14*x13
                x10.comp[i,j,:] = x14*x10.comp[i,j,:]
    
        x1 = self.analysis.donor[['donor', 'dsm_severity_score', 'dsm_severity_score_group', 'status']]
        x1 = x1.to_dataframe().drop_duplicates()
        x2 = x10.trans.to_dataframe().reset_index()
        x2 = x2.pivot_table(index=['t', 'donor'], columns='pc', values='trans').reset_index()
        x2 = x2.merge(x1, on='donor')
        x2 = x2[x2.dsm_severity_score_group!='']

        x3 = x10.comp.to_dataframe().reset_index()
        x3 = x3.pivot_table(index=['t', 'cytokine'], columns='pc', values='comp').reset_index()

        return x2, x3

    def plot4_data(self, n):
        x2, x6 = self.fit_level4_pca(n)

        x4 = np.quantile(x2.pc1, [0,1])    
        x5 = np.quantile(x2.pc2, [0,1])

        x3 = {t: x for t, x in x2.groupby('t')}

        x8 = {
            t: pd.Categorical(x.donor, categories=x.sort_values('pc1').donor).dtype
            for t,x in x2.groupby('t')
        }
        x9 = {
            t: pd.Categorical(x.cytokine, categories=x.sort_values('pc1').cytokine).dtype
            for t,x in x6.groupby('t')
        }
        x7 = self.fit_level4
        x7 = x7.sel(cytokine=x6.cytokine.drop_duplicates().to_list())
        x7 = x7.level.to_dataframe().reset_index()
        x7['level'] = x7.level - x7.groupby(['t', 'cytokine'])['level'].transform('mean')
        x7 = {
            t: x.assign(
                donor=lambda x: x.donor.astype(x8[t]),
                cytokine=lambda x: x.cytokine.astype(x9[t])
            )
            for t, x in x7.groupby('t') 
        }

        x11 = {
            t:x.assign(
                donor=lambda x: x.donor.astype(x8[t])
            )
            for t, x in x2[['t', 'donor', 'dsm_severity_score_group']].groupby('t') 
        }

        return [
            [x3, x4, x5],
            x7, x11
        ]


app4 = _app4()

#%%
if __name__ == '__main__':
    from plotnine import *
    from plotnine.animation import PlotnineAnimation
    import matplotlib.pyplot as plt

    self = app4

#%%
    self.cytokines

#%%
    x5 = self.plot1_data('ST2/IL-33R')

    for x in x5[0]:
        print(x)

    print(
        ggplot()+
            geom_line(
                data=x5[1],
                mapping=aes('t', 'level', group='donor', color='group'),
                alpha=0.5
            )+
            geom_line(
                data=x5[2],
                mapping=aes('t', 'level', group='donor', color='group'),
                size=2
            )+
            facet_grid('resp~.')+
            scale_color_manual(['red', 'blue', 'black', 'black'])
    )   

#%%
    x5 = self.plot1_data('ST2/IL-33R')
    x5 = x5[1]
    x5 = x5[x5.resp==0]
    c, d = 'ST2/IL-33R', x5.donor.iloc[3]
    x6, x7 = self.plot2_data(c, d)

    print(
        ggplot(x6)+
            aes('days_since_onset', 'level')+
            geom_point(aes(color=f'donor=="{d}"'))+
            geom_line(aes(group='donor', color=f'donor=="{d}"'))+
            geom_line(data=x7[x7.variable=='pred1'], linetype='dotted')+
            geom_line(data=x7[x7.variable=='pred2'], linetype='dashed')+
            geom_line(data=x7[x7.variable=='pred3'], linetype='solid')+
            theme(legend_position='none')+
            scale_color_manual(['lightgray', 'red'])
    )

#%%
    self.fit_level4_deg

#%% 
    x = self.plot3_data('ST2/IL-33R', 12)

    print(
        ggplot(x)+aes('dsm_severity_score_group', 'level')+
            geom_boxplot()+geom_point()
    )


#%%
    x3, x7, x11  = self.plot4_data(10)
    x3, x4, x5 = x3

    x3 = [
        ggplot(x)+
            aes(
                'pc1', 'pc2',
                color='dsm_severity_score_group'
            )+
            geom_point()+
            lims(x=x4, y=x5)
        for x in x3.values()
    ]
    x3[9]

    x10 = [
        ggplot(x)+
            aes('donor', 'cytokine', fill='np.clip(level, -5, 5)')+
            geom_tile()+
            scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0)+
            theme(axis_text_x=element_blank())
        for x in x7.values()
    ]
    x10[9]

    x12 = [
        ggplot(x)+
            aes('donor', y=1, fill='dsm_severity_score_group')+
            geom_tile()+
            theme(axis_text_x=element_blank(), axis_text_y=element_blank())
        for x in x11.values()
    ]
    x12[9]

    #from matplotlib import rc
    #rc('animation', html='html5')
    #%matplotlib ipympl

    p = PlotnineAnimation(
        x3, interval=100, repeat_delay=500
    )
    p.save('xxx.mp4')
    #p

#%%
    from scipy.cluster.hierarchy import dendrogram

    x1 = self.fit_level4
    x1 = x1.sel(cytokine='IL-1ra/IL-1F3')
    x1 = x1.squeeze().drop('cytokine') 

    x2 = hclust(x1.level, 2.9)

    plt.figure(figsize=(10, 4))
    p = dendrogram(x2.linkage.data, labels=None, color_threshold=x2.threshold)

    x5 = xa.merge([x1.level, x2.ord])
    x5 = x5.sel(donor=x5.ord.to_series().sort_values().reset_index().donor.to_list())
    plt.figure(figsize=(10, 15))
    q = plt.imshow(x5.level.data)
    plt.show()

    x6 = contingency(xa.merge([
        x2.clust, x1.dsm_severity_score_group
    ]).to_dataframe())
    for x in x6:
        print(x)

#%%
    x1 = self.fit_level4
    x1 = x1.sel(cytokine='IL-1ra/IL-1F3')
    x1 = x1.squeeze().drop('cytokine')

    x2 = hclust(x1.level, 2.9)

    x7 = xa.merge([x1[['level', 'dsm_severity_score_group']], x2]).\
        to_dataframe().reset_index()
    print(
        ggplot(x7)+aes('t', 'level', color='dsm_severity_score_group')+
            geom_line(aes(group='donor'))+
            facet_grid('clust~.', scales='free_y')
    )

#%%
    x1 = self.fit_level4
    x1 = x1.sel(cytokine='IL-1ra/IL-1F3')
    x1 = x1.squeeze().drop('cytokine')

    x3 = pd.Categorical(x1.dsm_severity_score_group.to_series(), categories=['DSM_low', 'DSM_high']).codes
    x5 = gmm_score1(x1.level, x3)
    x5 = xa.merge([x5, x1[['level', 'dsm_severity_score_group']]])
    x5 = x5.rename(dsm_severity_score_group='group')

#%%
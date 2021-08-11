from pathlib import Path
from common.defs import lazy_property, lazy_method
from common.dir import Dir, cached_property
import pandas as pd
import numpy as np
import scanpy
import matplotlib.pyplot as plt
import seaborn as sb
import xarray as xa
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

cache = Path('.cache')

def _cat(file):
    with file.open('rt') as file:
        for line in file:
            yield line.rstrip()

class Atlas:
    pass

def _():
    @lazy_property
    def storage(self):
        return Dir(cache/'data')
    Atlas.storage = storage

    def scanpy_read(self, name):
        return scanpy.read(Path(self.storage.path)/f'{name}.h5ad')
    Atlas.scanpy_read = scanpy_read

    Atlas.datasets = ['mgh', 'nih-adaptive', 'nih-innate']

    #@lazy_method(key=lambda name: name)
    def data(self, name):
        return self.scanpy_read(name)
    Atlas.data = data

    @lazy_property
    @cached_property(type=Dir.csv)
    def meta(self):
        meta = [(name, self.data(name).obs) for name in self.datasets]
        meta = [
            data.melt().value_counts().rename('count').reset_index().assign(dataset=name)
            for name, data in meta
        ]
        meta = pd.concat(meta)
        return meta
    Atlas.meta = meta

    @lazy_property
    @cached_property(type=Dir.csv)
    def meta_meta(self):
        meta = self.meta
        meta = meta[['variable', 'value', 'dataset']]
        meta = meta.drop_duplicates()
        meta = meta[['variable', 'dataset']]
        meta = meta.value_counts().rename('count').reset_index()
        meta = meta.query('count>1 & count<100')
        meta = meta.pivot_table(index='variable', columns=['dataset'], values='count', fill_value=0, aggfunc='sum')
        meta = meta.reset_index()
        return meta
    Atlas.meta_meta = meta_meta

    @lazy_property
    def data1(self):
        d5 = {}

        d5['mgh'] = self.data('mgh')
        d5['mgh'] = d5['mgh'][d5['mgh'].obs.time_point=='D0']
        d5['mgh'].obs['cell_type'] = d5['mgh'].obs['Annotation']
        d5['mgh'].obs['covid'] = d5['mgh'].obs.covid.astype('category')
        d5['mgh'].obs['donor'] = d5['mgh'].obs.patient_id
        d5['mgh'] = d5['mgh'][d5['mgh'].obs.cell_type != 'Plasmablast']

        d5['nih-innate'] = self.data('nih-innate')
        d5['nih-innate'] = d5['nih-innate'][d5['nih-innate'].obs.timepoint.isin(['T0', 'HC'])]
        d5['nih-innate'].obs['covid'] = np.where(d5['nih-innate'].obs.disease=='Normal', 0, 1)
        d5['nih-innate'].obs['covid'] = d5['nih-innate'].obs.covid.astype('category')

        d5['nih-adaptive'] = self.data('nih-adaptive')
        d5['nih-adaptive'] = d5['nih-adaptive'][d5['nih-adaptive'].obs.timepoint.isin(['T0', 'HC'])]
        d5['nih-adaptive'].obs['covid'] = np.where(d5['nih-adaptive'].obs.disease == 'Normal', 0, 1)
        d5['nih-adaptive'].obs['covid'] = d5['nih-adaptive'].obs.covid.astype('category')

        _ = [set(v.var.index) for v in d5.values()]
        _ = list(_[0] & _[1] & _[2])
        d5 = { k: v[:,_] for k, v in d5.items() }

        for k, v in d5.items():
            print(k)
            scanpy.pp.normalize_total(v, target_sum=1e4)
            scanpy.pp.log1p(v)

        return d5
    Atlas.data1 = data1

    @lazy_property
    @cached_property(type=Dir.csv)
    def cell_type_genes(self):
        def _(k, v):
            print(k)
            v = v[v.obs.covid==0,:]
            scanpy.tl.rank_genes_groups(v, 'cell_type', method='t-test', use_raw=False)
            d3 = pd.DataFrame(v.uns['rank_genes_groups']['names']).melt()
            d3['variable'] = [k[:3]+'.' + s for s in d3['variable']]
            d3['score'] = pd.DataFrame(v.uns['rank_genes_groups']['scores']).melt()['value']
            _ = np.arange(v.shape[1])
            _ = _.reshape((1, -1))
            _ = np.repeat(_, len(v.uns['rank_genes_groups']['names'].dtype), axis=0)
            _ = _.ravel()
            d3['rk'] = _

            return d3
        d3 = [_(k, v) for k, v in self.data1.items()]
        d3 = pd.concat(d3)
        return d3
    Atlas.cell_type_genes = cell_type_genes

    @lazy_property
    @cached_property(type=Dir.csv)
    def covid_genes(self):
        def _(k, v):
            def _(ct):
                print(k[:3] + '.' + ct)
                v1 = v[v.obs.cell_type==ct,:]
                scanpy.tl.rank_genes_groups(v1, 'covid', method='t-test', use_raw=False)
                d6 = pd.DataFrame(v1.uns['rank_genes_groups']['names']).loc[:,['1']]
                d6 = d6.rename(columns={'1': 'value'})
                d6['variable'] = k[:3] + '.' + ct
                d6['score'] = pd.DataFrame(v1.uns['rank_genes_groups']['scores']).melt()['value']
                d6['rk'] = np.arange(v1.shape[1])
                return d6
            d6 = [_(ct) for ct in v.obs.cell_type.drop_duplicates()]
            d6 = pd.concat(d6)
            return d6
        d6 = [_(k, v) for k, v in self.data1.items()]
        d6 = pd.concat(d6)
        return d6
    Atlas.covid_genes = covid_genes
_()

def _p(d2, **kwargs):
    sb.heatmap(d2, **kwargs)
    plt.subplots_adjust(left=0.2, bottom=0.3, right=0.99, top=0.99)
    ax = plt.gca()
    ax.set_xticks(np.arange(d2.shape[1])+0.5)
    ax.set_xticklabels(d2.columns)
    ax.set_yticks(np.arange(d2.shape[0])+0.5)
    ax.set_yticklabels(d2.index)
    ax.set_ylabel('')
    ax.set_xlabel('')
    plt.gca().tick_params(axis='both', which='major', labelsize=10)

def _p1(d):
    n = int(np.sqrt(d.variable.drop_duplicates().shape).round().item())
    for i, g in enumerate(d.groupby('variable')):
        #i, g = next(enumerate(_.groupby('variable')))
        k, v = g
        plt.subplot(n, n, i+1)
        sb.scatterplot(x='rk', y='score', data=v, color='white')
        y = [np.inf, -np.inf]
        for row in v.itertuples(index=False):
            t = plt.text(row.rk, row.score, row.value, fontsize=9, rotation=90)
            b = t.get_window_extent().transformed(plt.gca().transData.inverted())
            y = [min(y[0], b.ymin), max(y[1], b.ymax)]
        y = [y[0]-0.5, y[1]+0.5]
        ax = plt.gca()
        ax.set_ylim(y)
        ax.set_title(k)
        ax.set_xlabel('ranking' if i/n>=(n-1) else '')
    plt.subplots_adjust(wspace=0.2, hspace=0.3)

def _cor(d3, n):
    d4 = d3.query('rk<@n')
    d4 = d4[['variable', 'value']]
    d4 = d4.join(d4.set_index('value'), on='value', lsuffix='_1', rsuffix='_2')
    d4 = d4[['variable_1', 'variable_2']].value_counts().rename('n').reset_index()
    d4 = d4.pivot_table(index='variable_1', columns='variable_2', values='n', fill_value=0)

    _ = np.array(n-d4)
    _ = squareform(_)
    _ = linkage(_, method='ward')
    _ = leaves_list(_)
    d4 = d4.iloc[_, _]

    return d4

covid19_cell_atlas = Atlas()

self = covid19_cell_atlas

d = self.data1['nih-adaptive']
d = pd.DataFrame(dict(
    cell_type = d.obs.cell_type,
    umap_x = d.obsm['X_umap'][:,0],
    umap_y = d.obsm['X_umap'][:,1]
))
d['cell_type'] = d.cell_type.cat.remove_unused_categories()
sb.scatterplot(data=d, x='umap_x', y='umap_y', hue='cell_type', linewidth=0, size=5)
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()

d1 = self.data1
d1 = [v.obs[['donor', 'covid', 'cell_type']].assign(source=k[:3]) for k, v in d1.items()]
d1 = pd.concat(d1)
d1['cell_type'] = [x+'.'+y for x, y in zip(d1.source, d1.cell_type)]
d1 = d1.value_counts().rename('n').reset_index()
d1 = d1.join(d1.groupby(['source', 'donor', 'covid']).n.sum().rename('m'), on=['source', 'donor', 'covid'])
d1['f'] = d1.n/d1.m
d2 = d1.query('covid==0').groupby('cell_type').f.mean().sort_values().index.to_numpy()
d1['cell_type'] = d1.cell_type.astype('category').cat.reorder_categories(d2)

def _(d1, source):
    d1 = d1[d1.source==source].copy()
    d1['cell_type'] = d1.cell_type.cat.remove_unused_categories()
    sb.catplot(data=d1, kind='box', x='cell_type', y='f', hue='covid', ci='sd', legend=False)
    ax = plt.gca()
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=-45,
        horizontalalignment='left'
    )
    plt.subplots_adjust(left=0.1, bottom=0.4, right=0.7, top=0.9)
    ax.set_ylabel('Fraction of cells')
    ax.set_xlabel('')
    ax.set_ylim([0, 0.8])
    ax.set_title(source)
    ax.legend(loc='upper left', title='covid')
_(d1, 'mgh')
_(d1, 'nih')


_ = self.cell_type_genes
_ = _[_.variable.str.find('mgh')>=0].copy()
_['variable'] = _.variable.str.replace('mgh.', '', regex=False)
_ = _.query('rk<=20')
_p1(_)

_ = _cor(self.cell_type_genes, 100)
_p(_, vmin=0, vmax=100, cmap='Reds', annot=True, fmt='.0f', cbar=False)

_ = _cor(self.cell_type_genes, 100)
_1 = _.columns
_1 = _1.to_series().str.find('mgh')>=0
_ = _.iloc[_1.to_numpy(), (~_1).to_numpy()]
_p(_, vmin=0, vmax=100, cmap='Reds', annot=True, fmt='.0f', cbar=False)

_ = self.covid_genes
_ = _[_.variable.str.find('mgh')>=0].copy()
_['variable'] = _.variable.str.replace('mgh.', '', regex=False)
_ = _.query('rk<=20')
_p1(_)

_ = _cor(self.covid_genes, 100)
_p(_, vmin=0, vmax=100, cmap='Reds', annot=True, fmt='.0f', cbar=False)

d4 = self.covid_genes.query('rk<100')
d4 = d4[['variable', 'value']].copy()
d4['variable'] = [s[:3] for s in d4.variable]
d4 = d4.value_counts().rename('n').reset_index()
d4 = d4.pivot_table(index='value', columns='variable', values='n', fill_value=0)
d4['cross'] = d4.nih * d4.mgh
d4['cross1'] = d4.nih**2 + d4.mgh**2

plt.plot(
    d4.mgh + 0.15*np.random.randn(d4.shape[0]),
    d4.nih + 0.15*np.random.randn(d4.shape[0]),
    'o', markersize=2
)
ax = plt.gca()
ax.set_xlabel('Number of cell types in MGH')
ax.set_ylabel('Number of cell types in NIH')

d4.query('mgh==0')

_ = d4.query('cross>=5')
_ = _.sort_values('cross')
plt.barh(np.arange(_.shape[0]), _.cross)
ax = plt.gca()
ax.set_yticks(np.arange(_.shape[0]))
ax.set_yticklabels(_.index)
ax.set_xlabel('Number of MGH-NIH cell type pairs')
plt.subplots_adjust(left=0.3)

self.covid_genes.query('rk<100').query('value=="OAS3"').sort_values('rk')[['variable', 'rk']]

import ncbi.sql as ncbi
from scipy.stats import hypergeom

go = ncbi.query((
    'select a."Symbol" as gene, b."GO_term" as term'
    ' from gene_info as a' 
	' join gene2go as b on (a."GeneID"=b."GeneID")' 
    ' where b.tax_id=9606'
), schema='homo_sapiens')
go = go[go.gene.isin(self.covid_genes.value.drop_duplicates())]
d5 = d4[d4.index.isin(go.gene)].query('cross>0').reset_index()
d5 = d5.merge(go, how='left', left_on='value', right_on='gene')
d6 = d5.groupby('term').size().rename('x').to_frame()
d6 = d6.join(go.groupby('term').size().rename('m'))
d6['k'] = d5.value.drop_duplicates().shape[0]
d6['n'] = go.gene.drop_duplicates().shape[0]
d6['e'] = (d6.m*d6.k/d6.n).astype(int)
def _(x, n, m, k):
    hg = hypergeom(n, m, k)
    return hg.sf(x)+hg.pmf(x)
d6['p'] = [_(row.x, row.n, row.m, row.k) for row in d6.itertuples()]
d6 = d6.sort_values('p', ascending=True)

d5[d5.term.str.find('interferon')>=0][['value', 'mgh', 'nih', 'cross']].drop_duplicates().sort_values('cross', ascending=False)

def _():
    self.meta
    self.meta_meta

def _():
    d = self.data('mgh')

    d.obs['cluster_cat'] =  d.obs.cluster.astype('category')
    d.obs['covid_cat'] =  d.obs.covid.astype('category')

    scanpy.pp.highly_variable_genes(d, min_mean=0.0125, max_mean=3, min_disp=0.5)

    scanpy.pl.highly_variable_genes(d)

    scanpy.pl.umap(d, color=['covid_cat', 'GSTM1'])

    scanpy.pl.umap(d[d.obs.covid==0], color=['Annotation', 'GSTM1'])

    scanpy.pp.normalize_total(d, target_sum=1e4)

    scanpy.pp.log1p(d)

    d.obs.query('covid==0')

    d.raw = d

    scanpy.tl.rank_genes_groups(d, 'Annotation', method='t-test')

    scanpy.pl.rank_genes_groups(d, n_genes=25, sharey=False)

    pd.DataFrame(d.uns['rank_genes_groups']['names']).head(5)

    d.obs.Annotation.drop_duplicates()
    d1 = d[d.obs.Annotation=='Dendritic Cell']
    scanpy.tl.rank_genes_groups(d1, 'covid_cat', method='t-test')
    scanpy.pl.rank_genes_groups(d1, n_genes=25, sharey=False)
    x1 = pd.DataFrame(d1.uns['rank_genes_groups']['names']).head(5).melt()['value'].drop_duplicates()
    scanpy.pl.stacked_violin(d1, x1, groupby='covid_cat', rotation=90)


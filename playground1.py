from pathlib import Path
from common.defs import lazy_property, lazy_method
from common.dir import Dir, cached_property
import scanpy
import pandas as pd

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

    @lazy_method(key=lambda name: name)
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
        meta = meta[['variable', 'dataset']]
        meta = meta.drop_duplicates()
        meta['x'] = 1
        meta = meta.pivot_table(index='variable', columns=['dataset'], values='x', fill_value=0)
        return meta
    Atlas.meta_meta = meta_meta
_()

covid19_cell_atlas = Atlas()

self = covid19_cell_atlas

self.meta
self.meta_meta
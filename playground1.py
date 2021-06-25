from pathlib import Path
from common.defs import lazy_property
import scanpy

cache = Path('.cache')

def _cat(file):
    with file.open('rt') as file:
        for line in file:
            yield line.rstrip()

class Data:
    pass

def _():
    @lazy_property
    def mgh(self):
        return scanpy.read(cache/'mgh.h5ad')
    Data.mgh = mgh

    @lazy_property
    def nih_adaptive(self):
        return scanpy.read(cache/'nih-adaptive.h5ad')
    Data.nih_adaptive = nih_adaptive

    @lazy_property
    def nih_innate(self):
        return scanpy.read(cache/'nih-innate.h5ad')
    Data.nih_innate = nih_innate
_()

covid19_cell_atlas = Data()

self = covid19_cell_atlas

self.mgh.var

self.nih_innate.var

x = set(self.mgh.var.index) & set(self.nih_innate.var)
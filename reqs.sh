#!/bin/bash

#conda create -p ./env
#conda activate ./env
conda install python==3.7
conda install ipykernel
conda install zarr
conda install xarray
conda install scanpy
conda install plotnine
conda install sparse
conda install dask==2021.10.0
conda install mysql-connector-python

pip install rpy2
pip install --ignore-installe cffi #this might conflict with scanpy
#%%
from dask.distributed import LocalCluster, Client
from dask import delayed, compute

if __name__ == '__main__':
    cluster = LocalCluster(n_workers=2)
    client = Client(cluster)

    f = compute(delayed(lambda i:i)(i) for i in range(2))

    client.close()
    cluster.close()
# %%

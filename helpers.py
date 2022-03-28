#%%
from pathlib import Path

#%%
class _config:
    project='covid19-cell-atlas'
    cache = Path.home()/'.cache'/project
    root = Path(__file__).parent

config = _config()
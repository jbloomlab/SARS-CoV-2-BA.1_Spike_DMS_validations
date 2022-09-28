# Analysis of VSV-G pseudotyped virus neutalization

### Set up Analysis


```python
import itertools
import math
import os
import re
import warnings

from IPython.display import display, HTML

import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import natsort

import numpy as np
import pandas as pd
from plotnine import *
import seaborn

import neutcurve
from neutcurve.colorschemes import CBMARKERS, CBPALETTE

import yaml
```


```python
warnings.simplefilter('ignore')
```

Read config file.


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Set seaborn theme:


```python
theme_set(theme_seaborn(style='white', context='talk', font_scale=1))
plt.style.use('seaborn-white')
```


```python
resultsdir=config['resultsdir']
os.makedirs(resultsdir, exist_ok=True)
```

## Read in data


```python
frac_infect = pd.read_csv(config['mAb_neuts_VSVG'], index_col=0)
```

## Fit Hill curve to data using [`neutcurve`](https://jbloomlab.github.io/neutcurve/)


```python
fits = neutcurve.CurveFits(frac_infect, fixbottom=0, fixtop= True)
```


```python
fitparams = (
        fits.fitParams()
        # get columns of interest
        [['serum', 'ic50', 'ic50_bound','virus']]
        .assign(NT50=lambda x: 1/x['ic50'])        
        )
```


```python
fitparams['ic50_is_bound'] = fitparams['ic50_bound'].apply(lambda x: True if x!='interpolated' else False)
```


```python
fitparams
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serum</th>
      <th>ic50</th>
      <th>ic50_bound</th>
      <th>virus</th>
      <th>NT50</th>
      <th>ic50_is_bound</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>CC65.105</td>
      <td>2.753327</td>
      <td>interpolated</td>
      <td>BA.1</td>
      <td>0.363197</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>CC65.105</td>
      <td>300.000000</td>
      <td>lower</td>
      <td>VSV-G</td>
      <td>0.003333</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NTD_5-7</td>
      <td>0.427354</td>
      <td>interpolated</td>
      <td>BA.1</td>
      <td>2.339978</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NTD_5-7</td>
      <td>96.000000</td>
      <td>lower</td>
      <td>VSV-G</td>
      <td>0.010417</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Ly-CoV1404</td>
      <td>0.005371</td>
      <td>interpolated</td>
      <td>BA.1</td>
      <td>186.188050</td>
      <td>False</td>
    </tr>
    <tr>
      <th>5</th>
      <td>Ly-CoV1404</td>
      <td>5.000000</td>
      <td>lower</td>
      <td>VSV-G</td>
      <td>0.200000</td>
      <td>True</td>
    </tr>
    <tr>
      <th>6</th>
      <td>CC9.104</td>
      <td>2.335419</td>
      <td>interpolated</td>
      <td>BA.1</td>
      <td>0.428189</td>
      <td>False</td>
    </tr>
    <tr>
      <th>7</th>
      <td>CC9.104</td>
      <td>300.000000</td>
      <td>lower</td>
      <td>VSV-G</td>
      <td>0.003333</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



## Plot neut curves for mAbs


```python
# neut with Ly-CoV1404
fig, axes = fits.plotSera(sera = ['Ly-CoV1404', ],
                          viruses=['BA.1','VSV-G'],
                          xlabel='concentration (µg/ml)',
                          ncol=1,
                          widthscale=1,
                          heightscale=1,
                          titlesize=20, labelsize=20, ticksize=14,
                          legendfontsize=20, yticklocs=[0,0.5,1],
                          markersize=5, linewidth=1.5,
                          max_viruses_per_subplot = 8,
                          markers=['o','o']
                         )
plotfile = 'Ly-CoV1404_VSVG.pdf'
print(f"Saving to {plotfile}")
fig.savefig(f'{resultsdir}/{plotfile}', bbox_inches='tight')
```

    Saving to Ly-CoV1404_VSVG.pdf



    
![png](VSVG_neutralization_files/VSVG_neutralization_17_1.png)
    



```python
# neuts with NTD_5-7, CC9.104, and CC65.105
fig, axes = fits.plotSera(sera = ['NTD_5-7','CC9.104', 'CC65.105'],
                          viruses=['VSV-G', 'BA.1'],
                          xlabel='concentration (µg/ml)',
                          ncol=3,
                          sharex=False,
                          widthscale=1,
                          heightscale=1.2,
                          titlesize=20, labelsize=20, ticksize=14,
                          legendfontsize=20, yticklocs=[0,0.5,1],
                          markersize=5, linewidth=1.5,
                          max_viruses_per_subplot = 8,
                          markers=['o','o']
                         )
plotfile = 'NTS_S2_mAbs_VSVG.pdf'
print(f"Saving to {plotfile}")
fig.savefig(f'{resultsdir}/{plotfile}', bbox_inches='tight')
```

    Saving to NTS_S2_mAbs_VSVG.pdf



    
![png](VSVG_neutralization_files/VSVG_neutralization_18_1.png)
    



```python

```

#%%
import numpy as np
import pandas as pd
from plotnine import *
from covid19_cell_atlas._analysis3 import analysis3 as analysis
from covid19_cell_atlas._helpers import round_float
from covid19_cell_atlas.common.caching import compose, lazy
from scipy import stats

#%%
class _report:
    @compose(property, lazy)
    def single_cytokine_fit(self):
        x = analysis.fit_time1
        x = x.to_dataframe().reset_index()
        x = x.sort_values('Pr(>F)')
        x = x[x.q<0.05]
        for c in ['Pr(>F)', 'q']:
            x[c] = round_float(x[c], 2)
        for c in list(set(x.columns)-set(['cytokine', 'Pr(>F)', 'q'])):
            x[c] = round(x[c], 2)
        return x

    @compose(property, lazy)
    def multi_cytokine_fit(self):
        return analysis.fit_time2

    @compose(property, lazy)
    def data1(self):        
        x1 = analysis.cytokines
        x2 = self.single_cytokine_fit
        x1 = x1.sel(cytokine=x2.cytokine.to_list())
        x1 = x1[[
            'level', 'days_since_onset', 'donor', 'status', 
            'dsm_severity_score', 'dsm_severity_score_group'
        ]]
        x1['level'] = np.log1p(x1.level)/np.log(2)
        x1 = x1.to_dataframe().reset_index()
        x1 = x1[~x1['level'].isna()]
        return x1

    @compose(property, lazy)
    def data2(self):
        x1 = self.data1
        x2 = self.multi_cytokine_fit
        x3 = x1.groupby([
            'donor', 'days_since_onset', 'cytokine', 'status', 
            'dsm_severity_score', 'dsm_severity_score_group'
        ])['level'].mean()
        x3 = x3.reset_index()
        x3 = x3.pivot_table(
            index=[
                'days_since_onset', 'donor', 'status', 
                'dsm_severity_score', 'dsm_severity_score_group'
            ], 
            columns='cytokine', values='level'
        )
        x3 = x3.reset_index()
        x3['const'] = 1
        x3['days_since_onset_pred'] = x3[x2.coef.data].to_numpy() @ x2.coef_value.data
        x3['delta_days'] = x3.days_since_onset_pred-x3.days_since_onset
        x3 = x3[~x3.delta_days.isna()]
        return x3

report = _report()

#%%
# Good prediction of time using a single cytokine

report.single_cytokine_fit

#%%
# Cytokine level vs time plots

x1 = report.data1
for x in report.single_cytokine_fit.itertuples():
    x2 = x1[x1.cytokine==x.cytokine][['days_since_onset', 'level']]
    x3 = pd.DataFrame({
        'level': [x2.level.min(), x2.level.max()]
    })
    x3['days_since_onset'] = x.Intercept+x.level*x3['level']
    print(
        ggplot()+aes('level', 'days_since_onset')+
            geom_point(data=x2)+
            geom_line(data=x3)+
            labs(title=x.cytokine)
    )

#%%
# Mutivariate prediction of time using cytoking wich predict time well
print(report.multi_cytokine_fit.stat_value.to_series())

#%%    
print(
    ggplot(report.data2)+
        aes('days_since_onset', 'days_since_onset_pred', color='status')+
        geom_point()+
        geom_abline(slope=1, intercept=0)
)

#%%
print(
    ggplot(report.data2)+
        aes('status', 'delta_days')+
        geom_boxplot()+
        geom_point()
)

#%%
print(
    ggplot(report.data2)+
        aes('dsm_severity_score', 'delta_days', color='dsm_severity_score_group')+
        geom_point()+geom_smooth(method='lm')
)

#%%
x = report.data2
stats.ttest_ind(
    x.delta_days[x.dsm_severity_score_group=='DSM_high'],
    x.delta_days[x.dsm_severity_score_group=='DSM_low']
)

# %%

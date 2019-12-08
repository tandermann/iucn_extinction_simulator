#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:30:52 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)
import os


iucn_history_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/precompiled/iucn_history/MAMMALIA_iucn_history.txt'
iucn_history_df = pd.read_csv(iucn_history_file,sep='\t')
iucn_history_df = iucn_history_df.sort_values(by='species')
iucn_history_df = iucn_history_df.drop_duplicates()
iucn_history_df.index = np.arange(len(iucn_history_df))
iucn_history_df.to_csv(iucn_history_file,sep='\t',index=False)






from urllib.request import urlopen
from io import StringIO

url = 'https://raw.githubusercontent.com/tobiashofmann88/iucn_extinction_simulator/master/data/precompiled/iucn_history/MAMMALIA_iucn_history.txt'
outdir = '/Users/tobias/Desktop'
urlpath =urlopen(url)
string = urlpath.read().decode('utf-8')

string_input = StringIO(string)
pd.read_csv(string_input, sep="\t")



ext_rate_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/carnivora_output/future_simulations/extinction_prob_all_species.txt'
ext_rate_df = pd.read_csv(ext_rate_file,sep='\t')
plt.hist(np.log10(ext_rate_df.rate_e_mean),50)

plt.plot(sorted(ext_rate_df.rate_e_mean))


trans_rates_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/carnivora_output/transition_rates/sampled_status_change_rates.txt'
trans_rates_df = pd.read_csv(trans_rates_file,sep='\t')
states = ['LC','NT','VU','EN','CR','DD']
trans_rates_matrix = pd.DataFrame(data=np.zeros([6,6]),index=states,columns=states)
for i,s_string in enumerate(trans_rates_df.status_change.values):
    values = trans_rates_df.iloc[i,1:].values.astype(float)
    min_value = min(values)
    max_value = max(values)
    state_a,state_b = s_string.split('->')
    trans_rates_matrix.loc[state_a,state_b] = '%.4f-%.4f'%(min_value,max_value)
print(trans_rates_matrix.to_latex())

np.max(np.mean(trans_rates_df.iloc[:,1:].values,axis=1))



n100_sim_file = '/Users/tobias/Desktop/runsim_test/100sim_100years/extinction_prob_all_species.txt'
n10k_sim_file = '/Users/tobias/Desktop/runsim_test/10ksim_100years/extinction_prob_all_species.txt'
n100_t10k_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/carnivora_output/future_simulations_10k_years/extinction_prob_all_species.txt'

n100_sim = pd.read_csv(n100_sim_file,sep='\t').rate_e_mean.values
uncertainty_interval_n100_sim = np.abs(n100_sim-pd.read_csv(n100_sim_file,sep='\t').iloc[:,2:4].values.T)
n10k_sim = pd.read_csv(n10k_sim_file,sep='\t').rate_e_mean.values
uncertainty_interval_n10k_sim = np.abs(n10k_sim-pd.read_csv(n10k_sim_file,sep='\t').iloc[:,2:4].values.T)

n100_t10k_sim = pd.read_csv(n100_t10k_file,sep='\t').rate_e_mean.values

fig = plt.figure()
#plt.plot(np.log10(n100_sim),np.log10(n10k_sim),'o')
plt.errorbar(n100_sim,n10k_sim, xerr=uncertainty_interval_n100_sim, yerr=uncertainty_interval_n10k_sim, fmt='o',ecolor='black')
plt.xlim(0.00005,0.05)
plt.ylim(0.00005,0.05)
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('np.log(n100_sim)')
plt.ylabel('np.log(n10k_sim)')
plt.tight_layout()
fig.savefig('/Users/tobias/Desktop/plot_log.pdf',bbox_inches='tight', dpi = 500)








def calcHPD(data, level):
    assert (0 < level < 1)
    d = list(data)
    d.sort()
    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
        raise RuntimeError("not enough data")
    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)):
         rk = d[k+nIn-1] - d[k]
         if rk < r :
             r = rk
             i = k
    assert 0 <= i <= i+nIn-1 < len(d)
    return (d[i], d[i+nIn-1])


max_t = 100
x = np.array([33,12,45,100,100,54,100,33,44,6])
n_bins = 10000
n_samples=100

q_samples = np.linspace(0, 1, n_bins)
all_sampled_rates = []
for x_value in x:
    if x_value < max_t:
        #extinct        
        p = q_samples*np.exp(-q_samples*x_value)
    else:
        # not extinct
        p = np.exp(-q_samples*x_value)
    p = p/sum(p)
    sample_rates = np.random.choice(q_samples,size=n_samples,p=p,replace=1)
    all_sampled_rates.append(sample_rates)
all_rates = [item for sublist in all_sampled_rates for item in sublist]
log_rates = np.log(all_rates)

all_rates[all_rates<=0]


np.mean(all_rates)
calcHPD(all_rates,0.95)

plt.hist(all_rates,100)

np.random.random(n_bins)

# 10000 q-s from uniform between 0 and 1
# calculate posterior for all q-s and each given x-value
# resample 100 q-values with probability=posterior --> needs to be normalized, summing up to 1, replace=True






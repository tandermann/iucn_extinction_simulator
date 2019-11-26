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



from urllib.request import urlopen
import re

url = 'https://raw.githubusercontent.com/tobiashofmann88/iucn_extinction_simulator/master/data/precompiled/iucn_history/MAMMALIA_iucn_history.txt'
outdir = '/Users/tobias/Desktop'
urlpath =urlopen(url)
string = urlpath.read().decode('utf-8')
print(string)

pattern = re.compile('[A-Z]*_iucn_history.txt"') #the pattern actually creates duplicates in the list

filelist = pattern.findall(string)
print(filelist)

for filename in filelist:
    filename=filename[:-1]
    remotefile = urlopen(url + filename)
    localfile = open(os.path.join(outdir,filename),'wb')
    localfile.write(remotefile.read())
    localfile.close()
    remotefile.close()


n100_sim_file = '/Users/tobias/Desktop/extinction_prob_all_species_100.txt'
n10k_sim_file = '/Users/tobias/Desktop/extinction_prob_all_species_10k.txt'
n100_t10k_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/output_dir/simulations_100sim_10kyears/extinction_prob_all_species.txt'

n100_sim = pd.read_csv(n100_sim_file,sep='\t').rate_e_mean.values
n10k_sim = pd.read_csv(n10k_sim_file,sep='\t').rate_e_mean.values
n100_t10k_sim = pd.read_csv(n100_t10k_file,sep='\t').rate_e_mean.values

fig = plt.figure()
plt.plot(np.log10(n100_t10k_sim),np.log10(n10k_sim),'o')
plt.xlabel('np.log10(n100_t10k_sim)')
plt.ylabel('np.log10(n10k_sim)')
plt.tight_layout()
plt.xlim([-2.6,-1.0])
plt.ylim([-2.6,-1.0])
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






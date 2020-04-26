#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 14:40:30 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import os,sys
import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)

def update_multiplier(q,d=1.1):
    u = np.random.uniform(0,1)
    l = 2*np.log(d)
    m = np.exp(l*(u-.5))
    new_q = q * m
    return new_q, np.log(m)

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

def get_rate_estimate(ext_time_array,max_t,index,species_list,plot_posterior=0,pdf=0,n_gen = 100000,burnin = 1000):
    sys.stdout.write('\rProcessing species: %i/%i '%(index+1,len(species_list)))
    ext_time_array_new = ext_time_array.copy()
    ext_time_array_new[ext_time_array_new!=ext_time_array_new] = max_t
    ext_time_array_new = ext_time_array_new.astype(float)
    w_times = np.sum(ext_time_array_new)
    ext_events = len(ext_time_array_new[ext_time_array_new<max_t])
    post_samples = []
    q = 0.01
    likA = np.log(q)*ext_events -q*w_times    
    for i in range(n_gen):
        new_q, hast = update_multiplier(q)
        lik = np.log(new_q)*ext_events -new_q*w_times
        if lik-likA + hast >= np.log(np.random.random()):
            q = new_q
            likA = lik
        if i > burnin and i % 10==0:
            post_samples.append(q)
    mean_value = np.mean(post_samples)
    lower,upper = calcHPD(post_samples,0.95)
    if plot_posterior:
        plt.figure()
        plt.hist(post_samples,100)
        plt.xlabel('Extinction rate estimates')
        plt.ylabel('Counts')
        plt.title(species_list[index])
        plt.tight_layout()
        pdf.savefig()
        plt.close()
        #fig.savefig(os.path.join(posterior_plot_dir,'%s.pdf'%species_list[index]),bbox_inches='tight', dpi = 500)
    return [mean_value,lower,upper]


# ________________________________Simulating data______________________________
n_sim = 100
run_estimation = False

n_species = 1000
max_t = 100
n_gen = 10000
burnin = 1000
outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/simulations'
sim_scenario_string = 'n_%i'%n_sim
sim_species_list = ['Genus species%i'%i for i in np.arange(n_species)]
lc_yearly = 0.000000155728
cr_yearly = 0.066967008463
if run_estimation:
    true_rates = np.exp(np.random.uniform(np.log(lc_yearly),np.log(cr_yearly),n_species))
    all_estimated_rates = []
    for i,species in enumerate(sim_species_list):
        true_rate = true_rates[i]
        ext_time_array = (np.random.exponential(1./true_rate, n_sim)).astype(int)
        ext_time_array[ext_time_array>max_t] = max_t
        estimated_rates = get_rate_estimate(ext_time_array,max_t,i,sim_species_list,plot_posterior=0,pdf=0,n_gen = n_gen,burnin = burnin)
        all_estimated_rates.append(estimated_rates)
    all_joined = np.vstack([true_rates,np.array(all_estimated_rates).T]).T
    output_df = pd.DataFrame(all_joined,columns=['simulated_rates','estimated_rates_mean','estimated_rates_lower','estimated_rates_upper'])
    output_df.to_csv(os.path.join(outdir,'simulated_extinction_rates_%s.txt'%sim_scenario_string),sep='\t',index=False)

# _______________________________Plotting______________________________________
output_df = pd.read_csv(os.path.join(outdir,'simulated_extinction_rates_%s.txt'%sim_scenario_string),sep='\t')
simulated = output_df.values
true_rates =  simulated[:,0]
simulated_mean = simulated[:,1]
simulated_delta_lower = simulated_mean-simulated[:,2]
simulated_delta_upper = simulated[:,3]-simulated_mean

fig = plt.figure(figsize=[4,4])
plt.errorbar(true_rates,simulated_mean, yerr=np.array([simulated_delta_lower,simulated_delta_upper]), fmt='.',ecolor='grey',elinewidth=0.1,capsize=0.)
plt.plot(true_rates,true_rates,'-',color='red',zorder=10)
plt.axhline(0.00001,linestyle='--',color='black',alpha=0.5)
plt.xlim(min(true_rates)-0.4*min(true_rates),max(true_rates)+0.4*max(true_rates))
plt.ylim(min(true_rates)-0.4*min(true_rates),max(true_rates)+0.4*max(true_rates))
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('True extinction rates')
plt.ylabel('Estimated extinction rates')
plt.tight_layout()
fig.savefig(os.path.join(outdir,'extinction_rate_estimates_%s.pdf'%sim_scenario_string),bbox_inches='tight', dpi = 500)

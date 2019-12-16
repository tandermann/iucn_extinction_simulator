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
import os,glob




# Predicted extinctions figure
fig = plt.figure(figsize=(6,4))
data_files = glob.glob('/Users/tobias/GitHub/iucn_predictions/doc/figures/estimated_n_extinctions_arrays/*.txt')
data_array = [np.loadtxt(i) for i in data_files]
order = [2,0,3,1]
labels = np.array([os.path.basename(i).replace('.txt','').replace('_',' ').capitalize() for i in data_files])[order]
labels = ['GL + SC','no GL + SC','GL + no SC','no GL + no SC']
plt.boxplot([data_array[i][-1] for i in order],labels=labels,showfliers=False,patch_artist=True);
plt.ylabel('Extinct species')
plt.title('Simulated bird extinctions in next 100 years')
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/estimated_n_extinctions_arrays/predictions_extinction.pdf',bbox_inches='tight', dpi = 500)





# Figure 3
#_________________separate_______________________
status_through_time_file = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/figure_data/status_distribution_through_time.txt'
status_through_time_data = pd.read_csv(status_through_time_file,sep='\t')
#np.sum(status_through_time_data.iloc[:,1:].values,axis=1)
diversity_through_time_file = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/figure_data/future_diversity_trajectory.txt'
diversity_through_time = np.loadtxt(diversity_through_time_file)
# define time axis
time_axis = status_through_time_data.year.values
fig = plt.figure(figsize=(4,4))
y_values = diversity_through_time[0]
plt.plot(time_axis,y_values,color="#b80033", label='accounting for GL')
# get upper and lower confidence interval boundaries
min_hpd, max_hpd = diversity_through_time[1],diversity_through_time[2]
plt.fill_between(time_axis, min_hpd, max_hpd,
         color="#b80033", alpha=0.7)
#plt.legend()
#plt.ylabel('Total diversity')
#plt.xlabel('Years from present')
#ax = plt.gca()
#ax1 = ax.twinx()
# Set the limits of the new axis from the original axis limits
#ax1.set_ylim(ax.get_ylim())
#current_diversity = diversity_through_time[0,0]
#plt.yticks([np.mean(diversity_through_time[:,-1])],[int(current_diversity-np.mean(diversity_through_time[:,-1]))])
#plt.xticks(modified_q_matrix.year[::10],modified_q_matrix.year[::10])
#plt.ylabel('Lost species')
#plt.tight_layout()
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/div_plot.pdf',bbox_inches='tight', dpi = 500)

colors = np.array(["#227a00","#a5c279","#f3d248","#6956cb","#79262a","#b80033"])
labels = np.array(['LC', 'NT', 'VU', 'EN', 'CR', 'EX'])
def func(pct, allvals):
    absolute = int(np.round((pct/100.*np.sum(allvals))))
    return "{:d}".format(absolute)

fig = plt.figure(figsize=(4,4))
start_status_distr = status_through_time_data.iloc[0,1:].values
# status distribution beginning
wedges, texts, autotexts = plt.pie(start_status_distr[start_status_distr >0], colors= colors[start_status_distr >0], autopct=lambda pct: func(pct, start_status_distr[start_status_distr >0]), shadow=False,textprops=dict(color="w"))
#plt.legend(labels=final_labels[0],title="IUCN status\n(N=10965 sp.)",loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/pie_beginning.pdf',bbox_inches='tight', dpi = 500)

fig = plt.figure(figsize=(4,4))
end_status_distr = status_through_time_data.iloc[-1,1:].values
# status distribution end
wedges, texts, autotexts = plt.pie(end_status_distr[end_status_distr >0], colors= colors[end_status_distr >0], autopct=lambda pct: func(pct, end_status_distr[end_status_distr >0]), shadow=False,textprops=dict(color="w"))
final_labels = list(labels[end_status_distr[0] >0])
plt.legend(labels=final_labels[0],title="IUCN status\n(N=10965 sp.)",loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/pie_end.pdf',bbox_inches='tight', dpi = 500)


#_________________joined_______________________
fig, axs = plt.subplots(1, 3,figsize=(12,4))
y_values = diversity_through_time[0]
axs[0].plot(time_axis,y_values,color="#b80033", label='accounting for GL')
axs[0].set_xlabel('Years from present')
axs[0].set_ylabel('Species diveristy')
axs[0].set_ylim(10000,11100)
axs[0].grid(which='major', alpha=0.5)
#axs[0].set_title('Diversity trajectory Aves')
# get upper and lower confidence interval boundaries
min_hpd, max_hpd = diversity_through_time[1],diversity_through_time[2]
axs[0].fill_between(time_axis, min_hpd, max_hpd,color="#b80033", alpha=0.7)
start_status_distr = status_through_time_data.iloc[0,1:].values
end_status_distr = status_through_time_data.iloc[-1,1:].values
# status distribution beginning
wedges, texts, autotexts = axs[1].pie(start_status_distr[start_status_distr >0], colors= colors[start_status_distr >0], autopct=lambda pct: func(pct, start_status_distr[start_status_distr >0]), shadow=False,textprops=dict(color="w"))
final_labels = list(labels[start_status_distr[0] >0])
axs[1].set_title('Status distribution present')
# status distribution end
wedges, texts, autotexts = axs[2].pie(end_status_distr[end_status_distr >0], colors= colors[end_status_distr >0], autopct=lambda pct: func(pct, end_status_distr[end_status_distr >0]), shadow=False,textprops=dict(color="w"))
final_labels = list(labels[end_status_distr[0] >0])
axs[2].legend(labels=final_labels[0],title="IUCN status\n(N=10965 sp.)",loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
axs[2].set_title('Status distribution in 100 years')
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_3/all_joined.pdf',bbox_inches='tight', dpi = 500)




# plot Figure 2
post_sample_files_path = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_2/figure_data/posterior_samples/*.txt'
post_sample_files = glob.glob(post_sample_files_path)
species_file_dict = {}
for i in post_sample_files:
    filename = os.path.basename(i)
    species = '_'.join(filename.split('_')[:2])
    species_file_dict.setdefault(species,[])
    species_file_dict[species].append(i)

species,status = 'Cathartes_aura','LC'
#species,status = 'Sarcogyps_calvus','CR'
scenario = 'SC'
if scenario == 'SC':
    indeces = [0,1]
elif scenario == 'GL':
    indeces = [0,2]    
files = species_file_dict[species]
gl_ss_posterior = np.loadtxt([file for file in files if '%s_gl_status_change.txt'%species in file][0])
gl_posterior = np.loadtxt([file for file in files if '%s_gl_no_status_change.txt'%species in file][0])
ss_posterior = np.loadtxt([file for file in files if '%s_no_gl_status_change.txt'%species in file][0])
posterior = np.loadtxt([file for file in files if '%s_no_gl_no_status_change.txt'%species in file][0])
data = [gl_ss_posterior,gl_posterior,ss_posterior,posterior]
colors = ["#965da7","#84a955","#bc5d41"]
labels = ['GL + SC','GL + no SC','no GL + SC']
fig = plt.figure(figsize=(4,4))
plt.hist(data[indeces[0]], weights=np.ones(len(data[indeces[0]]))/len(data[indeces[0]]),density = False,color = colors[indeces[0]],label=labels[indeces[0]],bins=100,alpha=0.5);
plt.hist(data[indeces[1]], weights=np.ones(len(data[indeces[1]]))/len(data[indeces[1]]),density = False,color = colors[indeces[1]],label=labels[indeces[1]],bins=100,alpha=0.5);
#plt.hist(data[2], color = colors[2],label='GL + no SC',bins=100,alpha=0.5);
#sns.distplot(data[0], hist = True, color = '#79262a',kde = False,label='GL + SC')
#sns.distplot(data[1], hist = True, color = '#e34349',kde = False,label='GL + no SC')
#sns.distplot(data[2], hist = False, kde = True,kde_kws = {'shade': True, 'linewidth': 3},label='no GL + SC')    
#sns.distplot(data[3], hist = False, kde = True,kde_kws = {'shade': True, 'linewidth': 3},label='no GL + no SC')
plt.axvline(np.mean(data[indeces[0]]),linestyle='--',color=colors[indeces[0]],alpha=0.7)
plt.axvline(np.mean(data[indeces[1]]),linestyle='--',color=colors[indeces[1]],alpha=0.7)
#plt.axvline(np.mean(data[2]),linestyle='--',color=colors[2],alpha=0.7)
plt.title(species.replace('_',' ')+' (%s)'%status)
plt.legend(loc='upper right')
if status == 'LC' and scenario == 'GL':
    plt.xticks([0.0001,0.00015,0.0002],[0.0001,0.00015,0.0002])
    plt.xlim(0.00008,0.00022)
if status == 'LC' and scenario == 'SC':
    plt.xticks([0.0001,0.00015,0.0002],[0.0001,0.00015,0.0002])
    plt.xlim(0.00008,0.00022)
if status== 'CR' and scenario == 'SC':
    plt.xticks([0.015,0.0155,0.016],[0.015,0.0155,0.016])
if status== 'CR' and scenario == 'GL':
    plt.xlim(0.013,0.0255)
    plt.xticks([0.016,0.020,0.024],[0.016,0.020,0.024])
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_2/figure_data/%s.pdf'%(species+'_'+scenario),bbox_inches='tight', dpi = 500)











    
birds_ext_rates_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/birds_output/future_simulations_backup/extinction_prob_all_species.txt'
birds_ext_rates = pd.read_csv(birds_ext_rates_file,sep='\t')
birds_ext_rates_mean = birds_ext_rates.rate_e_mean.values

fig = plt.figure(figsize=[3,3])
plt.hist(np.log10(birds_ext_rates_mean),50);
plt.xticks(np.log10([0.0001,0.001,0.01]),['$10^{-4}$','$10^{-3}$','$10^{-2}$'])
#plt.xlabel('Estimated extinction rates (annual)')
#plt.ylabel('Count')
#plt.ylim(0,1000)
plt.tight_layout()
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/birds_estimated_extinction_rates_hist.pdf',bbox_inches='tight', dpi = 500)



# correct iucn history file, in case multiple rows of same species are present
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






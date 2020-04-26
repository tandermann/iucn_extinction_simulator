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
import os,glob,sys
import iucn_sim.functions as cust_func



def update_multiplier(q,d=1.1):
    u = np.random.uniform(0,1)
    l = 2*np.log(d)
    m = np.exp(l*(u-.5))
    new_q = q * m
    return new_q, np.log(m)

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
    lower,upper = cust_func.calcHPD(post_samples,0.95)
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




# Predicted extinctions figure
fig = plt.figure(figsize=(8,4))
data_files = glob.glob('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_6/*.txt')
data_array = [np.loadtxt(i) for i in data_files]
order = [4,0,5,1,3,2]
labels = np.array([os.path.basename(i).replace('.txt','').replace('_',' ').capitalize() for i in data_files])[order]
labels = ['GL+SC','no GL+SC','GL+no SC','no GL+no SC','EX mode 1','EX mode 1(PEX)']
plt.boxplot([data_array[i][-1] for i in order],labels=labels,showfliers=False,patch_artist=True);
plt.ylabel('Extinct species')
plt.title('Simulated bird extinctions in next 100 years')
fig.savefig('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_6/predictions_extinction.pdf',bbox_inches='tight', dpi = 500)




# q-matrix table
qmatrix_list_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/iucn_sim_output/aves/transition_rates_1/simulation_input_data.pkl'
input_data = cust_func.load_obj(qmatrix_list_file)
species_input_data, dd_probs = input_data
species_list = np.array([i[0] for i in species_input_data])
current_status_list = np.array([i[1] for i in species_input_data])
q_matrix_list = [i[2] for i in species_input_data]
len()
all_qmatrices = np.array([item for sublist in q_matrix_list for item in sublist])
np.mean(all_qmatrices,axis=0)




# average extinction rate - estimate the rates from the simulated extinction times of all species
te_array_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/iucn_sim_output/aves/future_sim_0/te_all_species.txt'
te_array = pd.read_csv(te_array_file,sep='\t',header=None).values
ext_date_data = te_array[:,1:].copy()

species_values = np.reshape(ext_date_data,[ext_date_data.shape[0]*ext_date_data.shape[1],])
sim_species_list = ['all birds']
i = 0
n_years = 100
rate_estimates = get_rate_estimate(species_values,n_years,i,sim_species_list,plot_posterior=True,pdf=0,n_gen=100000,burnin=1000)

rate_estimates


# EX mode 0
10961-10222.58
10961-10152.00
10961-10292.00

# EX mode 1
10961-10254.04
10961-10182.00
10961-10317.00

# EX mode 2, no PEX
10961 - 10905.31
10961 - 10868.00
10961 - 10940.00

# EX mode 2
10961-10834.41
10961-10779.00
10961-10879.00

EX mode 0, 500 years
10961-7570.70
10961-7245.00
10961-8018.00






# test rate and see how many exitncitons it produces
ext_species = []
for i in np.arange(10000):
    true_rate = 0.000217
    n_sim = 10961
    ext_time_array = (np.random.exponential(1./true_rate, n_sim)).astype(int)
    max_t=500
    ext_time_array[ext_time_array>max_t] = 0
    values,counts = np.unique(ext_time_array,return_counts=True)
    n_ext = 10961-counts[0]
    ext_species.append(n_ext)
np.mean(ext_species)
cust_func.calcHPD(ext_species,0.95)



# compare extinction rates:
highlight_dd_ne = True
plot_species = 'Strigops habroptila'
ex_mode_1_file  = '/Users/tobias/GitHub/iucn_extinction_simulator/data/iucn_sim_output/aves//future_sim_0/extinction_prob_all_species.txt'
ex_mode_2_file  = '/Users/tobias/GitHub/iucn_extinction_simulator/data/iucn_sim_output/aves//future_sim_1/extinction_prob_all_species.txt'
ex_mode_1 = pd.read_csv(ex_mode_1_file,sep='\t').rate_e_mean.values
ex_mode_2 = pd.read_csv(ex_mode_2_file,sep='\t').rate_e_mean.values
species_names = pd.read_csv(ex_mode_2_file,sep='\t').species.values
status_info_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/iucn_sim_output/aves/iucn_data/species_data.txt'
status_array = pd.read_csv(status_info_file,sep='\t',header=None).iloc[:,1].values
species_list_status_array = pd.read_csv(status_info_file,sep='\t',header=None).iloc[:,0].values

colors = np.array(['#1f77b4',[0, 0, 0], [1, 0.64, 0]])
col_index = np.zeros(len(species_names)).astype(int)
col_index[status_array=='DD'] = 1
#col_index[status_array=='NE'] = 2

fig = plt.figure(figsize=(4,4))
if highlight_dd_ne:
    filename = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_4/mode_1_vs_mode_2_with_dd_and_ne.pdf'
    plt.scatter(ex_mode_1,ex_mode_2,s=1,color=list(colors[col_index]))
else:
    filename = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_4/mode_1_vs_mode_2.pdf'
    plt.scatter(ex_mode_1,ex_mode_2,s=1)
plt.xlabel('$\mu_i$ (EX mode 0)')
plt.ylabel('$\mu_i$ (EX mode 1 + PEX taxa)')
min_value = min(ex_mode_2)-5e-6
max_value = max(ex_mode_1)+1e-2
plt.xlim(min_value,max_value)
plt.ylim(min_value,max_value)
plt.plot(np.linspace(min_value,max_value),np.linspace(min_value,max_value),color='red',linewidth=1)
if plot_species:
    plt.scatter(ex_mode_1[species_names==plot_species],ex_mode_2[species_names==plot_species],s=8,marker='o',color='green',edgecolors='black',linewidth=0.5)
    
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
#plt.title('$\mu$ estimates for individual bird species')
fig.savefig(filename,bbox_inches='tight', dpi = 500)




# Figure 2
indir = '/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_2/figure_data/mode_0'
#_________________separate_______________________
status_through_time_file = os.path.join(indir,'status_distribution_through_time.txt')
status_through_time_data = pd.read_csv(status_through_time_file,sep='\t')
#np.sum(status_through_time_data.iloc[:,1:].values,axis=1)
diversity_through_time_file = os.path.join(indir,'future_diversity_trajectory.txt')
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
fig.savefig(os.path.join(indir,'div_plot.pdf'),bbox_inches='tight', dpi = 500)
#
colors = np.array(["#227a00","#a5c279","#f3d248","#6956cb","#79262a","#b80033"])
labels = np.array(['LC', 'NT', 'VU', 'EN', 'CR', 'EX'])
def func(pct, allvals):
    absolute = int(np.round((pct/100.*np.sum(allvals))))
    return "{:d}".format(absolute)
#
fig = plt.figure(figsize=(4,4))
start_status_distr = status_through_time_data.iloc[0,1:].values
# status distribution beginning
wedges, texts, autotexts = plt.pie(start_status_distr[start_status_distr >0], colors= colors[start_status_distr >0], autopct=lambda pct: func(pct, start_status_distr[start_status_distr >0]), shadow=False,textprops=dict(color="w"))
#plt.legend(labels=final_labels[0],title="IUCN status\n(N=10965 sp.)",loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
fig.savefig(os.path.join(indir,'pie_beginning.pdf'),bbox_inches='tight', dpi = 500)
#
fig = plt.figure(figsize=(4,4))
end_status_distr = status_through_time_data.iloc[-1,1:].values
# status distribution end
wedges, texts, autotexts = plt.pie(end_status_distr[end_status_distr >0], colors= colors[end_status_distr >0], autopct=lambda pct: func(pct, end_status_distr[end_status_distr >0]), shadow=False,textprops=dict(color="w"))
final_labels = list(labels[end_status_distr[0] >0])
plt.legend(labels=final_labels[0],title="IUCN status\n(N=10965 sp.)",loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
fig.savefig(os.path.join(indir,'pie_end.pdf'),bbox_inches='tight', dpi = 500)
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
fig.savefig(os.path.join(indir,'all_joined.pdf'),bbox_inches='tight', dpi = 500)





# GL test
bird_et_al_file = '/Users/tobias/GitHub/iucn_predictions/data/raw/generation_length/birds/cobi13486-sup-0004-tables4.xlsx'
bird_data = pd.read_excel(bird_et_al_file)
bird_species_list = bird_data['Scientific name'].values
bird_gl_data = bird_data.GenLength.values

tobi_gl_data_file = '/Users/tobias/GitHub/iucn_extinction_simulator/data/precompiled/gl_data/andermann/aves_gl.txt'
tobi_data = pd.read_csv(tobi_gl_data_file,sep='\t',header=None)
tobi_species_list = tobi_data[0].values
tobi_species_list
tobi_gl_data = np.mean(tobi_data.iloc[:,1:].values,axis=1)

bird_gl_data_sorted = [bird_gl_data[bird_species_list==i][0] for i in tobi_species_list]
plt.scatter(tobi_gl_data,bird_gl_data_sorted)    


bird_gl_data_output_df = pd.DataFrame(np.array([tobi_species_list,bird_gl_data_sorted]).T)
bird_gl_data_output_df.to_csv('/Users/tobias/GitHub/iucn_extinction_simulator/data/precompiled/gl_data/aves_gl.txt',sep='\t',index=False,header=False)


zip()



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






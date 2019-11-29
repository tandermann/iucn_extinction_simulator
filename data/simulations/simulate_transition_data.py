#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:21:18 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import os
import numpy as np
np.set_printoptions(suppress=True)
np.random.seed(1234)
np.set_printoptions(suppress= 1) # prints floats, no scientific notation
np.set_printoptions(precision=5)
import matplotlib.pyplot as plt
import pandas as pd

def p_e_year(years,p_e):
    pe_year = 1-(1-float(p_e))**(1/years)
    return pe_year

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

def update_multiplier(q,d=1.1):
    u = np.random.uniform(0,1)
    l = 2*np.log(d)
    m = np.exp(l*(u-.5))
    new_q = q * m
    return new_q, np.log(m)

def sample_rate_mcmc(count, tot_time, n_gen = 100000,burnin = 1000):
    def get_loglik(count, dT, rate):
        return np.log(rate)*count - dT*rate
    post_samples = []
    q = 0.01
    likA = get_loglik(count,tot_time,q)    
    for i in range(n_gen):
        new_q, hast = update_multiplier(q)
        lik = get_loglik(count,tot_time,new_q) 
        if lik-likA + hast >= np.log(np.random.random()):
            q = new_q
            likA = lik
        if i > burnin and i % 10==0:
            post_samples.append(q)
    sampled_rates_mean = np.mean(post_samples)
    lower,upper = calcHPD(post_samples,0.95)
    return [sampled_rates_mean,lower,upper]



# ______________________________Simulate data__________________________________
size_ref_group = 100000
run_estimation = False

n_gen = 100000
burnin = 10000
outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/simulations'
sim_scenario_string = 'refgroupsize_%i'%size_ref_group

if run_estimation:
    # sample rates for simulation from empirical range (Mammalia)
    mammal_min_rate = 0.00004
    mammal_max_rate = 0.0145
    rate_samples = np.exp(np.random.uniform(np.log(mammal_min_rate),np.log(mammal_max_rate),36))
    true_rates = rate_samples.reshape([6,6])
    np.fill_diagonal(true_rates,0)
    
    # sample starting status from empirical status distribution
    states = ['LC','NT','VU','EN','CR','DD']
    mammalia_status_distribution = np.array([3308,348,531,500,203,861])
    #mammalia_status_distribution = np.array([1,1,1,1,1,1])
    probs = mammalia_status_distribution/sum(mammalia_status_distribution)
    starting_status_all_species = np.random.choice(states, size_ref_group, p = probs)
    
    # simulate IUCN history
    years = 20
    print('True rate matrix:')
    print(true_rates)
    transition_count = np.zeros([len(states),len(states)]).astype(int)
    years_in_each_category = {'LC':0,'NT':0,'VU':0,'EN':0,'CR':0,'DD':0}
    for species in np.arange(size_ref_group):
        status_series = []
        time_series = []
        starting_state = starting_status_all_species[species]
        status_index = states.index(starting_state)
        status_series.append(starting_state)
        t = 0
        while t <= years:
            status_index = states.index(starting_state)
            transition_rates = true_rates[status_index]
            #print(starting_state,transition_rates)
            waiting_time_until_next_event = np.random.exponential(1/sum(transition_rates))
            t += waiting_time_until_next_event
            if t <= years:
                probabilities = transition_rates/sum(transition_rates)
                new_state = np.random.choice(states, p=probabilities)
                new_status_index = states.index(new_state)
                transition_count[status_index][new_status_index] += 1 
                time_series.append(int(np.round(waiting_time_until_next_event)))
                starting_state = new_state
                status_series.append(new_state)
        time_series.append(years-sum(time_series))
        for i,status in enumerate(status_series):
            years_in_each_category[status] += time_series[i]
        #print(status_series)
        #print(time_series)
    #print('Simulated status transition occurrences:')
    #print(transition_count)
    print('Years spent in each category:')
    print(years_in_each_category)
    
    
    
    #_____________________Estimate status transition rates_________________________
    status_change_coutn_df = pd.DataFrame(data=transition_count,index = ['LC','NT','VU','EN','CR','DD'],columns=['LC','NT','VU','EN','CR','DD'])
    print(status_change_coutn_df)
    
    all_estimated_rates = []
    for status_a in status_change_coutn_df.columns:
        status_a_index = states.index(status_a)
        row = status_change_coutn_df.loc[status_a]
        for status_b in row.index.values:
            status_b_index = states.index(status_b)
            if not status_a == status_b:
                count = row[status_b]
                tot_time = years_in_each_category[status_a]
                rates = sample_rate_mcmc(count, tot_time, n_gen = n_gen, burnin = burnin)
                true_rate = true_rates[status_a_index,status_b_index]
                all_estimated_rates.append([true_rate] + rates)
    all_estimated_rates = np.array(all_estimated_rates)
    output_df = pd.DataFrame(all_estimated_rates,columns=['simulated_rates','estimated_rates_mean','estimated_rates_lower','estimated_rates_upper'])
    output_df.to_csv(os.path.join(outdir,'simulated_transition_rates_%s.txt'%sim_scenario_string),sep='\t',index=False)
    

# ___________________________________Plotting__________________________________
output_df = pd.read_csv(os.path.join(outdir,'simulated_transition_rates_%s.txt'%sim_scenario_string),sep='\t')
simulated = output_df.values
true_rates =  simulated[:,0]
simulated_mean = simulated[:,1]
simulated_delta_lower = simulated_mean-simulated[:,2]
simulated_delta_upper = simulated[:,3]-simulated_mean

fig = plt.figure(figsize=[4,4])
plt.errorbar(true_rates,simulated_mean, yerr=np.array([simulated_delta_lower,simulated_delta_upper]), fmt='.',ecolor='grey',elinewidth=0.1,capsize=0.)
plt.plot(true_rates,true_rates,'-',color='red',zorder=10)
#plt.axhline(0.00009729,linestyle='--',color='black',alpha=0.5)
plt.xlim(min(true_rates)-0.4*min(true_rates),max(true_rates)+0.4*max(true_rates))
plt.ylim(min(true_rates)-0.4*min(true_rates),max(true_rates)+0.4*max(true_rates))
ax=plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
plt.xlabel('True transition rates')
plt.ylabel('Estimated transition rates')
plt.tight_layout()
fig.savefig(os.path.join(outdir,'transition_rate_estimates_%s.pdf'%sim_scenario_string),bbox_inches='tight', dpi = 500)















# __________________________CODE SCRAPS________________________________________

#def sample_rate(count, tot_time, n_samples = 1, range_factor = 100, n_bins = 10000):
#    def get_loglik(count, dT, rate):
#        return np.log(rate)*count - dT*rate
#    mle = count/tot_time
#    if count == 0:
#        minRange = 1*np.exp(-10)
#        maxRange = range_factor/tot_time
#    else:
#        minRange = mle/range_factor
#        maxRange = mle*range_factor
#    rates = np.linspace(minRange, maxRange,n_bins)
#    lik = get_loglik(count,tot_time,rates)
#    lik = lik - np.max(lik)
#    sample_rates = np.random.choice(rates,size=n_samples,p=np.exp(lik)/np.sum(np.exp(lik)),replace=1)
#    #plt.plot(np.log(rates),lik)
#    #np.log(mle)
#    return sample_rates


## test the function
#estimated_rate_matrix = []
#for i,row in enumerate(transition_count):
#    starting_cat = states[i]
#    time_in_starting_cat = years_in_each_category[starting_cat]
#    sampled_rates = []
#    for j,count in enumerate(row):
#        #sampled_rate = sample_rate(count,time_in_starting_cat, n_samples = 1, range_factor = 100, n_bins = 10000)
#        sampled_rate_mcmc = sample_rate_mcmc(count,time_in_starting_cat,n_samples=1)
#        #print(states[i],states[j],np.round(sampled_rate,4),np.round(np.array([sampled_rate_mcmc]),4))
#        sampled_rates.append(sampled_rate_mcmc[0])
#    estimated_rate_matrix.append(np.array(sampled_rates))
#print('Estimated rates:')
#estimated_rate_matrix = np.array(estimated_rate_matrix)
#print(estimated_rate_matrix)



#
## get current status for all species we want to simulate___________________
#mammalia_status_distribution = np.array([3308,348,531,500,203,861])
#probs = mammalia_status_distribution/sum(mammalia_status_distribution)
#current_status_list = np.random.choice(states, sim_species,p=probs)
#final_df_current_status = pd.DataFrame(data=np.array([sim_species_list,current_status_list]).T,columns=['species','current_status'])
##final_df_current_status.to_csv(os.path.join(iucn_outdir,'current_status_all_species.txt'),sep='\t',index=False)
#
#en_risks = np.array([[p_e_year(20,0.2)]*sim_reps]*len(sim_species_list))
#
#if en_risks.shape[1] == 1:
#    en_risks = np.array([en_risks for i in range(sim_reps)])[:,:,0].T        
#en_risks_df = pd.DataFrame(np.zeros((len(sim_species_list),sim_reps+1)))
#en_risks_df.columns = ['species']+ ['en_extinction_risk_yearly_%i'%i for i in np.arange(0,sim_reps)]
#en_risks_df.species = sim_species_list
#en_risks_df.iloc[:,1:] = en_risks
##en_risks_df.to_csv(os.path.join(outdir,'en_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
#
#cr_risks = np.array([[p_e_year(10,0.5)]*sim_reps]*len(sim_species_list))
#if cr_risks.shape[1] == 1:
#    cr_risks = np.array([cr_risks for i in range(sim_reps)])[:,:,0].T
#cr_risks_df = pd.DataFrame(np.zeros((len(sim_species_list),sim_reps+1)))
#cr_risks_df.columns = ['species']+ ['cr_extinction_risk_yearly_%i'%i for i in np.arange(0,sim_reps)]
#cr_risks_df.species = sim_species_list
#cr_risks_df.iloc[:,1:] = cr_risks
##cr_risks_df.to_csv(os.path.join(outdir,'cr_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
##______________________________________________________________________________
#
#
#
##________________________________________RUN SIM_______________________________
#n_sim = 10000
#n_years = 100
#all_lc=False
#status_change=True
#plot_posterior=False
#outdir='/Users/tobias/Desktop/runsim_test/10ksim_100years'
#
#species_list_status = final_df_current_status
#transition_rates = pd.DataFrame(data = sampled_rates_df.iloc[:,1:].values,index = sampled_rates_df.status_change,columns = sampled_rates_df.columns[1:])
#species_list = en_risks_df.iloc[:,0].values
#en_ext_data = en_risks_df.iloc[:,1:].values
#cr_ext_data = cr_risks_df.iloc[:,1:].values
#
## calculate all q-matrices (all species and sim replicates)________________
#n_rates = transition_rates.shape[1]
#print("\nCalculating species-specific q-matrices ...")
#qmatrix_dict_list = []
#for i in np.arange(n_rates):
#    rates = transition_rates.iloc[:,i]
#    sys.stdout.write('\rProgress: %i %%'%int(((i+1)/n_rates)*100))
#    en_risks_rep = en_ext_data.T[i]
#    cr_risks_rep = cr_ext_data.T[i]
#    q_matrix_dict = {}
#    for j,species in enumerate(species_list):
#        en_risk = en_risks_rep[j]
#        cr_risk = cr_risks_rep[j]
#        status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category
#        q_matrix = cust_func.qmatrix(rates, status_specific_p_e)
#        q_matrix_dict[species] = q_matrix
#    qmatrix_dict_list.append(q_matrix_dict)
#
## get transition rates for DD______________________________________________
#dd_changes = []
#dd_rates = []
#for row_id,change_type in enumerate(transition_rates.index.values):
#    states = change_type.split('->')
#    if states[0] == 'DD':
#        dd_changes.append('-'.join(states))
#        rates = transition_rates[transition_rates.index==change_type].values
#        dd_rates.append(rates[0])
#dd_probs = dd_rates/sum(np.array(dd_rates))
#
#if n_sim == 0:
#    n_rep = n_rates
#else:
#    n_rep = n_sim
#
## simulations______________________________________________________________    
## add dd frequencies for additional simulation replicates
#if n_rep-n_rates >= 0:
#    resampling_rates_indexes = np.random.choice(np.arange(n_rates),(n_rep-n_rates))
#    append_this_dd = np.array([i[resampling_rates_indexes] for i in dd_probs])
#    final_dd_probs = np.concatenate([dd_probs,append_this_dd],axis=1)
#    # redraw n samples of transition rates to fill all simulation replicates
#    append_this = np.array(qmatrix_dict_list)[resampling_rates_indexes]
#    final_qmatrix_dict_list = list(qmatrix_dict_list) + list(append_this)
#else:
#    resampling_rates_indexes = np.random.choice(np.arange(n_rates),n_rep)
#    final_dd_probs = np.array([i[resampling_rates_indexes] for i in dd_probs])
#    final_qmatrix_dict_list = np.array(qmatrix_dict_list)[resampling_rates_indexes]
#
#delta_t = n_years
#dynamic_qmatrix=True
#print('\nStarting simulations ...')
#diversity_through_time,te_array,status_through_time = cust_func.run_multi_sim(n_rep,delta_t,species_list_status,final_dd_probs,final_qmatrix_dict_list,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)
#
#
#
#
## calculate some extinction stats__________________________________________
#sim_species_list = te_array[:,0].copy()
#ext_date_data = te_array[:,1:].copy()
#extinction_occs = np.array([len(row[~np.isnan(list(row))]) for row in ext_date_data])
#extinction_prob = extinction_occs/ext_date_data.shape[1]
## estimate extinction rates scaled by year
#print('\nEstimating extinction rates from simulation output...')
##ext_date_data = ext_date_data[:10,:]
#if plot_posterior:
#    with PdfPages(os.path.join(outdir,'posterior_ext_rate_histograms.pdf')) as pdf:
#        sampled_rates = np.array([get_rate_estimate(species_values,n_years,i,sim_species_list,plot_posterior=plot_posterior,pdf=pdf,n_gen=n_gen,burnin=burnin) for i,species_values in enumerate(ext_date_data)])
#else:
#    sampled_rates = np.array([get_rate_estimate(species_values,n_years,i,sim_species_list,plot_posterior=plot_posterior,pdf=0,n_gen=n_gen,burnin=burnin) for i,species_values in enumerate(ext_date_data)])
## export extinction stats to file
#column_names = ['species','rate_e_mean','rate_e_lower','rate_e_upper','simulated_p_e_in_%i_years'%delta_t]
#extinction_prob_df = pd.DataFrame(np.array([sim_species_list,sampled_rates[:,0],sampled_rates[:,1],sampled_rates[:,2],extinction_prob]).T,columns=column_names)
#extinction_prob_df[column_names[1:]] = extinction_prob_df[column_names[1:]].astype(float)
#extinction_prob_df.to_csv(os.path.join(outdir,'extinction_prob_all_species.txt'),sep='\t',index=False,float_format='%.5f')
#print('\n')
#
#
#
## plot diversity trajectory of species list________________________________
##colors = ["#9a002e","#df4a3d","#fecd5f","#5cd368","#916200"]
#colors = ["#b80033","#bf7400","#008349"]
## define time axis
#time_axis = np.array(range(len(diversity_through_time[0])))
#fig = plt.figure()
#plt.plot(time_axis,np.mean(diversity_through_time, axis =0),color=colors[0], label='accounting for GL')
#for rep in diversity_through_time:
#    plt.plot(time_axis,rep,color=colors[0],alpha=.08)  
##plt.legend()
#plt.ylabel('Total diversity')
#plt.xlabel('Years from present')
#ax = plt.gca()
#ax1 = ax.twinx()
## Set the limits of the new axis from the original axis limits
#ax1.set_ylim(ax.get_ylim())
#current_diversity = diversity_through_time[0,0]
#plt.yticks([np.mean(diversity_through_time[:,-1])],[int(current_diversity-np.mean(diversity_through_time[:,-1]))])
##plt.xticks(modified_q_matrix.year[::10],modified_q_matrix.year[::10])
#plt.ylabel('Lost species')
#plt.tight_layout()
#fig.savefig(os.path.join(outdir,'future_diversity_trajectory.pdf'),bbox_inches='tight', dpi = 500)
#
#
#
#
##plt.plot(extinction_prob_df.rate_e_mean,extinction_prob_df.rate_e_mean)



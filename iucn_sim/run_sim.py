#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate q-matrices and run simulations

Created on Wed Oct 30 20:59:28 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
np.random.seed(1234)
import os
import sys
import iucn_sim.functions as cust_func


def add_arguments(parser):
    parser.add_argument(
        '--indir',
        required=True,
        help="Path to directory created by get_rates function."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--n_years',
        default=100,
        help="How many years to simulate into the future."
    )
    parser.add_argument(
        '--n_sim',
        default=0,
        help="How many simulation replicates to run. By default (value 0) as many simulation replicates are being produced as there are available rate estimates, resulting from the get_rates function (set by --n_rep flag in get_rates). If the number of simulation replicates exceeds the number of available transition rate estimates, these rates will be randomely resampled for the remaining simulations."
    )
    parser.add_argument(
        '--status_change',
        default=1,
        help="Model IUCN status changes in future simulations. 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--conservation_increase_factor',
        default=1,
        help="The transition rates leading to improvements in IUCN conservation status are multiplied by this factor."
    )
    parser.add_argument(
        '--threat_increase_factor',
        default=1,
        help="Opposite of conservation_increase_factor, multiplies the transition rates leading to worsening in IUCN conservation status."
    )
    parser.add_argument(
        '--model_unknown_as_lc',
        default=0,
        help="Model new status for all DD and NE species as LC (best case scenario). 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--extinction_rates',
        default=1,
        help="Estimation of extinction rates from simulation results: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--n_gen',
        default=100000,
        help="Number of generations for MCMC for extinction rate estimation (default=100000)."
    )
    parser.add_argument(
        '--burnin',
        default=1000,
        help="Burn-in for MCMC for extinction rate estimation (default=1000)."
    )
    parser.add_argument(
        '--plot_diversity_trajectory',
        default=1,
        help="Plots the simulated diversity trajectory: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--plot_histograms',
        default=0,
        help="Plots histograms of simulated extinction times for each species: 0=off, 1=on (default=0)."
    )
    parser.add_argument(
        '--plot_posterior',
        default=0,
        help="Plots histograms of posterior rate estimates for each species: 0=off, 1=on (default=0)."
    )

def p_e_year(years,p_e):
    pe_year = 1-(1-p_e)**(1/years)
    return pe_year

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

## test rate estimator
#true_rate = 0.01
#n_sim = 100
#ext_time_array = (np.random.exponential(1./true_rate, n_sim)).astype(int)
#max_t=100
#ext_time_array[ext_time_array>max_t] = max_t
#get_rate_estimate(ext_time_array,max_t)


def main(args):
    
    indir = args.indir
    outdir = args.outdir
    n_years = int(args.n_years)
    n_sim = int(args.n_sim)
    allow_status_change = int(args.status_change)
    conservation_increase_factor = int(args.conservation_increase_factor)
    threat_increase_factor = int(args.threat_increase_factor)
    extinction_rates = int(args.extinction_rates)
    n_gen = int(args.n_gen)
    burnin = int(args.burnin)
    plot_diversity_trajectory = int(args.plot_diversity_trajectory)
    plot_histograms = int(args.plot_histograms)
    plot_posterior = int(args.plot_posterior)
    model_unknown_as_lc = int(args.model_unknown_as_lc)

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    species_list_status_file = os.path.join(indir,'iucn_data/current_status_all_species.txt')
    transition_rates_file = os.path.join(indir,'sampled_status_change_rates.txt')
    en_ext_risk_file = os.path.join(indir,'en_extinction_risks_all_species.txt')
    cr_ext_risk_file = os.path.join(indir,'cr_extinction_risks_all_species.txt')
    species_list_status = pd.read_csv(species_list_status_file,sep='\t')
    transition_rates = pd.read_csv(transition_rates_file,sep='\t',index_col='status_change')
    species_list = pd.read_csv(en_ext_risk_file,sep='\t').iloc[:,0].values
    en_ext_data = pd.read_csv(en_ext_risk_file,sep='\t').iloc[:,1:].values
    cr_ext_data = pd.read_csv(cr_ext_risk_file,sep='\t').iloc[:,1:].values
    
    
    # calculate all q-matrices (all species and sim replicates)________________
    n_rates = transition_rates.shape[1]
    print("\nCalculating species-specific q-matrices ...")
    qmatrix_dict_list = []
    for i in np.arange(n_rates):
        rates = transition_rates.iloc[:,i]
        sys.stdout.write('\rProgress: %i %%'%int(((i+1)/n_rates)*100))
        en_risks_rep = en_ext_data.T[i]
        cr_risks_rep = cr_ext_data.T[i]
        q_matrix_dict = {}
        for j,species in enumerate(species_list):
            en_risk = en_risks_rep[j]
            cr_risk = cr_risks_rep[j]
            status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category
            q_matrix = cust_func.qmatrix(rates, status_specific_p_e)
            if conservation_increase_factor != 1:
                indeces_lower_triangle = np.tril_indices(q_matrix.shape[0],-1)
                q_matrix[indeces_lower_triangle] = q_matrix[indeces_lower_triangle] * conservation_increase_factor
                np.fill_diagonal(q_matrix,0)
                np.fill_diagonal(q_matrix, -np.sum(q_matrix,axis=1))
            if threat_increase_factor != 1:
                indeces_upper_triangle = np.triu_indices(q_matrix.shape[0],1)
                q_matrix[indeces_upper_triangle] = q_matrix[indeces_upper_triangle] * threat_increase_factor
                np.fill_diagonal(q_matrix,0)
                np.fill_diagonal(q_matrix, -np.sum(q_matrix,axis=1))                
            q_matrix_dict[species] = q_matrix
        qmatrix_dict_list.append(q_matrix_dict)


    # get transition rates for DD______________________________________________
    dd_changes = []
    dd_rates = []
    for row_id,change_type in enumerate(transition_rates.index.values):
        states = change_type.split('->')
        if states[0] == 'DD':
            dd_changes.append('-'.join(states))
            rates = transition_rates[transition_rates.index==change_type].values
            dd_rates.append(rates[0])
    dd_probs = dd_rates/sum(np.array(dd_rates))

    if n_sim == 0:
        n_rep = n_rates
    else:
        n_rep = n_sim

    # simulations______________________________________________________________    
    # add dd frequencies for additional simulation replicates
    if n_rep-n_rates >= 0:
        resampling_rates_indexes = np.random.choice(np.arange(n_rates),(n_rep-n_rates))
        append_this_dd = np.array([i[resampling_rates_indexes] for i in dd_probs])
        final_dd_probs = np.concatenate([dd_probs,append_this_dd],axis=1)
        # redraw n samples of transition rates to fill all simulation replicates
        append_this = np.array(qmatrix_dict_list)[resampling_rates_indexes]
        final_qmatrix_dict_list = list(qmatrix_dict_list) + list(append_this)
    else:
        resampling_rates_indexes = np.random.choice(np.arange(n_rates),n_rep)
        final_dd_probs = np.array([i[resampling_rates_indexes] for i in dd_probs])
        final_qmatrix_dict_list = np.array(qmatrix_dict_list)[resampling_rates_indexes]

    #current_year = datetime.datetime.now().year 
    #final_year = current_year+int(n_years)
    delta_t = n_years
    if model_unknown_as_lc:
        print('\nSetting all DD and NE species to LC.')
        all_lc=True
    else:
        all_lc=False
    if allow_status_change:
        status_change=True
    else:
        print('\nNot simulating future status changes!')
        status_change=False
    dynamic_qmatrix=True
    print('\nStarting simulations ...')
    diversity_through_time,te_array,status_through_time = cust_func.run_multi_sim(n_rep,delta_t,species_list_status,final_dd_probs,final_qmatrix_dict_list,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)


    if plot_diversity_trajectory:
        # plot diversity trajectory of species list________________________________
        #colors = ["#9a002e","#df4a3d","#fecd5f","#5cd368","#916200"]
        colors = ["#b80033","#bf7400","#008349"]
        # define time axis
        time_axis = np.array(range(len(diversity_through_time[0])))
        fig = plt.figure()
        plt.plot(time_axis,np.mean(diversity_through_time, axis =0),color=colors[0], label='accounting for GL')
        for rep in diversity_through_time:
            plt.plot(time_axis,rep,color=colors[0],alpha=.08)  
        #plt.legend()
        plt.ylabel('Total diversity')
        plt.xlabel('Years from present')
        ax = plt.gca()
        ax1 = ax.twinx()
        # Set the limits of the new axis from the original axis limits
        ax1.set_ylim(ax.get_ylim())
        current_diversity = diversity_through_time[0,0]
        plt.yticks([np.mean(diversity_through_time[:,-1])],[int(current_diversity-np.mean(diversity_through_time[:,-1]))])
        #plt.xticks(modified_q_matrix.year[::10],modified_q_matrix.year[::10])
        plt.ylabel('Lost species')
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,'future_diversity_trajectory.pdf'),bbox_inches='tight', dpi = 500)

    
    if extinction_rates:
        # calculate some extinction stats__________________________________________
        sim_species_list = te_array[:,0].copy()
        ext_date_data = te_array[:,1:].copy()
        extinction_occs = np.array([len(row[~np.isnan(list(row))]) for row in ext_date_data])
        extinction_prob = extinction_occs/ext_date_data.shape[1]
        # estimate extinction rates scaled by year
        print('\nEstimating extinction rates from simulation output...')
        #ext_date_data = ext_date_data[:10,:]
        if plot_posterior:
            with PdfPages(os.path.join(outdir,'posterior_ext_rate_histograms.pdf')) as pdf:
                sampled_rates = np.array([get_rate_estimate(species_values,n_years,i,sim_species_list,plot_posterior=plot_posterior,pdf=pdf,n_gen=n_gen,burnin=burnin) for i,species_values in enumerate(ext_date_data)])
        else:
            sampled_rates = np.array([get_rate_estimate(species_values,n_years,i,sim_species_list,plot_posterior=plot_posterior,pdf=0,n_gen=n_gen,burnin=burnin) for i,species_values in enumerate(ext_date_data)])
        # export extinction stats to file
        column_names = ['species','rate_e_mean','rate_e_lower','rate_e_upper','simulated_p_e_in_%i_years'%delta_t]
        extinction_prob_df = pd.DataFrame(np.array([sim_species_list,sampled_rates[:,0],sampled_rates[:,1],sampled_rates[:,2],extinction_prob]).T,columns=column_names)
    else:
        column_names = ['species','simulated_p_e_in_%i_years'%delta_t]
        extinction_prob_df = pd.DataFrame(np.array([sim_species_list,extinction_prob]).T,columns=column_names)        
    extinction_prob_df[column_names[1:]] = extinction_prob_df[column_names[1:]].astype(float)
    extinction_prob_df.to_csv(os.path.join(outdir,'extinction_prob_all_species.txt'),sep='\t',index=False,float_format='%.5f')
    print('\n')    


    if plot_histograms:
        # plot histograms of extinction times
        with PdfPages(os.path.join(outdir,'extinction_time_histograms.pdf')) as pdf:
            for i,species in enumerate(te_array[:,0]):
                sys.stdout.write('\rPlotting extinction histogram for species %i/%i'%(i+1,len(te_array[:,0])))
                plt.figure()
                species_te_array = te_array[:,1:][i]
                not_na_values = species_te_array[~np.isnan(list(species_te_array))]
                heights, bins = np.histogram(not_na_values,np.arange(0,delta_t+10,10))
                percent = heights/n_rep
                plt.bar(bins[:-1],percent,width=10, align="edge")
                plt.ylim(0,0.5)
                #survival_prob = 1-sum(percent)
                #if survival_prob >= 0.5:
                #    text_color = 'green'
                #else:
                #    text_color = 'red'
                ax = plt.gca()
                #plt.text(0.05, 0.7, 'survival probability: %.2f'%survival_prob,color=text_color, horizontalalignment='left',verticalalignment='baseline', transform=ax.transAxes)
                # annotate last bar            
                if ax.patches[-1].get_height() > 0:
                    ax.text(ax.patches[-1].get_x()+3, np.round(ax.patches[-1].get_height()+0.001,4), '**', fontsize=12, color='black')
                plt.title('%s - Extinct in %i years: %i/%i'%(species,delta_t,sum(heights),n_rep))
                plt.xlabel('Years from present')
                plt.ylabel('Fraction of simulations')
                plt.tight_layout()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
    print('\n')

        

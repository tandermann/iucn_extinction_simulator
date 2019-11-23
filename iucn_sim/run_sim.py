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
import os, glob
import datetime
import sys
import iucn_sim.functions as cust_func

def p_e_year(years,p_e):
    pe_year = 1-(1-p_e)**(1/years)
    return pe_year

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
        '--plot_diversity_trajectory',
        default=1,
        help="0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--plot_histograms',
        default=0,
        help="0=off, 1=on (default=0)."
    )


def get_rate_estimate(ext_time_array,max_t,index,final_index,n_bins = 10000,n_samples = 100):
    sys.stdout.write('\rProcessing species: %i/%i '%(index,final_index))
    ext_time_array_new = ext_time_array.copy()
    ext_time_array_new[ext_time_array_new!=ext_time_array_new] = max_t
    ext_time_array_new = ext_time_array_new.astype(float)
    q_samples = np.random.random(n_bins)
    # extinction occurrences
    n_extinct = len(ext_time_array_new[ext_time_array_new< float(max_t)])
    n_alive = len(ext_time_array_new[ext_time_array_new== float(max_t)])
    p = np.zeros([ext_time_array_new.shape[0],n_bins])
    if n_extinct>0:
        p[ext_time_array_new< float(max_t),:] = q_samples*np.exp(-np.array([q_samples for i in np.arange(n_extinct)]).T*ext_time_array_new[ext_time_array_new< float(max_t)]).T
    # non extinct occurrences
    p[ext_time_array_new==float(max_t),:] = np.exp(-np.array([q_samples for i in np.arange(n_alive)]).T*ext_time_array_new[ext_time_array_new==float(max_t)]).T
    #p = np.array([q_samples*np.exp(-q_samples*sim_value) if sim_value < float(max_t) else np.exp(-q_samples*sim_value) for sim_value in ext_time_array_new])
    p = (p.T/np.sum(p,axis=1))
    all_sampled_rates = [np.random.choice(q_samples,size=n_samples,p=p_rep,replace=1) for p_rep in p.T]
    all_rates = [item for sublist in all_sampled_rates for item in sublist]
    mean_value = np.mean(all_rates)
    lower,upper = cust_func.calcHPD(all_rates,0.95)
    return [mean_value,lower,upper]

#def get_rate_estimate(ext_time_array,max_t,index,final_index,n_bins = 10000,n_samples = 100):
#    sys.stdout.write('\rProcessing species: %i/%i '%(index,final_index))
#    ext_time_array_new = ext_time_array.copy()
#    ext_time_array_new[ext_time_array_new!=ext_time_array_new] = max_t
#    ext_time_array_new = ext_time_array_new.astype(float)
#    q_samples = np.random.random(n_bins)
#    p = np.array([q_samples*np.exp(-q_samples*sim_value) if sim_value < float(max_t) else np.exp(-q_samples*sim_value) for sim_value in ext_time_array_new])
#    p = (p.T/np.sum(p,axis=1))
#    all_sampled_rates = [np.random.choice(q_samples,size=n_samples,p=p_rep,replace=1) for p_rep in p.T]
#    all_rates = [item for sublist in all_sampled_rates for item in sublist]
#    mean_value = np.mean(all_rates)
#    lower,upper = cust_func.calcHPD(all_rates,0.95)
#    return [mean_value,lower,upper]

def main(args):
    indir = args.indir
    outdir = args.outdir
    n_years = args.n_years
    plot_diversity_trajectory = args.plot_diversity_trajectory
    plot_histograms = args.plot_histograms
    n_sim = int(args.n_sim)

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
    print("\nCalculating species-specific q-matrices for all %i simulation replicates ..."%n_rates)
    qmatrix_dict_list = []
    for i in np.arange(n_rates):
        rates = transition_rates.iloc[:,i]
        sys.stdout.write('\rCalculating q-matrices for rate-replicate %i/%i '%(i,n_rates))
        en_risks_rep = en_ext_data.T[i]
        cr_risks_rep = cr_ext_data.T[i]
        q_matrix_dict = {}
        for j,species in enumerate(species_list):
            en_risk = en_risks_rep[j]
            cr_risk = cr_risks_rep[j]
            status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category
            q_matrix = cust_func.qmatrix(rates, status_specific_p_e)
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
    resampling_rates_indexes = np.random.choice(np.arange(n_rates),(n_rep-n_rates))
    append_this_dd = np.array([i[resampling_rates_indexes] for i in dd_probs])
    final_dd_probs = np.concatenate([dd_probs,append_this_dd],axis=1)
    # redraw n samples of transition rates to fill all simulation replicates
    append_this = np.array(qmatrix_dict_list)[resampling_rates_indexes]
    final_qmatrix_dict_list = list(qmatrix_dict_list) + list(append_this)

    current_year = datetime.datetime.now().year 
    final_year = current_year+int(n_years)
    all_lc=False
    status_change=True
    dynamic_qmatrix=True
    print('\nStarting simulations ...')
    
    diversity_through_time,te_array,status_through_time = cust_func.run_multi_sim(n_rep,final_year,current_year,species_list_status,final_dd_probs,final_qmatrix_dict_list,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)
    # calculate some extinction stats__________________________________________
    sim_species_list = te_array[:,0].copy()
    ext_date_data = te_array[:,1:].copy()
    extinction_occs = np.array([len(row[~np.isnan(list(row))]) for row in ext_date_data])
    extinction_prob = extinction_occs/ext_date_data.shape[1]
    # estimate extinction rates scaled by year
    print('\nEstimating extinction rates from simulation output...')
    sampled_rates = np.array([get_rate_estimate(species_values,n_years,i+1,ext_date_data.shape[0]) for i,species_values in enumerate(ext_date_data)])

    # export extinction stats to file
    column_names = ['species','rate_e_mean','rate_e_lower','rate_e_upper','simulated_p_e_by_%iAD'%final_year]
    extinction_prob_df = pd.DataFrame(np.array([sim_species_list,sampled_rates[:,0],sampled_rates[:,1],sampled_rates[:,2],extinction_prob]).T,columns=column_names)
#    extinction_prob_df.iloc[:,1:] = extinction_prob_df.iloc[:,1:].astype(float).round(decimals=6)
    extinction_prob_df[column_names[1:]] = extinction_prob_df[column_names[1:]].astype(float)
    extinction_prob_df.to_csv(os.path.join(outdir,'extinction_prob_all_species.txt'),sep='\t',index=False,float_format='%.9f')
    
    if plot_diversity_trajectory:
        # plot diversity trajectory of species list________________________________
        #colors = ["#9a002e","#df4a3d","#fecd5f","#5cd368","#916200"]
        colors = ["#b80033","#bf7400","#008349"]
        # define time axis
        time_axis = np.array(range(len(diversity_through_time[0])))+current_year
        fig = plt.figure()
        plt.plot(time_axis,np.mean(diversity_through_time, axis =0),color=colors[0], label='accounting for GL')
        for rep in diversity_through_time:
            plt.plot(time_axis,rep,color=colors[0],alpha=.08)  
        #plt.legend()
        plt.ylabel('Total diversity')
        plt.xlabel('Years AD')
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

    if plot_histograms:
        # plot histograms of extinction times
        with PdfPages(os.path.join(outdir,'extinction_time_histograms.pdf')) as pdf:
            for i,species in enumerate(te_array[:,0]):
                sys.stdout.write('\rPlotting extinction histogram for species %i/%i'%(i,len(te_array[:,0])))
                plt.figure()
                species_te_array = te_array[:,1:][i]
                not_na_values = species_te_array[~np.isnan(list(species_te_array))]
                heights, bins = np.histogram(not_na_values,np.arange(0,(final_year-current_year)+10,10))
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
                plt.title('%s - Extinct by year %i: %i/%i'%(species,final_year,sum(heights),n_rep))
                plt.xlabel('Years from present')
                plt.ylabel('Fraction of simulations')
                plt.tight_layout()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()


        
        
        

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
import .functions as cust_func

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
        '--github_repo',
        required=True,
        help="Provide path to a copy of the iucn_extinction_simulator GitHub repo. This is needed to search for pre-compiled files and updated functions. Download link: https://github.com/tobiashofmann88/iucn_extinction_simulator/archive/master.zip (make sure to unzip the downloaded repo)."
    )
    parser.add_argument(
        '--n_years',
        default=100,
        help="How many years to simulate into the future."
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
    
def main(args):
#    import argparse
#    p = argparse.ArgumentParser()
#    args = p.parse_args()    
#
#    args.indir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/output'
#    args.outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/output/simulation_results'
#    args.github_repo = '/Users/tobias/GitHub/iucn_extinction_simulator/'
#    args.n_years = 81
#    args.plot_diversity_trajectory = 1
#    args.plot_histograms = 1

    indir = args.indir
    outdir = args.outdir
    github_repo = args.github_repo
    n_years = args.n_years
    plot_diversity_trajectory = args.plot_diversity_trajectory
    plot_histograms = args.plot_histograms
    
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
    n_sim = transition_rates.shape[1]
    print("Calculating species-specific q-matrices for all %i simulation replicates ..."%n_sim)
    qmatrix_dict_list = []
    for i in np.arange(n_sim):
        rates = transition_rates.iloc[:,i]
        print(i)
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


    # simulations______________________________________________________________
    n_rep = n_sim
    current_year = datetime.datetime.now().year 
    final_year = current_year+n_years
    all_lc=False
    status_change=True
    dynamic_qmatrix=True
    print('Starting simulations ...')
    diversity_through_time,te_array,status_through_time = cust_func.run_multi_sim(n_rep,final_year,current_year,species_list_status,dd_probs,qmatrix_dict_list,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)
    
    
    # export list with extinction probability by year x________________________
    sim_species_list = te_array[:,0]
    ext_date_data = te_array[:,1:]
    extinction_occs = np.array([len(row[~np.isnan(list(row))]) for row in ext_date_data])
    extinction_prob = extinction_occs/ext_date_data.shape[1]
    extinction_prob_df = pd.DataFrame(np.array([sim_species_list,extinction_prob]).T,columns=['species','p_extinct_by_%iAD'%final_year])
    extinction_prob_df.to_csv(os.path.join(outdir,'extinction_prob_all_species.txt'),sep='\t',index=False)
    
    
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

    # plot histograms of extinction times
    with PdfPages(os.path.join(outdir,'extinction_time_histograms.pdf')) as pdf:
        for i,species in enumerate(te_array[:,0]):
            print('Plotting extinction histogram for',species)
            plt.figure()
            species_te_array = te_array[:,1:][i]
            not_na_values = species_te_array[~np.isnan(list(species_te_array))]
            heights, bins = np.histogram(not_na_values,np.arange(0,(final_year-current_year)+10,10))
            percent = heights/n_rep
            plt.bar(bins[:-1],percent,width=10, align="edge")
            plt.ylim(0,0.5)
            survival_prob = 1-sum(percent)
            if survival_prob >= 0.5:
                text_color = 'green'
            else:
                text_color = 'red'
            ax = plt.gca()
            plt.text(0.05, 0.7, 'survival probability: %.2f'%survival_prob,color=text_color, horizontalalignment='left',verticalalignment='baseline', transform=ax.transAxes)
            # annotate last bar            
            if ax.patches[-1].get_height() > 0:
                ax.text(ax.patches[-1].get_x()+3, np.round(ax.patches[-1].get_height()+0.001,4), '**', fontsize=12, color='black')
            plt.title('Simulated extinction dates - %s'%species)
            plt.xlabel('Years from present')
            plt.ylabel('Fraction of simulations')
            plt.tight_layout()
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()


        
        
        

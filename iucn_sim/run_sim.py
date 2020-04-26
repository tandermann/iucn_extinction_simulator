#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run future simulations based on IUCN data and status transition rates

Created on Wed Oct 30 20:59:28 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
import iucn_sim.functions as cust_func


def add_arguments(parser):
    parser.add_argument(
        '--input_data',
        required=True,
        help="Path to 'simulation_input_data.pkl' file created by esimate_rates function."
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
        default=10000,
        help="How many simulation replicates to run. If the number of simulation replicates exceeds the number of available transition rate estimates (produced by the 'transition_rates' function), these rates will be randomely resampled for the remaining simulations."
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
        '--plot_status_trajectories',
        default=1,
        help="Plots the simulated IUCN status trajectory: 0=off, 1=on (default=0)."
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
    parser.add_argument(
        '--plot_status_piechart',
        default=1,
        help="Plots pie charts of status distribution: 0=off, 1=on (default=1)."
    )
    parser.add_argument(
        '--seed',
        default=None,
        help="Set random seed for future simulations."
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

def select_target_species(species,species_list_status,species_list,en_ext_data,cr_ext_data):
    target_species = species
    target_index = species_list_status[species_list_status.species==target_species].index.values[0]
    species_list_status = species_list_status.iloc[target_index,:]
    species_list = np.array([species_list[target_index]])
    en_ext_data = np.array([en_ext_data[target_index]])
    cr_ext_data = np.array([cr_ext_data[target_index]])
    return pd.DataFrame(species_list_status).T,species_list,en_ext_data,cr_ext_data

def get_rate_estimate_posterior(ext_time_array,max_t,index,species_list,n_gen = 100000,burnin = 1000):
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
    return post_samples

## test rate estimator
#true_rate = 0.01
#n_sim = 100
#ext_time_array = (np.random.exponential(1./true_rate, n_sim)).astype(int)
#max_t=100
#ext_time_array[ext_time_array>max_t] = max_t
#get_rate_estimate(ext_time_array,max_t)

def main(args):
    
    # get user input___________________________________________________________
    random_seed = args.seed
    try:
        random_seed = int(random_seed)
        np.random.seed(random_seed)
        print('Simulating with starting seed %i.'%random_seed)
    except:
        print('Simulating without starting seed.')
    infile = args.input_data
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
    plot_status_trajectories = int(args.plot_status_trajectories)
    plot_status_piechart = int(args.plot_status_piechart)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    input_data = cust_func.load_obj(infile)
    species_input_data, dd_probs = input_data
    species_list = np.array([i[0] for i in species_input_data])
    current_status_list = np.array([i[1] for i in species_input_data])
    q_matrix_list = [i[2] for i in species_input_data]
    n_rates = dd_probs.shape[1]
    #__________________________________________________________________________    




    # modify q-matrices, if set by user________________________________________
    final_qmatrix_list = []
    q_matrix_list_copy = np.array(q_matrix_list).copy()
    for q_matrix_list_i in q_matrix_list_copy:
        q_matrix_list_temp = []
        for q_matrix in q_matrix_list_i:
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
            q_matrix_list_temp.append(q_matrix)
        final_qmatrix_list.append(q_matrix_list_temp)
    # turn into a dict with the n qmatrices for each species
    final_qmatrix_dict = dict(zip(species_list,final_qmatrix_list))    
    #__________________________________________________________________________    
    



    # if more n_rep are set than there are q-matrices, resample________________
    if n_sim <= n_rates:
        sample_columns = np.random.choice(np.arange(n_rates),size=n_sim,replace=False)
    # since there are only as many cr and en p(ex) estimates as there are provided GL values, we may have to resample some (but make sure all are present at least once)
    else:
        sample_columns1 = np.random.choice(np.arange(n_rates),size=n_rates,replace=False)
        sample_columns2 = np.random.choice(np.arange(n_rates),size=(n_sim-n_rates),replace=True)
        sample_columns = np.concatenate([sample_columns1,sample_columns2])    
    #__________________________________________________________________________    




    # run simulations__________________________________________________________
    delta_t = n_years
    model_ne_as_dd = False
    dynamic_qmatrix = True

    if model_ne_as_dd:
        current_status_list[current_status_list=='NE'] = 'DD'    
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

    print('\nStarting simulations ...')
    diversity_through_time,te_array,status_through_time = cust_func.run_multi_sim(n_sim,delta_t,species_list,current_status_list,dd_probs,final_qmatrix_dict,sample_columns,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix)
    # summarize simulation results
    sim_species_list = te_array[:,0].copy()
    ext_date_data = te_array[:,1:].copy()
    extinction_occs = np.array([len(row[~np.isnan(list(row))]) for row in ext_date_data])
    extinction_prob = extinction_occs/ext_date_data.shape[1]
    # produce output file for status distribution through time
    mean_status_through_time = np.mean(status_through_time,axis=2)
    year = np.arange(delta_t+1).astype(int)
    status_df_data = np.round(np.vstack([year,mean_status_through_time])).astype(int)
    status_df = pd.DataFrame(data = status_df_data.T,columns=['year','LC','NT','VU','EN','CR','EX'])
    status_df.to_csv(os.path.join(outdir,'status_distribution_through_time.txt'),sep='\t',index=False)
    np.savetxt(os.path.join(outdir,'simulated_extinctions_array.txt'),status_through_time[-1],fmt='%i')
    pd.DataFrame(data=te_array).to_csv(os.path.join(outdir,'te_all_species.txt'),sep='\t',header=False,index=False)
    #__________________________________________________________________________   




    #__________________________________________________________________________       
#    if target_species:
#        posterior = get_rate_estimate_posterior(ext_date_data[0],n_years,0,sim_species_list,n_gen=n_gen,burnin=burnin)
#        np.savetxt('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_2/figure_data/posterior_samples/%s_gl_no_status_change.txt'%target_species.replace(' ','_'),posterior,fmt='%.8f')
#        print('\nPrinted posterior')
    #__________________________________________________________________________   




    #__________________________________________________________________________           
    if plot_diversity_trajectory:
        # plot diversity trajectory of species list________________________________
        #colors = ["#9a002e","#df4a3d","#fecd5f","#5cd368","#916200"]
        # define time axis
        time_axis = np.array(range(len(diversity_through_time[0])))
        fig = plt.figure()
        y_values = np.mean(diversity_through_time, axis =0)
        plt.plot(time_axis,y_values,color="#b80033", label='accounting for GL')
        # get upper and lower confidence interval boundaries
        min_hpd, max_hpd = np.array([cust_func.calcHPD(i,0.95) for i in diversity_through_time.T]).T
        mean_min_max = np.vstack([y_values,min_hpd,max_hpd])
        np.savetxt(os.path.join(outdir,'future_diversity_trajectory.txt'),mean_min_max,fmt='%.2f')
        plt.fill_between(time_axis, min_hpd, max_hpd,
                 color="#b80033", alpha=0.2)
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
    #__________________________________________________________________________   




    #__________________________________________________________________________   
    if plot_status_trajectories:
        # color palette
        colors = ["#227a00","#a5c279","#f3d248","#6956cb","#79262a","#e34349"]
        # define time axis
        time_axis = np.array(range(len(diversity_through_time[0])))
        # plot results
        def plot_mean_and_interval(div,color,label,fig):
            plt.plot(time_axis,np.mean(div,axis=0),color=color,label=label);
            min_hpd, max_hpd = np.array([cust_func.calcHPD(i,0.95) for i in div.T]).T
            plt.fill_between(time_axis, min_hpd, max_hpd, color=color, alpha=0.2);
            return fig
        fig = plt.figure(figsize=(10,10))
        plot_mean_and_interval(status_through_time[0,:,:].T,colors[0],'LC',fig)
        plot_mean_and_interval(status_through_time[1,:,:].T,colors[1],'NT',fig)
        plot_mean_and_interval(status_through_time[2,:,:].T,colors[2],'VU',fig)
        plot_mean_and_interval(status_through_time[3,:,:].T,colors[3],'EN',fig)
        plot_mean_and_interval(status_through_time[4,:,:].T,colors[4],'CR',fig)
        plot_mean_and_interval(status_through_time[5,:,:].T,colors[5],'EX',fig)
        # add title, legend and axis-labels
        plt.legend(loc='best',fancybox=True)
        plt.title('Diversity trajectory IUCN categories - status change') #10x higher conservation
        plt.ylabel('Number species in category')
        plt.xlabel('Years from present')
        ax = plt.gca()
        ax1 = ax.twinx()
        # Set the limits of the new axis from the original axis limits
        ax1.set_ylim(ax.get_ylim())
        # annotate final counts with labels
        right_ticks = [int(np.round(np.mean(status_through_time[i,-1,:]))) for i in range(status_through_time.shape[0])]
        plt.yticks(right_ticks,right_ticks)
        #plt.xticks(modified_q_matrix.year[::10],modified_q_matrix.year[::10])
        plt.tight_layout()
        fig.savefig(os.path.join(outdir,'future_status_trajectory.pdf'),bbox_inches='tight', dpi = 500)
    #__________________________________________________________________________   




    #__________________________________________________________________________       
    if plot_status_piechart:
        statuses, counts = np.unique(current_status_list,return_counts=True)
        init_status_dict = dict(zip(statuses, counts))
        init_status_dict['EX'] = 0
        iucn_status_code = {0:'LC', 1:'NT', 2:'VU', 3:'EN', 4:'CR', 5:'EX', 6:'DD'}
        status_count_list = []
        for status_id in np.arange(status_through_time.shape[0]+1):
            status = iucn_status_code[status_id]
            if status in init_status_dict.keys():
                pre_dd_modeling_count = init_status_dict[status]
            else:
                pre_dd_modeling_count = 0
            if not status == 'DD':
                present_status_count = int(np.round(np.mean(status_through_time[status_id][0])))
                final_status_count = int(np.round(np.mean(status_through_time[status_id][-1])))
            else:
                present_status_count = 0
                final_status_count = 0
            status_count_list.append([pre_dd_modeling_count,present_status_count,final_status_count])
        status_count_list = np.array(status_count_list).T
        colors = np.array(["#227a00","#a5c279","#f3d248","#6956cb","#79262a","#b80033",'black'])
        labels = np.array(['LC', 'NT', 'VU', 'EN', 'CR', 'EX', 'DD'])
        def func(pct, allvals):
            absolute = int(np.round((pct/100.*np.sum(allvals))))
            return "{:d}".format(absolute)
        fig, axs = plt.subplots(1, 3,figsize=(12,10))
        # status distribution beginning
        wedges, texts, autotexts =axs[1].pie(status_count_list[1][status_count_list[1] >0], colors= colors[status_count_list[1] >0], autopct=lambda pct: func(pct, status_count_list[1][status_count_list[1] >0]), shadow=False,textprops=dict(color="w"))
        # status distribution end
        wedges, texts, autotexts =axs[2].pie(status_count_list[2][status_count_list[2] >0], colors= colors[status_count_list[2] >0], autopct=lambda pct: func(pct, status_count_list[2][status_count_list[2] >0]), shadow=False,textprops=dict(color="w"))
        ext = wedges[-1]
        # status distribution pre-dd
        wedges, texts, autotexts =axs[0].pie(status_count_list[0][status_count_list[0] >0], colors= colors[status_count_list[0] >0], autopct=lambda pct: func(pct, status_count_list[0][status_count_list[0] >0]), shadow=False,textprops=dict(color="w"))
        axs[0].set_title('Current (including DD)')
        axs[1].set_title('Current (DD corrected)')
        axs[2].set_title('Final (%i years)'%delta_t)
        final_labels = list(labels[status_count_list[0] >0]) + ['EX']
        plt.legend(wedges+[ext], final_labels,title="IUCN status\n(N=%i sp.)"%status_count_list[2].sum(),loc="center left",bbox_to_anchor=(1, 0, 0.5, 1))
        fig.savefig(os.path.join(outdir,'status_pie_chart.pdf'),bbox_inches='tight', dpi = 500)
    #__________________________________________________________________________   




    #__________________________________________________________________________   
    if extinction_rates:
        # calculate some extinction stats
        # estimate extinction rates scaled by year
        print('\nRunning %i MCMCs to estimate species-specific extinction rates from simulation output...'%len(species_list))
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
    extinction_prob_df.to_csv(os.path.join(outdir,'extinction_prob_all_species.txt'),sep='\t',index=False,float_format='%.8f')
    print('\n')
    #__________________________________________________________________________   




    #__________________________________________________________________________   
    if plot_histograms:
        # plot histograms of extinction times
        with PdfPages(os.path.join(outdir,'extinction_time_histograms.pdf')) as pdf:
            for i,species in enumerate(te_array[:,0]):
                sys.stdout.write('\rPlotting extinction histogram for species %i/%i'%(i+1,len(te_array[:,0])))
                plt.figure()
                species_te_array = te_array[:,1:][i]
                not_na_values = species_te_array[~np.isnan(list(species_te_array))]
                heights, bins = np.histogram(not_na_values,np.arange(0,delta_t+10,10))
                percent = heights/n_sim
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
                plt.title('%s - Extinct in %i years: %i/%i'%(species,delta_t,sum(heights),n_sim))
                plt.xlabel('Years from present')
                plt.ylabel('Fraction of simulations')
                plt.tight_layout()
                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()
    print('\n')
    #__________________________________________________________________________   


        

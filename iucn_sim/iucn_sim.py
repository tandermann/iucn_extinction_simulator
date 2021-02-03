#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 10:56:17 2021

@author: Tobias Andermann (tobiasandermann88@gmail.com)
"""

import numpy as np
import pandas as pd
import math
import pickle
import csv
import os, sys
from collections import Counter
from urllib.request import urlopen
from io import StringIO
import subprocess
import datetime
import warnings
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class get_iucn_data():
    def __init__(self,
                 target_species_list,
                 outdir,
                 reference_group,
                 reference_rank = None,
                 iucn_key = None,
                 no_online_sync = False,
                 from_file = True):
        
        self._reference_group = reference_group
        self._reference_rank = reference_rank
        self._target_species_list = target_species_list
        self._outdir = outdir
        self._iucn_key = iucn_key
        self._no_online_sync = no_online_sync
        self._from_file = from_file
        self.run()

    def run(self):
        # get user input
        taxon_reference_group = self._reference_group
        reference_rank = self._reference_rank
        species_list = self._target_species_list
        load_from_file = self._from_file
        outdir = self._outdir
        iucn_key = self._iucn_key
        avoid_precompiled_iucn_data = self._no_online_sync
        
        # get iucn history of reference group
        iucn_history_file = get_iucn_history(reference_group=taxon_reference_group,
                                                reference_rank=reference_rank,
                                                iucn_key=iucn_key,
                                                avoid_precompiled_iucn_data=avoid_precompiled_iucn_data,
                                                outdir=outdir)
        self._iucn_history_file = iucn_history_file

        # get most recent status for each taxon in target species list
        extant_taxa_current_status = get_most_recent_status_target_species(species_list=species_list,
                                                                            iucn_history_file=iucn_history_file,
                                                                            iucn_key=iucn_key,
                                                                            load_from_file = load_from_file,
                                                                            outdir=outdir)
        self._species_data = extant_taxa_current_status

        # get info about possibly extinct taxa
        possibly_extinct_taxa = get_possibly_extinct_iucn_info(iucn_history_file,
                                                                outdir=outdir)
        self._possibly_extinct_taxa = possibly_extinct_taxa


class transition_rates():
    def __init__(self,
                 species_iucn_status,
                 iucn_history_file,
                 outdir,
                 extinction_probs_mode=0,
                 possibly_extinct_list=[],
                 species_specific_regression=False,
                 rate_samples=100,
                 n_gen=100000,
                 burnin=1000,
                 seed=None,
                 load_from_file=True):
        
        self._species_data = species_iucn_status
        self._iucn_history = iucn_history_file
        self._outdir = outdir
        self._extinction_probs_mode = extinction_probs_mode
        self._possibly_extinct_list = possibly_extinct_list
        self._species_specific_regression = species_specific_regression
        self._rate_samples = rate_samples
        self._n_gen = n_gen
        self._burnin = burnin
        self._seed = seed
        self._load_from_file = load_from_file

        self.run()

    def run(self):
        # get user input___________________________________________________________
        iucn_history = self._iucn_history
        outdir = self._outdir
        try:
            extinction_probs_mode = int(self._extinction_probs_mode)
        except:
            print('\nInvalid extinction_probs_mode provided. Please choose between the currenlty available options 0 or 1')
            quit()
        possibly_extinct_list = self._possibly_extinct_list
        n_rep = int(self._rate_samples)
        n_gen = int(self._n_gen)
        burnin = int(self._burnin)
    
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    
        seed = self._seed
        try:
            self._random_seed = False
            seed = int(seed)
        except:
            seed = np.random.randint(999999999)
            self._random_seed = True
        self._seed = seed
        np.random.seed(seed)
        np.savetxt(os.path.join(outdir,'starting_seed.txt'),np.array([seed]),fmt='%i')
    
        # get input data (either from file or from memory)
        if self._load_from_file:            
            species_data_input = pd.read_csv(self._species_data,sep='\t',header=None)
        else:
            species_data_input = self._species_data
            
        species_data_input = species_data_input.dropna()
        invalid_status_taxa = species_data_input[~species_data_input.iloc[:,1].isin(['LC','NT','VU','EN','CR','DD','NE'])]
        if len(invalid_status_taxa)>0:
            print('\nFound invalid IUCN statuses:',list(invalid_status_taxa[1].values),'\n\nMake sure that the second column of your --species_data input contains the current IUCN status of your target species, which must be one of the following valid extant statuses: LC, NT, VU, EN, CR, DD, NE')
            # if this effects only a minority of taxa, continue after removing these
            if len(invalid_status_taxa)/len(species_data_input) < 0.5:
                print('\nAutomatically dropping the following taxa because of invalid IUCN status information:', list(invalid_status_taxa[0].values))
                species_data_input = species_data_input[species_data_input.iloc[:,1].isin(['LC','NT','VU','EN','CR','DD','NE'])]
            else:
                quit('\nPlease fix your species_data input file. Check presence of current IUCN status information and column order.')
            
        # get the list of species
        species_list = species_data_input.iloc[:,0].values.astype(str)
        # replace underscores in species name in case they are present
        species_list = np.array([i.replace('_',' ') for i in species_list])
        # Check if all species names are binomial
        for species in species_list:
            if len(species.split(' ')) != 2:
                print('ERROR','*'*50,'\nABORTED: All provided species names provided under --species_data flag must be binomial! Found non binomial name:\n%s\n'%species,'*'*50)
                quit()
        
        # get the current IUCN status of all species    
        current_status = species_data_input.iloc[:,1].values.astype(str)        
    
        # get GL data if provided
        gl_data_available = False
        if species_data_input.shape[1] > 2:
            gl_matrix = species_data_input.iloc[:,2:].values
            gl_data_available = True
        #__________________________________________________________________________
            
        # process the IUCN history data____________________________________________
        stat_time_df = process_iucn_history(self._iucn_history)

        # if possibly_extinct_list provided, read that list and set the status for those taxa to extinct, starting at provided year
        if len(possibly_extinct_list) > 0:
            master_stat_time_df = set_taxa_as_extinct(stat_time_df,possibly_extinct_list,from_file=self._load_from_file)
        else:
            master_stat_time_df = stat_time_df.copy()

        # extract most recent valid status for each taxon
        valid_status_dict,most_recent_status_dict,status_series,taxon_series = extract_valid_statuses(master_stat_time_df)
    
        # extinction prob mode 0: remove all currently extinct taxa
        if extinction_probs_mode == 0:
            ext_indices = np.array([num for num,i in enumerate(most_recent_status_dict.keys()) if most_recent_status_dict[i] == 'EX'])
            master_stat_time_df = master_stat_time_df.drop(ext_indices)
            master_stat_time_df.index = np.arange(len(master_stat_time_df))
            # replace any occurrence of 'EX' as a past status with NaN to avoid problems with counting types of transitions (treating these assessments as invalid)
            master_stat_time_df.replace('EX', np.nan,inplace=True)
        # extinciton prob mode 1: remove only taxa that have been extinct all along, keeping those that have recorded transition to extinct within time frame
        elif extinction_probs_mode == 1:
            ext_indices = np.array([num for num,i in enumerate(master_stat_time_df.iloc[:,1:].values.astype(str)) if 'EX' in np.unique(i) and len(np.unique(i))==2])
            master_stat_time_df = master_stat_time_df.drop(ext_indices)
            master_stat_time_df.index = np.arange(len(master_stat_time_df))
            
        # write IUCN history df to file
        master_stat_time_df.to_csv(os.path.join(outdir,'formatted_iucn_history_reference_taxa.txt'),sep='\t')
    
        # extract most recent valid status for each taxon
        valid_status_dict,most_recent_status_dict,status_series,taxon_series = extract_valid_statuses(master_stat_time_df)
        # count current status distribution
        unique, counts = np.unique(status_series, return_counts=True)
        print('\nCurrent IUCN status distribution in reference group:',dict(zip(unique, counts)))
        # count how often each status change occurs
        change_type_dict = count_status_changes(master_stat_time_df,valid_status_dict)
        print('Summing up years spend in each category ...')
        years_in_each_category = get_years_spent_in_each_category(master_stat_time_df,valid_status_dict)
    
        # write the status change data to file
        final_years_count_array = np.array([list(years_in_each_category.keys()),list(years_in_each_category.values())]).T
        np.savetxt(os.path.join(outdir,'years_spent_in_each_category.txt'),final_years_count_array,fmt='%s\t%s')
        change_type_dict_array = np.array([list(change_type_dict.keys()),list(change_type_dict.values())]).T
        np.savetxt(os.path.join(outdir,'change_type_dict.txt'),change_type_dict_array,fmt='%s\t%s')   
        #__________________________________________________________________________
        # sample transition rates for all types of changes_________________________    
        if extinction_probs_mode == 0:
            status_change_coutn_df = pd.DataFrame(data=np.zeros([6,6]).astype(int),index = ['LC','NT','VU','EN','CR','DD'],columns=['LC','NT','VU','EN','CR','DD'])
        elif extinction_probs_mode == 1:
            status_change_coutn_df = pd.DataFrame(data=np.zeros([7,7]).astype(int),index = ['LC','NT','VU','EN','CR','DD','EX'],columns=['LC','NT','VU','EN','CR','DD','EX'])
            
        for status_change in change_type_dict.keys():
            states = status_change.split('->')
            original_state = states[0]
            new_state = states[1]
            count = change_type_dict[status_change]
            status_change_coutn_df.loc[original_state,new_state] = count
        status_change_coutn_df.to_csv(os.path.join(outdir,'status_change_counts.txt'),sep='\t',index=True)
        print('Counted the following transition occurrences in IUCN history of reference group:')
        print(status_change_coutn_df)
        if not self._random_seed:
            print('Running MCMC with user-set starting seed %i ...'%seed)
        else:
            print('Running MCMC with randomely generated starting seed %i ...'%seed)    
        sampled_rates_df = pd.DataFrame(columns = ['status_change']+ ['rate_%i'%i for i in np.arange(0,n_rep)])
        for status_a in status_change_coutn_df.columns:
            row = status_change_coutn_df.loc[status_a]
            for status_b in row.index.values:
                if not status_a == status_b:
                    count = row[status_b]
                    total_time = years_in_each_category[status_a]
                    rates = sample_rate_mcmc(count, total_time, n_samples = n_rep, n_gen = n_gen, burnin = burnin)
                    sampled_rates_df = sampled_rates_df.append(pd.DataFrame(data=np.matrix(['%s->%s'%(status_a,status_b)]+list(rates)),columns = ['status_change']+ ['rate_%i'%i for i in np.arange(0,n_rep)]),ignore_index=True)
        sampled_rates_df[['rate_%i'%i for i in np.arange(0,n_rep)]] = sampled_rates_df[['rate_%i'%i for i in np.arange(0,n_rep)]].apply(pd.to_numeric)
        sampled_rates_df.to_csv(os.path.join(outdir,'sampled_status_change_rates.txt'),sep='\t',index=False,float_format='%.8f')
        print('Sampled %i rates from MCMC posterior for each transition type.'%n_rep)
        #__________________________________________________________________________
        # if mode 0, calculate extinction probabilities for EN and CR with GL data_________________________
        if extinction_probs_mode == 0:
            # calculate yearly extinction risks for categories EN and CR
            if gl_data_available:
                dims = gl_matrix.shape[1]
                en_risks = []
                for gl_array in gl_matrix:
                    if dims == 1:
                        gl_array = np.array(gl_array)
                    #replace all nan values with the standard en extinction risk
                    en_risks_species = p_e_year(np.minimum(np.maximum([20]*len(gl_array),5*gl_array),100),0.2)
                    n_nan = len(en_risks_species[en_risks_species!=en_risks_species])
                    en_risks_species[en_risks_species!=en_risks_species] = [p_e_year(20,0.2)]*n_nan
                    en_risks.append(en_risks_species)
                en_risks = np.array(en_risks)
            else:
                print('Warning: No generation length (GL) data found. Extinction risks for status EN and CR are calculated without using GL data.')
                dims = 1
                en_risks = np.array([[p_e_year(20,0.2)]]*len(species_list))
            en_risks_df = make_empty_rate_df(species_list,dims,'EN')
            en_risks_df.iloc[:,1:] = en_risks
            en_risks_df.to_csv(os.path.join(outdir,'en_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
    
            if gl_data_available:
                dims = gl_matrix.shape[1]
                cr_risks = []
                for gl_array in gl_matrix:
                    if dims == 1:
                        gl_array = np.array(gl_array)
    
                    #replace all nan values with the standard en extinction risk
                    cr_risks_species = p_e_year(np.minimum(np.maximum([10]*len(gl_array),3*gl_array),100),0.5)
                    n_nan = len(cr_risks_species[cr_risks_species!=cr_risks_species])
                    cr_risks_species[cr_risks_species!=cr_risks_species] = [p_e_year(10,0.5)]*n_nan
                    cr_risks.append(cr_risks_species)
                cr_risks = np.array(cr_risks)
            else:
                dims = 1
                cr_risks = np.array([[p_e_year(10,0.5)]]*len(species_list))
            cr_risks_df = make_empty_rate_df(species_list,dims,'CR')
            cr_risks_df.iloc[:,1:] = cr_risks
            cr_risks_df.to_csv(os.path.join(outdir,'cr_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
            
        if self._species_specific_regression:
            # make regression for all other categories based on EN and CR risks
            print('Fitting species-specific regression function to determine LC, NT, and VU extinction probabilities ...')
            vu_risks_df = make_empty_rate_df(species_list,dims,'VU')
            nt_risks_df = make_empty_rate_df(species_list,dims,'NT')
            lc_risks_df = make_empty_rate_df(species_list,dims,'LC')
            for i,species in enumerate(cr_risks_df.species.values):
                en_risks = en_risks_df.iloc[i,1:].values
                cr_risks = cr_risks_df.iloc[i,1:].values
                vu_risks = []
                nt_risks = []
                lc_risks = []
                for j,_ in enumerate(en_risks):
                    en_prob = en_risks[j]
                    cr_prob = cr_risks[j]
                    x = [4.,5.]
                    y = [en_prob,cr_prob]
                    # fitting the power function to the 2 data points of each species (EN and CR risk)
                    with warnings.catch_warnings(): 
                        # this is to avoid printing the warning from curve_fit when trying to fit function to only 2 points: "OptimizeWarning: Covariance of the parameters could not be estimated"
                        warnings.filterwarnings("ignore")
                        a_b = curve_fit(power_function,x,y);
                    # extracting the values for a and b from the curve fit function
                    a = a_b[0][0]
                    b = a_b[0][1]
                    # get values for LC, NT, and VU
                    p_year_LC = power_function(1,a,b)
                    p_year_NT = power_function(2,a,b)
                    p_year_VU = power_function(3,a,b)
                    vu_risks.append(p_year_VU)
                    nt_risks.append(p_year_NT)
                    lc_risks.append(p_year_LC)
                vu_risks_df.iloc[vu_risks_df[vu_risks_df.species == species].index.values[0],1:] = np.array(vu_risks)
                nt_risks_df.iloc[nt_risks_df[nt_risks_df.species == species].index.values[0],1:] = np.array(nt_risks)
                lc_risks_df.iloc[lc_risks_df[lc_risks_df.species == species].index.values[0],1:] = np.array(lc_risks)
            vu_risks_df.to_csv(os.path.join(outdir,'vu_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
            nt_risks_df.to_csv(os.path.join(outdir,'nt_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
            lc_risks_df.to_csv(os.path.join(outdir,'lc_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
        #__________________________________________________________________________
        # populate q-matrices______________________________________________________
        print("\nPopulating species-specific q-matrices ...")
        sampled_rates_df.index = sampled_rates_df.status_change.values
        
        if extinction_probs_mode == 0:
            transition_rates = sampled_rates_df.iloc[:,1:]
            
            # randomely sample cr and en extinction probs to be used in q-matrices.
            if n_rep <= dims:
                sample_columns = np.random.choice(np.arange(dims),size=n_rep,replace=False)
            # since there are only as many cr and en p(ex) estimates as there are provided GL values, we may have to resample some (but make sure all are present at least once)
            else:
                sample_columns1 = np.random.choice(np.arange(dims),size=dims,replace=False)
                sample_columns2 = np.random.choice(np.arange(dims),size=(n_rep-dims),replace=True)
                sample_columns = np.concatenate([sample_columns1,sample_columns2])
            # get the corresponding en and cr ex-risk columns
            cr_risks_selection = cr_risks_df.iloc[:,1:].values[:,sample_columns]
            en_risks_selection = en_risks_df.iloc[:,1:].values[:,sample_columns]
            # smae for the vu, nt, and lc cats if species_specific_regression is activated
            if self._species_specific_regression:
                vu_risks_selection = vu_risks_df.iloc[:,1:].values[:,sample_columns]
                nt_risks_selection = nt_risks_df.iloc[:,1:].values[:,sample_columns]
                lc_risks_selection = lc_risks_df.iloc[:,1:].values[:,sample_columns]
            
        elif extinction_probs_mode == 1:
            target_keys = [i for i in sampled_rates_df.status_change.values if i[-2:] == 'EX']
            ex_probs = sampled_rates_df[sampled_rates_df.status_change.isin(target_keys)].iloc[:-1,1:].values.T
            transition_rates = sampled_rates_df[~sampled_rates_df.status_change.isin(target_keys)].iloc[:30,1:]
        
        for i in np.arange(n_rep):
            rates_i = transition_rates.iloc[:,i]
            sys.stdout.write('\rProgress: %i %%'%int(((i+1)/n_rep)*100))
    
            # for each rep (i), create list of q-matrices, 1 for each species
            if extinction_probs_mode == 0:
                cr_risks_rep = cr_risks_selection[:,i]
                en_risks_rep = en_risks_selection[:,i]
                if self._species_specific_regression:
                    vu_risks_rep = vu_risks_selection[:,i]
                    nt_risks_rep = nt_risks_selection[:,i]
                    lc_risks_rep = lc_risks_selection[:,i]                
                q_matrix_list_i = []
                for j,__ in enumerate(species_list):
                    en_risk = en_risks_rep[j]
                    cr_risk = cr_risks_rep[j]
                    if self._species_specific_regression:
                        lc_nt_vu = [lc_risks_rep[j],nt_risks_rep[j],vu_risks_rep[j]]
                    else:
                        lc_nt_vu = [0.000000155728,0.000041551152,0.001053050310]
                    status_specific_p_e = np.array(lc_nt_vu+[en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category    
                    q_matrix = qmatrix(rates_i, status_specific_p_e)
                    q_matrix_list_i.append([q_matrix])
            elif extinction_probs_mode == 1:
                q_matrix_list_i = []
                status_specific_p_e = ex_probs[i]
                q_matrix = qmatrix(rates_i, status_specific_p_e)
                q_matrix_list_i = []
                for spec in species_list:
                    q_matrix_list_i.append([q_matrix])
            
            q_matrix_list_i_copy = q_matrix_list_i.copy()
            if i == 0:
                qmatrix_list_dict = dict(zip(list(species_list),q_matrix_list_i_copy)).copy()
            else:
                update_dict = [qmatrix_list_dict[species].append(q_matrix_list_i_copy[i][0]) for i, species in enumerate(list(species_list))]
        print('\n')
        #__________________________________________________________________________
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
        #__________________________________________________________________________
        # Finally write all the compiled info to a pickle file_____________________
        species_specific_data = [[species,current_status[i],qmatrix_list_dict[species]]for i,species in enumerate(species_list)]
        final_output_data = [species_specific_data,dd_probs]
        simdata_outfile = os.path.join(outdir,'simulation_input_data.pkl')
        save_obj(final_output_data,simdata_outfile)
        self._simulation_input_data = final_output_data
        self._simdata_outfile = simdata_outfile


class run_sim():
    def __init__(self,
                 input_data,
                 outdir,
                 n_years=100,
                 n_sim=10000,
                 status_change=1,
                 conservation_increase_factor=1,
                 threat_increase_factor=1,
                 model_unknown_as_lc=0,
                 until_n_taxa_extinct=0,
                 extinction_rates=1,
                 n_gen=100000,
                 burnin=1000,
                 plot_diversity_trajectory=1,
                 plot_status_trajectories=1,
                 plot_histograms=0,
                 plot_posterior=0,
                 plot_status_piechart=1,
                 seed=None,
                 load_from_file = True):
        
        self._input_data = input_data
        self._outdir = outdir
        self._n_years = n_years
        self._n_sim = n_sim
        self._status_change = status_change
        self._conservation_increase_factor = conservation_increase_factor
        self._threat_increase_factor = threat_increase_factor
        self._model_unknown_as_lc = model_unknown_as_lc
        self._until_n_taxa_extinct = until_n_taxa_extinct
        self._extinction_rates = extinction_rates
        self._n_gen = n_gen
        self._burnin = burnin
        self._plot_diversity_trajectory = plot_diversity_trajectory
        self._plot_histograms = plot_histograms
        self._plot_posterior = plot_posterior
        self._plot_status_trajectories = plot_status_trajectories
        self._plot_status_piechart = plot_status_piechart
        self._seed = seed
        self._load_from_file = load_from_file
        self.run()

    def run(self):
    
        # get user input___________________________________________________________
        seed = self._seed
        try:
            random_seed = int(seed)
            print('Simulating with user-set starting seed %i.'%random_seed)
        except:
            random_seed = np.random.randint(999999999)
            print('Simulating with randomely generated starting seed %i.'%random_seed)
        np.random.seed(random_seed)
        self._seed = random_seed
        

        outdir = self._outdir
        n_years = int(self._n_years)
        n_sim = int(self._n_sim)
        n_extinct_taxa = int(self._until_n_taxa_extinct)
        allow_status_change = int(self._status_change)
        conservation_increase_factor = int(self._conservation_increase_factor)
        threat_increase_factor = int(self._threat_increase_factor)
        extinction_rates = int(self._extinction_rates)
        n_gen = int(self._n_gen)
        burnin = int(self._burnin)
        plot_diversity_trajectory = int(self._plot_diversity_trajectory)
        plot_histograms = int(self._plot_histograms)
        plot_posterior = int(self._plot_posterior)
        model_unknown_as_lc = int(self._model_unknown_as_lc)
        plot_status_trajectories = int(self._plot_status_trajectories)
        plot_status_piechart = int(self._plot_status_piechart)
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        np.savetxt(os.path.join(outdir,'starting_seed.txt'),np.array([random_seed]),fmt='%i')
            
        if self._load_from_file:
            infile = self._input_data
            input_data = load_obj(infile)
        else:
            input_data = self._input_data
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
        diversity_through_time,te_array,status_through_time,time_until_n_extinctions_list = run_multi_sim(n_sim,delta_t,species_list,current_status_list,dd_probs,final_qmatrix_dict,sample_columns,outdir,all_lc=all_lc,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix,n_extinct_taxa=n_extinct_taxa)
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
        self._status_through_time_mean = status_df
        #self._status_through_time = status_through_time
        status_df.to_csv(os.path.join(outdir,'status_distribution_through_time.txt'),sep='\t',index=False)
        self._simulated_extinctions_through_time = status_through_time[-1]
        np.savetxt(os.path.join(outdir,'simulated_extinctions_array.txt'),status_through_time[-1],fmt='%i')
        self._extinction_times = te_array
        pd.DataFrame(data=te_array).to_csv(os.path.join(outdir,'te_all_species.txt'),sep='\t',header=False,index=False)
        #__________________________________________________________________________   
    
        #__________________________________________________________________________       
    #    if target_species:
    #        posterior = get_rate_estimate_posterior(ext_date_data[0],n_years,0,sim_species_list,n_gen=n_gen,burnin=burnin)
    #        np.savetxt('/Users/tobias/GitHub/iucn_predictions/doc/figures/Figure_2/figure_data/posterior_samples/%s_gl_no_status_change.txt'%target_species.replace(' ','_'),posterior,fmt='%.8f')
    #        print('\nPrinted posterior')
        #__________________________________________________________________________   
        
        if n_extinct_taxa:        
            np.savetxt(os.path.join(outdir,'time_until_%i_extinctions.txt'%n_extinct_taxa),time_until_n_extinctions_list,fmt='%.2f')
            self._time_until_n_extinctions = time_until_n_extinctions_list
            fig = plt.figure()
            plt.hist(time_until_n_extinctions_list)
            plt.xlabel('Time in years')
            plt.ylabel('N')
            plt.title('Simulated years until %i extinctions (%i simulations)'%(n_extinct_taxa,n_sim))
            fig.savefig(os.path.join(outdir,'time_until_%i_extinctions.pdf'%n_extinct_taxa),bbox_inches='tight', dpi = 500)
    
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
            min_hpd, max_hpd = np.array([calcHPD(i,0.95) for i in diversity_through_time.T]).T
            mean_min_max = np.vstack([y_values,min_hpd,max_hpd])
            self._future_div_mean_min_max = mean_min_max
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
                min_hpd, max_hpd = np.array([calcHPD(i,0.95) for i in div.T]).T
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
        self._extinction_probs = extinction_prob_df
        
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


def get_iucn_history(reference_group=None,reference_rank=None,iucn_key=None,avoid_precompiled_iucn_data=False,outdir=''):
    # create the r-scripts to be used later on:
    write_r_scripts(outdir,script_id = 'history')

    if reference_group:
        taxon_reference_groups = reference_group.split(',')
    else:
        print('No reference group provided. Provide the name of a taxonomic group as reference_group, e.g. Mammalia')
        quit()
    if reference_rank:
        reference_ranks = reference_rank.split(',')
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    # get IUCN history_________________________________________________________  
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    precompiled_taxon_groups = []
    precompiled_taxon_group_files = []
    if not avoid_precompiled_iucn_data:
        for taxon_group in taxon_reference_groups:         
            try:
                # look for precompiled files online    
                url = 'https://raw.githubusercontent.com/tobiashofmann88/iucn_extinction_simulator/master/data/precompiled/iucn_history/%s_iucn_history.txt'%taxon_group.upper()
                urlpath =urlopen(url)
                string = urlpath.read().decode('utf-8')        
                string_input = StringIO(string)
                ref_group_data = pd.read_csv(string_input, sep="\t")
                hist_outfile = os.path.join(outdir,os.path.basename(url))
                ref_group_data.to_csv(hist_outfile,sep='\t',index=False)
                precompiled_taxon_groups.append(str.lower(taxon_group))
                precompiled_taxon_group_files.append(hist_outfile)
            except:
                pass

    iucn_history_files = []
    for i,taxon_group in enumerate(taxon_reference_groups):
        if str.lower(taxon_group) in precompiled_taxon_groups and not avoid_precompiled_iucn_data:
            print('Using precompiled IUCN history data for %s.'%taxon_group)                    
            iucn_history_files.append([file for file in precompiled_taxon_group_files if os.path.basename(file).startswith(str.upper(taxon_group))][0])
        else:
            print('Fetching IUCN history using rredlist')
            rank = reference_ranks[i]
            if rank == 'genus':
                iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_iucn_status_data_and_species_list.r'), str.capitalize(taxon_group), str.lower(rank), iucn_key, outdir]
            else:
                iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_iucn_status_data_and_species_list.r'), str.upper(taxon_group), str.lower(rank), iucn_key, outdir]
                
            if not iucn_key:
                print('ERROR','*'*50,'\nIUCN-KEY ERROR: Need to download IUCN history for specified reference group: %s. Please provide a valid IUCN key (using the --iucn_key flag). Alternatively choose a precompiled reference group (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/)'%(taxon_group),'*'*55)
                quit()
            #iucn_error_file = os.path.join(outdir,'get_iucn_status_data_and_species_list_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            run_iucn_cmd = subprocess.Popen(iucn_cmd)
            run_iucn_cmd.wait()
            iucn_history_files.append(os.path.join(outdir,'%s_iucn_history.txt'%str.upper(taxon_group)))
            
    if len(iucn_history_files) > 1:
        df_previous = pd.DataFrame()
        outfile_stem = ''
        for i in iucn_history_files:
            outfile_stem += os.path.basename(i).split('_')[0]+'_'
            df_new = pd.read_csv(i,sep='\t')
            if len(df_previous) > 0:
                df_previous = pd.concat([df_previous,df_new],sort=False,ignore_index=True)
            else:
                df_previous = df_new
            os.remove(i)
        df_previous.to_csv(os.path.join(outdir,'%siucn_history.txt'%str.upper(outfile_stem)),sep='\t',index=False)
    else:
        outfile_stem = os.path.basename(iucn_history_files[0]).split('_')[0]+'_'

    # save the file name of the output file for further actions that need to read this file
    iucn_history_file = os.path.join(outdir,'%siucn_history.txt'%str.upper(outfile_stem))
    return(iucn_history_file)


def process_iucn_history(iucn_history_file,iucn_start_year=2001,final_year=None):
    # get current IUCN status all target species_______________________________
    # process the IUCN history data
    if not final_year:
        current_year = datetime.datetime.now().year  
    else:
        current_year = final_year
    master_stat_time_df = pd.DataFrame(columns=['species']+list(np.arange(iucn_start_year,current_year+1).astype(str)))
    statuses_through_time = pd.read_csv(iucn_history_file, delimiter = '\t')
    target_columns = [column for column in master_stat_time_df.columns if column in statuses_through_time.columns]
    master_stat_time_df[target_columns] = statuses_through_time[target_columns]

    # check if we have sufficient number of species for rate estimation
    if len(master_stat_time_df) < 1000:
        print('\n\n%s'%('#'*50),'\nWarning: Only %i species in reference dataset. This may not be sufficient for proper estimation of status transition rates. It is recommended to choose a larger reference group encompassing >1000 species at the least!\n\nContinuing processing IUCN history data of reference group ...\n%s\n\n'%(len(master_stat_time_df),'#'*50))
    # treat EW as EX
    master_stat_time_df.replace('EW', 'EX',inplace=True)
    # replace occurrences of NR (not recognized) with nan
    master_stat_time_df.replace('NR', np.nan,inplace=True)
    
    # clean and sort df
    master_stat_time_df = master_stat_time_df.sort_values(by='species')
    master_stat_time_df = master_stat_time_df.drop_duplicates()
    master_stat_time_df.index = np.arange(len(master_stat_time_df))   
    
    # set most recent status to NE for species without any IUCN status information
    na_row_indeces = np.where(master_stat_time_df.iloc[:,1:].T.isnull().all().values)
    for index in na_row_indeces:
        master_stat_time_df.iloc[index,-1] = 'NE'
    return(master_stat_time_df)

def set_taxa_as_extinct(stat_time_df,possibly_extinct_list,from_file=True): # provide 2-column df as possibly_extinct_list with species names and supposed year of extinction
    master_stat_time_df = stat_time_df.copy()
    if from_file:
        pex_data = pd.read_csv(possibly_extinct_list,sep='\t')
    else:
        pex_data = pd.DataFrame(possibly_extinct_list)
    pex_species_list = pex_data.iloc[:,0].values.astype(str)
    pex_year = pex_data.iloc[:,1].values.astype(int)
    column_names = master_stat_time_df.columns.values
    row_names = master_stat_time_df.species.values
    #df_selection = master_stat_time_df[master_stat_time_df.species.isin(pex_species_list)]
    for i,species in enumerate(pex_species_list):
        row_index = np.where(row_names==species)[0][0]
        assessment_year = pex_year[i]
        column_index = np.where(column_names==str(assessment_year))[0][0]
        master_stat_time_df.iloc[row_index,column_index:] = 'EX'
    return(master_stat_time_df)

def get_most_recent_status_target_species(species_list=[],iucn_history_file=None,iucn_key=None,load_from_file=True,outdir=''):
    if iucn_history_file:
        # process IUCN history file
        master_stat_time_df = process_iucn_history(iucn_history_file)
        valid_status_dict_refgroup,most_recent_status_dict_refgroup,status_series_refgroup,taxon_series_refgroup = extract_valid_statuses(master_stat_time_df)
        provided_history_file = True
    else:
        provided_history_file = False
    # see if we want to only stick with the reference taxa
    if len(species_list) == 0:
        gl_data = []
        if provided_history_file:
            joined_df = pd.DataFrame([taxon_series_refgroup,status_series_refgroup]).T
        else:
            print('ERROR','*'*50,'\nABORTED: You need to provide either a target species list or a IUCN history file (or both). Found neither.')
    else:
        if not load_from_file:
            species_list_data = pd.DataFrame(species_list)
            species_list_data.columns = np.arange(len(species_list_data.columns)) # make sure column naming is integers, as in header=None scenario
        else:
            # get list of species we want to simulate
            species_list_data = pd.read_csv(species_list,sep='\t',header=None)
        species_list_data = species_list_data.drop_duplicates(0)
        species_list_data = species_list_data.sort_values(0)
        species_list_data.index = np.arange(len(species_list_data))
        species_list = np.unique(species_list_data.iloc[:,0].values.astype(str))
        species_list = np.array([i.replace('_',' ') for i in species_list])
        # Check if all species names are binomial
        for species in species_list:
            if len(species.split(' ')) != 2:
                print('ERROR','*'*50,'\nABORTED: All provided species names provided under --target_species_list flag must be binomial! Found non binomial name:\n%s\n'%species,'*'*50)
                quit()        
        try:
            gl_data = species_list_data.iloc[:,1:].values
            print('Found generation length data, which will be written to output file for downstream processing.')
        except:
            gl_data = []
        # see which species are not in reference group and hence need current IUCN status extraction
        if provided_history_file:
            remaining_species_wo_iucn_status = sorted(list(set(species_list)-set(taxon_series_refgroup)))
            fraction_overlap = np.round(1-(len(remaining_species_wo_iucn_status)/len(species_list)),3)      
            if fraction_overlap <= 0.5:
                print('\n\n%s'%('#'*50),'\nWarning: Only %.3f of target species found in reference group. Is this intended? You may want to reconsider chosing a reference group encompassing more of your target species.\n\nContinuing with downloading current status information for all remaining species...\n%s\n\n'%(fraction_overlap,'#'*50))
            else:
                print('\n\nA fraction of %.3f of the specified target species is present in reference group.\n\n'%fraction_overlap)
        else:
            remaining_species_wo_iucn_status = sorted(list(species_list))
        tmp_outdir = os.path.join(outdir,'other_files')
        if not os.path.exists(tmp_outdir):
            os.makedirs(tmp_outdir)
        if len(remaining_species_wo_iucn_status) > 0:
            # write r-script for downloading missing taxa
            write_r_scripts(outdir,script_id = 'current')
            # write the missing taxa to file for the r-script
            species_list_out_file = os.path.join(tmp_outdir,'target_species_not_in_reference_group.txt')
            np.savetxt(species_list_out_file,remaining_species_wo_iucn_status,fmt='%s')
            # extract the current status for those missing species
            print('Extracting current status for target species...')
            iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_current_iucn_status_missing_species.r'), species_list_out_file, iucn_key, tmp_outdir]
            #iucn_error_file = os.path.join(iucn_outdir,'get_current_iucn_status_missing_species_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            if not iucn_key:
                quit('***IUCN-KEY ERROR:*** Trying to download current status for target species from IUCN. Please provide a valid IUCN key (using the --iucn_key flag) to access IUCN data. Alternatively you can turn off this functionality by setting "--target_species_list 0". In that case you need to compile your own current IUCN status list manually. In that case store the status data in a tab-separated format with the header "species current_status".')
            print('Downloading current IUCN status information for %i target species that are not present in reference group:'%len(remaining_species_wo_iucn_status))
            run_iucn_cmd = subprocess.Popen(iucn_cmd)
            run_iucn_cmd.wait()
            # get the IUCN data and combine with recent status info from refgroup taxa
            current_status_target_species_file = os.path.join(tmp_outdir,'current_status_missing_species.txt')
            current_status_missing_taxa = pd.read_csv(current_status_target_species_file,sep='\t',header=None)
            # print info which species were not found in IUCN
            current_status_missing_list = current_status_missing_taxa[1].values.astype(str)
            nan_taxa = current_status_missing_taxa[0].values[current_status_missing_list=='nan']
            if len(nan_taxa) > 0:
                print('\n\nNo IUCN information found for the following %i species. This could be due to taxonomic issues. Make sure that all species names match with the most recent IUCN taxonomy.\n\nFor now, these species will be coded as Not Evaluated (NE)...\n\n%s\n\n'%(len(nan_taxa),str(nan_taxa)))        
            if provided_history_file:
                target_reference_taxa = list(set(species_list)-set(current_status_missing_taxa[0].values))
                status_series_refgroup = np.array(status_series_refgroup).astype(str)
                taxon_series_refgroup = np.array(taxon_series_refgroup).astype(str)
                status_reference_taxa = [status_series_refgroup[taxon_series_refgroup == i][0] for i in target_reference_taxa]
                current_status_reference_taxa = pd.DataFrame(data = np.array([target_reference_taxa,status_reference_taxa]).T)
                joined_df = pd.concat([current_status_missing_taxa,current_status_reference_taxa],ignore_index=True).sort_values(by=[0])
            else:
                joined_df = current_status_missing_taxa.sort_values(by=[0])
        else:
            target_reference_taxa = list(species_list)
            status_series_refgroup = np.array(status_series_refgroup).astype(str)
            taxon_series_refgroup = np.array(taxon_series_refgroup).astype(str)
            status_reference_taxa = [status_series_refgroup[taxon_series_refgroup == i][0] for i in target_reference_taxa]
            joined_df = pd.DataFrame(data = np.array([target_reference_taxa,status_reference_taxa]).T)
    # fix NE taxa and sort
    joined_df = joined_df.replace(np.nan,'NE')        
    joined_df.index = np.arange(len(joined_df))
    if len(gl_data) > 0:
        joined_df = pd.concat([joined_df,pd.DataFrame(np.round(gl_data,3))],axis=1)
    # remove all extinct taxa from future steps
    ext_boolean = joined_df.iloc[:,1]=='EX'
    extinct_taxa = joined_df.iloc[ext_boolean.values,:]
    print('The following taxa listed in your target species list are extinct according to IUCN:\n')
    print(extinct_taxa.iloc[:,0].values)
    print('\nThese taxa will be removed from the list for downstream processing.')
    extant_taxa = joined_df.iloc[~ext_boolean.values,:]
    extant_taxa.to_csv(os.path.join(outdir,'species_data.txt'),sep='\t',index=False,header=False)
    return(extant_taxa)

def get_possibly_extinct_iucn_info(iucn_history_file,outdir=''):
    # get info about possibly extinct taxa_____________________________________
    url = 'https://raw.githubusercontent.com/tobiashofmann88/iucn_extinction_simulator/master/data/precompiled/pex_taxa/2020_1_RL_Stats_Table_9.txt'
    urlpath =urlopen(url)
    string = urlpath.read().decode('utf-8')        
    string_input = StringIO(string)
    pe_data = pd.read_csv(string_input, sep="\t")
    # which species are we looking for?
    reference_taxa = pd.read_csv(iucn_history_file,sep='\t').species.values.astype(str)    
    # extract all these species and write to file
    reference_taxa_listed_as_pe = pe_data[pe_data['Scientific name'].isin(reference_taxa)]
    reference_taxa_listed_as_pe = reference_taxa_listed_as_pe.iloc[:,[0,3]]
    pe_data_outfile = os.path.join(outdir,'possibly_extinct_reference_taxa.txt')
    reference_taxa_listed_as_pe.to_csv(pe_data_outfile,sep='\t',index=False)
    return(reference_taxa_listed_as_pe)

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

def p_e_year(years,p_e):
    pe_year = 1-(1-float(p_e))**(1/years)
    return pe_year

def update_multiplier(q,d=1.1):
    u = np.random.uniform(0,1)
    l = 2*np.log(d)
    m = np.exp(l*(u-.5))
    new_q = q * m
    return new_q, np.log(m)

def sample_rate_mcmc(count, tot_time, n_samples = 1, n_gen = 100000,burnin = 1000):
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
    sampled_rates = np.random.choice(post_samples,n_samples,replace=False)
    return sampled_rates

def power_function(x,a,b):
    # defining the power function
    y = float(a)*x**float(b)
    return y

def make_empty_rate_df(species_list,rate_columns,status_label):
    rate_df = pd.DataFrame(np.zeros((len(species_list),rate_columns+1)))
    rate_df.columns = ['species']+ ['%s_p_ext_%i'%(status_label,i) for i in np.arange(0,rate_columns)]
    rate_df.species = species_list
    return rate_df

def format_iucn_df(statuses_through_time_file,start_year,current_year):
    new_df = pd.DataFrame()
    new_df['species'] = statuses_through_time_file['species']
    for number in range(start_year,current_year+1):
        year = str(number)
        if year in statuses_through_time_file.columns:
            pass
        else:
            statuses_through_time_file[year] = np.nan
        new_df[year] = statuses_through_time_file[year]
    new_df.replace({'NE': np.nan},inplace=True)
    print('Completed: Adding missing columns and selecting target years.')
    return new_df

def exclude_extinct(iucn_dataframe, status_column_name):
    new_df = iucn_dataframe.copy()
    new_df = new_df[(new_df[status_column_name] != 'EW') & (new_df[status_column_name] != 'EX')]
    # a few odd species have come back from extinction, meaning that even though their current status is neither 'EX' or 'EW' they had 'EX' as a previous status
    new_df.replace({'EX': 'NaN'},inplace=True)
    new_df.drop(new_df.columns[[-2, -1]], axis=1, inplace=True)
    return new_df

def count_status_changes(iucn_dataframe,valid_status_dict):
    change_types = {}
    for row in iucn_dataframe.iterrows():
        species_index = row[0]
        species_name = iucn_dataframe.species[species_index]
        current_status = ''
        if species_name in valid_status_dict:
            valid_status = list(valid_status_dict[species_name])
            for i in range(len(valid_status)):        
                # compile all change types and count
                current_status = valid_status[i]
                #this takes care of when we are getting out of the list index for the following steps and when we are at the end of the list anyways
                if i+1 >= len(valid_status):
                    pass
                elif valid_status[i+1] != current_status:
                    next_status = valid_status[i+1]
                    change = "%s->%s" %(current_status,next_status)
                    change_types.setdefault(change,[])
                    change_types[change].append(1)
    for category in change_types:
        change_types[category] = sum(change_types[category])
    return change_types

def get_years_spent_in_each_category(extant_iucn_df,valid_status_dict):
    species_list = extant_iucn_df.species.values
    only_year_columns = extant_iucn_df.iloc[:,1:].copy()
    max_index = only_year_columns.shape[-1]
    years_in_status = [np.diff(np.append(np.where(only_year_columns[extant_iucn_df.species==species].notna().values[0])[0],max_index)) for species in species_list]
    status_master_list = np.array([[item for sublist in [[j]*years_in_status[index_i][index_j] for index_j,j in enumerate(valid_status_dict[i])] for item in sublist]+['NA']*(max_index-sum(years_in_status[index_i])) for index_i,i in enumerate(species_list)])
    status_master_flat_array = np.concatenate(status_master_list)
    statuses, counts = np.unique(status_master_flat_array,return_counts=True)
    years_in_each_category = dict(zip(statuses, counts))
    years_in_each_category.pop('NA')
    return years_in_each_category

def replace_iucn_status_with_int(change_types,sum_years):
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    new_change_types = {}
    new_sum_years = {}
    for element in change_types:
        states = element.split('->')
        new_states = []
        new_states.append(str(iucn_code[states[0]]))
        new_states.append(str(iucn_code[states[1]]))
        new_states = '->'.join(new_states)
        new_change_types.setdefault(new_states,change_types[element])
    for status in sum_years:
        new_status = iucn_code[status]
        new_sum_years.setdefault(new_status,sum_years[status])
    return new_change_types,new_sum_years


def save_obj(obj, file_name):
    with open(file_name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)

def qmatrix(rates, status_specific_p_e):   
    Q_matrix = np.zeros((5,6))
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    for row_id,change_type in enumerate(rates.index.values):
        states = change_type.split('->')
        if not 'DD' in states:
            original_state = iucn_code[states[0]]
            derived_state = iucn_code[states[1]]
            rate = rates[row_id]
            Q_matrix[original_state,derived_state] = rate
            Q_matrix[original_state,5] = status_specific_p_e[original_state]   
    # This fills the diagonal line of the matrix with negative values, so that the line has a sum of 0
    np.fill_diagonal(Q_matrix, -np.sum(Q_matrix,axis=1))
    np.set_printoptions(precision=8)
    return Q_matrix

def treat_dd_species(iucn_dataframe,change_type_dict,all_lc = False):
    new_df = iucn_dataframe.copy()
    #get the change frequencies from DD species to a real category and based on that draw new, valid categories for DD data
    dd_change_vector = np.zeros(5)
    status_index_dict = {'LC':0,'NT':1,'VU':2,'EN':3,'CR':4}
    for key in list(change_type_dict.keys()):
        if key.startswith('DD'):
            old_status,new_status = key.split('->')
            new_status_index = status_index_dict[new_status]
            new_value = change_type_dict[key]
            dd_change_vector[new_status_index] = new_value
    dd_change_probabilities_vector = dd_change_vector/sum(dd_change_vector)
    #extant_iucn_df.replace('DD', list(status_index_dict.keys())[list(status_index_dict.values()).index(random_choice_P(dd_change_vector)[1])])
    dd_coords = np.where(new_df == 'DD')
    new_draws = np.random.choice(['LC','NT','VU','EN','CR'], size=len(dd_coords[0]), replace=True, p=dd_change_probabilities_vector)
    if all_lc:
        new_draws = np.array(['LC']*len(dd_coords[0]))
  #  iucn_dataframe.iloc[dd_coords[0],dd_coords[1]]
    new_df = np.array(new_df)
    new_df[new_df=="DD"] = new_draws
    new_df = pd.DataFrame(new_df)
    new_df.columns = np.array(list(iucn_dataframe.columns))
    return new_df

def extract_valid_statuses(formatted_status_through_time_file):
    taxon_list = list(formatted_status_through_time_file.species.values)
    #valid_status_matrix = [list(line[line.isin(['NE','EX','EW','DD','CR','EN','VU','NT','LC'])].values) for it,line in formatted_status_through_time_file.iterrows()]
    df_array = formatted_status_through_time_file.values
    valid_status_matrix = [list(line[np.isin(line,['DD','EX','CR','EN','VU','NT','LC'])]) for line in df_array]
    # if taxon has no single valid status in IUCN history, model as NE at present
    valid_status_matrix = [['NE'] if len(i) == 0 else i for i in valid_status_matrix]
    valid_status_dict = dict(zip(taxon_list, valid_status_matrix))
    current_status_list = [line[-1] for line in valid_status_matrix]
    most_recent_status_dict = dict(zip(taxon_list, current_status_list))
    return valid_status_dict,most_recent_status_dict,current_status_list,taxon_list

def random_choice_P(vector): # randomly sample an element of 'vector' based on their values, input needs to be float
    probDeath=vector/sum(vector) # (larger values have higher prob of being sampled)
    r=np.random.choice(probDeath, p=probDeath)
    ind=np.where(probDeath==r)[0][0]
    return [vector[ind], ind]

def round_up(value):
    try:
        return math.ceil(float(value))
    except ValueError:
        return value

def simulate_extinction_and_status_change(delta_t,list_of_all_current_species_statuses,species_list,outdir,qmatrix_dict,status_change=False,dynamic_qmatrix=True,n_extinct_taxa=0):
	# write the species name and the index to a separate txt file
    # this is necessary because two years get lost when summarizing the simulated data into year bins
    #final_year = final_year+2
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    species_cat_arg = [iucn_code[x] for x in list_of_all_current_species_statuses]
    species_cat = np.array(species_cat_arg)+0
    current_state=species_cat
    root = -delta_t
    n_columns = abs(root)
    status_array = np.tile(np.array([current_state]).transpose(), (1, n_columns))
    species_list_log = open("%s/species_list.txt" %outdir, "w")
    species_list_file=csv.writer(species_list_log, delimiter='\t')
    #print(q_matrix)
    extinction_time_list = []
    for taxon in range(len(species_list)):     
        species = species_list[taxon]
        if dynamic_qmatrix:
            q_matrix = qmatrix_dict[species].copy()
        else:
            q_matrix = qmatrix_dict
        #print(q_matrix)
        if not status_change:
            # if the status change is deactivated set all status shift rates in the q-matrix to 0
            q_matrix[0:5,0:5] = 0.0
            np.fill_diagonal(q_matrix, -np.sum(q_matrix,axis=1))
        species_list_file.writerow([taxon,species])
        current_state = int(status_array[taxon,0])
        D = -np.diagonal(q_matrix)
        extinction = 'extant'
    	#start at 1 so that the current status will always be protected/kept stable
        t=0
        while t < abs(root):
    		# draw randomely when the next change is going to happen, based on the negative inverse diagonal value of the q-matrix (corrected by scale)
            delta_t =np.random.exponential(1./(D[current_state]))
    		# fill all columns between t and delta_t with the current_status (gets updated after this)
            status_array[taxon,np.arange(int(t),min(int(t+delta_t), abs(root)))] = current_state
            t = min(t+delta_t, abs(root))
    		# get the q_row of the current status containing all change rates
            Qrow = q_matrix[current_state]+0
    		# eliminate the negative value in the diagonal and set to 0 since we don't want to sample this one
            Qrow[Qrow<0]= 0 #=1-sum(Qrow[Qrow>0])
            [rate, ind] = random_choice_P(Qrow)
            current_state = ind
            if ind == 5:
                # go extinct
                time_of_extinction = delta_t+t
                if time_of_extinction > abs(root):
                    break
                else:
                    status_array[taxon,np.arange(int(time_of_extinction), abs(root))] = 5
                    extinction = time_of_extinction
                    break
        extinction_time_list.append(extinction)
    a = np.array(species_list)
    b = np.array(extinction_time_list)
    time_until_n_extinctions = None
    if n_extinct_taxa:
        ext_taxa = b[b!='extant'].astype(float)
        if len(ext_taxa) < n_extinct_taxa:
            print('\n\nWARNING: The --n_extinct_taxa flag is being ignored in this analysis.\n\nCan not find %i extinct taxa in this simulation replicate. You may have to increase the simulated time frame by adjusting the --n_years value.'%n_extinct_taxa)
        else:
            # select the n smallest values from the extinction times
            sorted_times = np.sort(ext_taxa)
            time_until_n_extinctions = np.max(sorted_times[:n_extinct_taxa])
    extinction_array = np.stack((a, b))
    list_ex = list(map(round_up, extinction_array[1]))
    extinctions_per_year = Counter(list_ex)
    extinctions_per_year = dict(extinctions_per_year)
    return status_array, extinction_array, extinctions_per_year, time_until_n_extinctions
    

def get_dtt_array_from_extinction_per_year_dict(extinction_dict_sim_out,current_year,final_year):
    year_bins = np.arange(current_year,final_year+1)
    ext_count= [extinction_dict_sim_out[i] if i in extinction_dict_sim_out.keys() else 0 for i in year_bins]
    diversity_through_time = sum(extinction_dict_sim_out.values()) - np.cumsum(ext_count)
    diversity_through_time
    return diversity_through_time


def run_multi_sim(n_rep,delta_t,species_list,input_status_list,dd_probs,qmatrix_dict,rate_index_list,outdir,all_lc=False,status_change=True,dynamic_qmatrix=True,n_extinct_taxa=0):
    iucn_code = {'LC':0, 'NT':1, 'VU':2, 'EN':3, 'CR':4}
    extinct_per_year_array = np.zeros([n_rep,delta_t+1])
    te_array = np.zeros((len(species_list),n_rep+1)).astype(object)
    te_array[:,0]=species_list
    status_through_time_dict = {}
    status_through_time = np.zeros([6,delta_t+1,n_rep])
    time_until_n_extinctions_list = []
    for n in range(n_rep):
        sys.stdout.write('\rRunning simulation rep %i/%i' %(n+1,n_rep))
        target_column = rate_index_list[n]
            
        # new modeling of DD species every rep
        current_status_list_new_dd = input_status_list.copy()
        dd_indices = np.where([current_status_list_new_dd=='DD'])[1]    
        dd_prob_vector = dd_probs.T[target_column]
        if all_lc:
	        new_draws = np.array(['LC']*len(dd_indices))
        else:
            new_draws = np.random.choice(['LC','NT','VU','EN','CR'], size=len(dd_indices), replace=True, p=dd_prob_vector)
        current_status_list_new_dd[dd_indices] = new_draws
        
        # new modeling of NE species every rep
        status_count_dict = Counter(input_status_list)
        counts = np.array([status_count_dict[key] for key in status_count_dict.keys() if key not in ['DD','NE']])
        ne_probs = counts/sum(counts)
        status_array_count = [key for key in status_count_dict.keys() if key not in ['DD','NE']]        
        ne_indices = np.where([current_status_list_new_dd=='NE'])[1]
        if all_lc:
	        new_draws = np.array(['LC']*len(ne_indices))
        else:
            new_draws = np.random.choice(status_array_count, size=len(ne_indices), replace=True, p=ne_probs)
        current_status_list_new_dd[ne_indices] = new_draws

        # sample a different q-matrix for each rep, according to the rate_index_list input
        if dynamic_qmatrix:
            q_matrix_list_rep = [qmatrix_dict[species][target_column] for species in species_list]
            qmatrix_dict_rep = dict(zip(species_list,q_matrix_list_rep))
        # the dynamic_qmatrix flag can be turned off to pass only a single q-matrix that will be used for all taxa and sim reps
        else:
            qmatrix_dict_rep = qmatrix_dict
    
        future_status_array, extinction_array, extinction_dict, time_until_n_extinctions = simulate_extinction_and_status_change(delta_t,current_status_list_new_dd,species_list,outdir,qmatrix_dict_rep,status_change=status_change,dynamic_qmatrix=dynamic_qmatrix,n_extinct_taxa=n_extinct_taxa)
        time_until_n_extinctions_list.append(time_until_n_extinctions)
        # diversity_through time array
        for year in extinction_dict.keys():
            if not year == 'extant':
                pos = int(year)
                extinct_per_year_array[n,pos] = extinction_dict[year]
        # status through time array
        year = 1
        for column in future_status_array.T:
            the_dict = dict(Counter(list(column)))
            for status in range(6):
                key = '%s_%s'%(str(year),str(status))
                status_through_time_dict.setdefault(key,[])
                if status in the_dict.keys():
                    status_through_time_dict[key].append(the_dict[status])
                else:
                    status_through_time_dict[key].append(0)
            year += 1
        status_list = [iucn_code[status] for status in current_status_list_new_dd]
        counts_present = dict(Counter(status_list))
        for i in counts_present.keys():
            status_through_time[i,0,n]=counts_present[i]
        # fill te_array
        for species in range(len(extinction_array[1])):
            extinction_date = ''
            if extinction_array[1][species] == 'extant':
                extinction_date = np.nan
            else:
                extinction_date = np.round(float(extinction_array[1][species]),2)
            te_array[species,n+1]=extinction_date
    # finish diversity through time array
    diversity_through_time = sum(extinction_dict.values()) - np.cumsum(extinct_per_year_array,axis=1)
    # finish status_through_time_array
    for key in status_through_time_dict:
        x = int(key.split('_')[1])
        y = int(key.split('_')[0])
        status_through_time[x,y,:] = status_through_time_dict[key]
    # write data-objects to output folder
    #with open(os.path.join(outdir,'diversity_through_time.pkl'), 'wb') as f:
    #    pickle.dump(diversity_through_time, f, pickle.HIGHEST_PROTOCOL)
    #with open(os.path.join(outdir,'te_array.pkl'), 'wb') as f:
    #    pickle.dump(te_array, f, pickle.HIGHEST_PROTOCOL)
    #with open(os.path.join(outdir,'status_through_time.pkl'), 'wb') as f:
    #    pickle.dump(status_through_time, f, pickle.HIGHEST_PROTOCOL)
    return diversity_through_time,te_array,status_through_time,time_until_n_extinctions_list


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


def write_r_scripts(output_folder,script_id = 'both'):
	
    script_1_content = """
    library(rredlist)
    
    args = commandArgs(trailingOnly = TRUE)
    taxon_group = args[1]
    group_rank = args[2]
    iucn_key = args[3]
    outdir = args[4]
    exclude_extinct = FALSE
    
    # load all IUCN data
    data = c()
    for (i in seq(0, 20, 1)){
      data = c(data,c(rl_sp(key=iucn_key,page = i)))
    }
    
    # get taxon list, class list and status list from data
    taxon_list = c()
    group_list = c()
    status_list = c()
    taxon_id = c()
    for (page in data){
      if (length(page) > 1){
        target_column = which(startsWith(colnames(page),group_rank))
        taxon_list = c(taxon_list,page$scientific_name)
        group_list = c(group_list,page[,target_column])
        status_list = c(status_list,page$category)
        taxon_id = c(taxon_id,page$taxonid)
      }
    }
    
    # exclude extinct taxa if needed
    if (exclude_extinct){
      boolean = !grepl('EX|EW',status_list)
      taxon_list = taxon_list[boolean]
      group_list = group_list[boolean]
      status_list = status_list[boolean]
      taxon_id = taxon_id[boolean]
    }

    # remove all non-species level identifications
    boolean = !grepl('subsp.|ssp.|subpopulation|Subpopulation',taxon_list)
    taxon_list = taxon_list[boolean]
    group_list = group_list[boolean]
    status_list = status_list[boolean]
    taxon_id = taxon_id[boolean]
    
    # select target taxa
    selected_taxon_list = taxon_list[group_list==taxon_group]
    selected_ids = taxon_id[group_list==taxon_group]
    final_sorted_taxon_list = selected_taxon_list
    #final_taxon_list = as.data.frame(cbind(selected_taxon_list,selected_ids))
    #final_sorted_taxon_list = final_taxon_list[order(final_taxon_list$selected_taxon_list),]
    write.table(sort(final_sorted_taxon_list),file=paste0(outdir,'/',taxon_group,"_species_list.txt"), quote=F,row.names=F,sep='	',col.names = FALSE)
    
    
    # get historic data __________________________
    # create new dataframe with species as first column
    historic_assessments = selected_taxon_list
    historic_assessments = as.data.frame(historic_assessments)
    colnames(historic_assessments) = c('species')
    # find historic assessments and fill into dataframe
    log_frequency = 1000
    counter = 1
    for (i in seq(1, length(selected_taxon_list), 1)){
      species = selected_taxon_list[i]
      species_id = selected_ids[i]
      print(paste0('Downloading IUCN history: species ',counter, ' of ',length(selected_taxon_list)))
      #print(species)
      row_id = which(historic_assessments$species == species)
      hist_data <- rl_history(id=species_id,key=iucn_key)
      for (year in hist_data$result$year){
        id = which(hist_data$result$year == year)
        #some species have multiple assignments for some years
        if (length(hist_data$result$code[id])>1){
          historic_assessments[row_id,year] <- hist_data$result$code[id][1]
        }
        else{
          historic_assessments[row_id,year] <- hist_data$result$code[id]
        }
      }
    if (counter %% log_frequency == 0){
      write.table(historic_assessments,file=paste0(outdir,'/',taxon_group,"_iucn_history.txt"), quote=F,row.names=F,sep='	')
    }
      counter = counter+1
    }
    write.table(historic_assessments,file=paste0(outdir,'/',taxon_group,"_iucn_history.txt"), quote=F,row.names=F,sep='	')
    #___________________________________    
    """
    
    script_2_content = """
    library(rredlist)

    args = commandArgs(trailingOnly = TRUE)
    species_list_file = args[1]
    iucn_key = args[2]
    outdir = args[3]
    
    data = read.csv(species_list_file,header = FALSE)
    species_list = data$V1
    status_list = c()
    for (i in 1:length(species_list)){
      species = as.character(species_list[i])
      print(paste0('Extracting current status for ', species,' (',i,' of ',length(species_list),')'))
      iucn_info = rl_search(species,key = iucn_key)
      category = iucn_info$result$category
      if (is.null(category)){
        category = NaN
      }
      status_list = c(status_list,category)
    }
    
    species_cat_df = cbind(as.character(species_list),status_list)
    write.table(species_cat_df,file=paste0(outdir,'/current_status_missing_species.txt'),quote=FALSE,sep = '\t',col.names = FALSE,row.names=FALSE)
    """
    
    rscript_out = os.path.join(output_folder,'rscripts')
    if not os.path.exists(rscript_out):
        os.makedirs(rscript_out)

    if script_id == 'both':
        with open(os.path.join(rscript_out,'get_iucn_status_data_and_species_list.r'),'w') as file:
            file.write(script_1_content)
            file.close()
        with open(os.path.join(rscript_out,'get_current_iucn_status_missing_species.r'),'w') as file:
            file.write(script_2_content)
            file.close()
    elif script_id == 'history':
        with open(os.path.join(rscript_out,'get_iucn_status_data_and_species_list.r'),'w') as file:
            file.write(script_1_content)
            file.close()
    elif script_id == 'current':    
        with open(os.path.join(rscript_out,'get_current_iucn_status_missing_species.r'),'w') as file:
            file.write(script_2_content)
            file.close()
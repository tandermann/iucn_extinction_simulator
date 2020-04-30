#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MCMC-estimation of status transition rates from IUCN record

Created on Mon Oct 28 14:43:44 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import os,sys
import datetime
import iucn_sim.functions as cust_func

# get extinction probs_________________________________________________________
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

def add_arguments(parser):
    parser.add_argument(
        '--species_data',
        required=True,
        metavar='<path>',
        help="File containing species list and current IUCN status of species, as well as generation length (GL) data estimates if available. GL data is only used for '--extinction_probs_mode 0' ('species_data.txt' output from get_iucn_data function).",
    )
    parser.add_argument(
        '--iucn_history',
        required=True,
        metavar='<path>',
        help="File containing IUCN history of the reference group for transition rate estimation ('*_iucn_history.txt' output of get_iucn_data function)."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        metavar='<path>',
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--extinction_probs_mode',
        default=0,
        metavar='N',
        help="Set to '0' to use IUCN defined extinction probabilities (e.g. Mooers et al, 2008 approach), also using available GL data to estimate species-specific extinction probabilities. Set to '1' to simulate extinctions based on recorded extinctions in IUCN history (e.g. Monroe et al, 2019 approach, no GL data is being used)."
    )
    parser.add_argument(
        '--possibly_extinct_list',
        default=0,
        metavar='<path>',
        help="File containing list of taxa that are likely extinct, but that are listed as extant in IUCN, including the year of their assessment as possibly extinct ('possibly_extinct_reference_taxa.txt' output from get_iucn_data function). These species will then be modeled as extinct by the esimate_rates function, which will effect the estimated extinction probabilities when chosing `--extinction_probs_mode 1`",
    )
    parser.add_argument(
        '--rate_samples',
        default=100,
        metavar='N',
        help="How many rates to sample from the posterior transition rate estimates. These rates will be used to populate transition rate q-matrices for downstream simulations. Later on you can still chose to run more simulation replicates than the here specified number of produced transition rate q-matrices, in which case the `run_sim` function will randomely resample from the available q-matrices (default=100, this is ususally sufficient, larger numbers can lead to very high output file size volumes)."
    )
    parser.add_argument(
        '--n_gen',
        default=100000,
        metavar='N',
        help="Number of generations for MCMC for transition rate estimation (default=100000)."
    )
    parser.add_argument(
        '--burnin',
        default=1000,
        metavar='N',
        help="Burn-in for MCMC for transition rate estimation (default=1000)."
    )
    parser.add_argument(
        '--seed',
        default=None,
        help="Set random seed for the MCMC."
    )

def main(args):
    # get user input___________________________________________________________
    input_data = args.species_data
    iucn_history = args.iucn_history
    outdir = args.outdir
    try:
        extinction_probs_mode = int(args.extinction_probs_mode)
    except:
        print('Invalid extinction_probs_mode provided. Please choose between the currenlty available options 0 or 1')
        quit()
    possibly_extinct_list = args.possibly_extinct_list
    n_rep = int(args.rate_samples)
    n_gen = int(args.n_gen)
    burnin = int(args.burnin)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    seed = args.seed
    try:
        random_seed = False
        seed = int(seed)
    except:
        seed = np.random.randint(999999999)
        random_seed = True
    np.random.seed(seed)
    np.savetxt(os.path.join(outdir,'starting_seed.txt'),np.array([seed]),fmt='%i')
        


    # get input data
    species_data_input = pd.read_csv(input_data,sep='\t',header=None).dropna()
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
    iucn_start_year = 2001
    current_year = datetime.datetime.now().year  
    master_stat_time_df = pd.DataFrame(columns=['species']+list(np.arange(iucn_start_year,current_year+1).astype(str)))
    statuses_through_time = pd.read_csv(iucn_history, delimiter = '\t')
    target_columns = [column for column in master_stat_time_df.columns if column in statuses_through_time.columns]
    master_stat_time_df[target_columns] = statuses_through_time[target_columns]
    # treat EW as EX
    master_stat_time_df.replace('EW', 'EX',inplace=True)
    # replace occurrences of NR (not recognized) with nan
    master_stat_time_df.replace('NR', np.nan,inplace=True)
    
    # clean and sort df
    master_stat_time_df = master_stat_time_df.sort_values(by='species')
    master_stat_time_df = master_stat_time_df.drop_duplicates()
    master_stat_time_df.index = np.arange(len(master_stat_time_df))   
    
    # set the assessment at the current year to NE for species without any assessments
    na_row_indeces = np.where(master_stat_time_df.iloc[:,1:].T.isnull().all().values)
    for index in na_row_indeces:
        master_stat_time_df.iloc[index,-1] = 'NE'
    
    # if possibly_extinct_list provided, read that list and set the status for those taxa to extinct, starting at provided year
    if possibly_extinct_list:
        pex_data = pd.read_csv(possibly_extinct_list,sep='\t')
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
    
    # extract most recent valid status for each taxon
    valid_status_dict,most_recent_status_dict,status_series,taxon_series = cust_func.extract_valid_statuses(master_stat_time_df)

    # extinciton prob mode 0: remove all currently extinct taxa
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
    valid_status_dict,most_recent_status_dict,status_series,taxon_series = cust_func.extract_valid_statuses(master_stat_time_df)
    # count current status distribution
    unique, counts = np.unique(status_series, return_counts=True)
    print('Current IUCN status distribution in reference group:',dict(zip(unique, counts)))
    # count how often each status change occurs
    change_type_dict = cust_func.count_status_changes(master_stat_time_df,valid_status_dict)
    print('Summing up years spend in each category...')
    years_in_each_category = cust_func.get_years_spent_in_each_category(master_stat_time_df,valid_status_dict)

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
    if not random_seed:
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
        en_risks_df = pd.DataFrame(np.zeros((len(species_list),dims+1)))
        en_risks_df.columns = ['species']+ ['EN_p_ext_%i'%i for i in np.arange(0,dims)]
        en_risks_df.species = species_list
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
        cr_risks_df = pd.DataFrame(np.zeros((len(species_list),dims+1)))
        cr_risks_df.columns = ['species']+ ['CR_p_ext_%i'%i for i in np.arange(0,dims)]
        cr_risks_df.species = species_list
        cr_risks_df.iloc[:,1:] = cr_risks
        cr_risks_df.to_csv(os.path.join(outdir,'cr_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')
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
        en_risks_selection = en_risks[:,sample_columns]
        cr_risks_selection = cr_risks[:,sample_columns]
        
    elif extinction_probs_mode == 1:
        target_keys = [i for i in sampled_rates_df.status_change.values if i[-2:] == 'EX']
        ex_probs = sampled_rates_df[sampled_rates_df.status_change.isin(target_keys)].iloc[:-1,1:].values.T
        transition_rates = sampled_rates_df[~sampled_rates_df.status_change.isin(target_keys)].iloc[:30,1:]
    
    for i in np.arange(n_rep):
        rates_i = transition_rates.iloc[:,i]
        sys.stdout.write('\rProgress: %i %%'%int(((i+1)/n_rep)*100))

        # for each rep (i), create list of q-matrices, 1 for each species
        if extinction_probs_mode == 0:
            en_risks_rep = en_risks_selection[:,i]
            cr_risks_rep = cr_risks_selection[:,i]  
            q_matrix_list_i = []
            for j,__ in enumerate(species_list):
                en_risk = en_risks_rep[j]
                cr_risk = cr_risks_rep[j]
                status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category    
                q_matrix = cust_func.qmatrix(rates_i, status_specific_p_e)
                q_matrix_list_i.append([q_matrix])
        elif extinction_probs_mode == 1:
            q_matrix_list_i = []
            status_specific_p_e = ex_probs[i]
            q_matrix = cust_func.qmatrix(rates_i, status_specific_p_e)
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
    cust_func.save_obj(final_output_data,os.path.join(outdir,'simulation_input_data.pkl'))
    #__________________________________________________________________________










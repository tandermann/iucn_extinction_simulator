#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 14:43:44 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
#import matplotlib.pyplot as plt
np.random.seed(1234)
import subprocess
import os, glob
import datetime
import sys
import iucn_sim.functions as cust_func


# get extinction probs_________________________________________________________
def p_e_year(years,p_e):
    pe_year = 1-(1-float(p_e))**(1/years)
    return pe_year

def sample_rate(count, tot_time, n_samples = 1, range_factor = 100, n_bins = 10000):
    def get_loglik(count, dT, rate):
        return np.log(rate)*count - dT*rate
    mle = count/tot_time
    if count == 0:
        minRange = 1*np.exp(-10)
        maxRange = range_factor/tot_time
    else:
        minRange = mle/range_factor
        maxRange = mle*range_factor
    rates = np.linspace(minRange, maxRange,n_bins)
    lik = get_loglik(count,tot_time,rates)
    lik = lik - np.max(lik)
    sample_rates = np.random.choice(rates,size=n_samples,p=np.exp(lik)/np.sum(np.exp(lik)),replace=1)
    #plt.plot(np.log(rates),lik)
    #np.log(mle)
    return sample_rates

def add_arguments(parser):
    parser.add_argument(
        '--input_data',
        required=True,
        help="Path to generation length (GL) data: first column taxon list, followed by n columns of GL values."
    )
    parser.add_argument(
        '--reference_group',
        default=0,
        help="Name of taxonomic group (or list of groups) to be used for calculating status transition rates (e.g. 'Mammalia' or 'Rodentia,Chiroptera'). Alternatively provide path to text file containing a list of species names, compatible with IUCN taxonomy (>1000 species recommended). If none provided, the input species list with GL data will be used for calculating transition rates. Tip: Use precompiled group for significantly faster processing (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/)"
    )
    parser.add_argument(
        '--reference_rank',
        default=0,
        help="Provide the taxonomic rank of the provided reference group(s). E.g. in case of 'Mammalia', provide 'class' for this flag, in case of 'Rodentia,Chiroptera' provide 'order,order'. Has to be at least 'Family' or above. This flag is not needed if species list is provided as reference_group or if reference group is already pre-compiled."
    )
    parser.add_argument(
        '--n_rep',
        default=0,
        help="How many independent simulation replicates to run (default == number of provided GL value columns per species under 'input_data' flag)."
    )
    parser.add_argument(
        '--iucn_key',
        default=0,
        help="Provide your IUCN API key (see https://apiv3.iucnredlist.org/api/v3/token) for downloading IUCN history of your provided reference group. Not required if using precompiled reference group."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--github_repo',
        default=0,
        help="Provide path to a copy of the iucn_extinction_simulator GitHub repo. This is needed to search for pre-compiled files that significantly shorten the computation time. Download link: https://github.com/tobiashofmann88/iucn_extinction_simulator/archive/master.zip (make sure to unzip the downloaded repo)."
    )
    parser.add_argument(
        '--status_list',
        default=0,
        help="Provide a text file containing a valid IUCN status (LC,NT,VU,EN,CR,DD) for each species, separated by newline (same order as species names provided under --input_data)."
    )
    parser.add_argument(
        '--allow_precompiled_iucn_data',
        default=1,
        help="Set this flag to 0 if you want to avoid using precompiled IUCN history data. By default (1) this data is used if available for your specified reference organism group."
    )
    parser.add_argument(
        '--rate_sampling_range',
        default=100,
        help="Set the range of status change rate sampling. Default = 100 will sample from a uniform distribution centered at the Maximum Likelihood Estimate (MLE) with a range of MLE/100 and MLE*100."
    )
    parser.add_argument(
        '--rate_bins',
        default=10000,
        help="Set number of bins from which to sample the status transition rates (default=10000)."
    )
    
    
def main(args):
        
    # get user input
    input_data = args.input_data
    taxon_reference_group = args.reference_group
    reference_rank = args.reference_rank
    n_rep = int(args.n_rep)
    iucn_key = args.iucn_key
    outdir = args.outdir
    github_repo = args.github_repo
    allow_precompiled_iucn_data = args.allow_precompiled_iucn_data
    rate_sampling_range = args.rate_sampling_range
    rate_bins = args.rate_bins
    status_list = args.status_list
        
    taxon_reference_groups = taxon_reference_group.split(',')
    reference_ranks = reference_rank.split(',')
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # get gl data to calculate en and cr extinction risk for all species
    gl_data = pd.read_csv(input_data,sep='\t',header=None)
    species_list = gl_data.iloc[:,0].values
    gl_data_available = False
    if gl_data.shape[1] > 1:
        gl_matrix = gl_data.iloc[:,1:].values
        gl_data_available = True
    
    status_list_data = list(pd.read_csv(status_list,sep='\t',header=None).iloc[:,0].values)

    if len(species_list) != len(status_list_data):
        print('Error: Length of provided status list does not match length of species list provided as --input_data!')
        quit()

    # get IUCN history_________________________________________________________
    iucn_outdir = os.path.join(outdir,'iucn_data')    
    if not os.path.exists(iucn_outdir):
        os.makedirs(iucn_outdir)
    precompiled_taxon_groups=[]
    if github_repo:
        precompiled_taxon_group_files = glob.glob('%s/data/precompiled/iucn_history/*_iucn_history.txt'%github_repo)
        precompiled_taxon_groups = [str.lower(os.path.basename(filename).replace('_iucn_history.txt','')) for filename in precompiled_taxon_group_files]
    iucn_history_files = []
    for i,taxon_group in enumerate(taxon_reference_groups):
        if str.lower(taxon_group) in precompiled_taxon_groups and allow_precompiled_iucn_data:
            print('Loading precompiled IUCN history data from %s/data/precompiled/iucn_history/'%github_repo)
            iucn_history_files.append([file for file in precompiled_taxon_group_files if os.path.basename(file).startswith(str.upper(taxon_group))][0])
        else:
            print('Fetching IUCN history using rredlist')
            rank = reference_ranks[i]
            iucn_cmd = ['Rscript','./iucn_sim/r_scripts/get_iucn_status_data_and_species_list.r', str.upper(taxon_group), str.lower(rank), iucn_key, iucn_outdir]
            #iucn_error_file = os.path.join(iucn_outdir,'get_iucn_status_data_and_species_list_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            seqtk = subprocess.Popen(iucn_cmd)
            seqtk.wait()
            iucn_history_files.append(os.path.join(iucn_outdir,'%s_iucn_history.txt'%str.upper(taxon_group)))

    
    # process the IUCN history data____________________________________________
    iucn_start_year = 2000
    current_year = datetime.datetime.now().year  
    master_stat_time_df = pd.DataFrame(columns=['species']+list(np.arange(iucn_start_year,current_year+1).astype(str)))
    for iucn_history in iucn_history_files:
        statuses_through_time = pd.read_csv(iucn_history, delimiter = '\t')
        data = np.zeros([len(statuses_through_time),(current_year-iucn_start_year)+2]).astype(object)
        data[:] = np.nan
        data[:,0] = statuses_through_time.species.values
        for year in statuses_through_time.columns[1:]:
            array_col = (int(year)-iucn_start_year)+1
            if array_col >= 0:
                data[:,array_col] = statuses_through_time[year].values
        master_stat_time_df = master_stat_time_df.append(pd.DataFrame(data,columns = ['species']+list(np.arange(iucn_start_year,current_year+1).astype(str))),ignore_index=True)
    # check if we have sufficient number of species for rate estimation
    if len(master_stat_time_df) < 1000:
        print('\n\n','#'*50,'\nWarning: Only %i species in reference dataset. This may not be sufficient for proper estimation of status transition rates. It is recommended to choose a larger reference group encompassing >1000 species!'%len(master_stat_time_df),'#'*50,'\n\n')
    # treat EW as EX
    master_stat_time_df.replace('EW', 'EX',inplace=True)
    # replace EX with NaN
    master_stat_time_df.replace('EX', np.nan,inplace=True)
    # write IUCN history df to file
    master_stat_time_df.to_csv(os.path.join(iucn_outdir,'iucn_history_reference_taxa.txt'),sep='\t')
    # extract most recent valid status for each taxon
    valid_status_dict,most_recent_status_dict,status_series,taxon_series = cust_func.extract_valid_statuses(master_stat_time_df)
    # count current status distribution
    unique, counts = np.unique(status_series, return_counts=True)
    print('Current IUCN status distribution in specified reference group:',dict(zip(unique, counts)))
    # count how often each status change occurs
    change_type_dict = cust_func.count_status_changes(master_stat_time_df,valid_status_dict)
    print('Observed the following status changes (type: count)')
    print(change_type_dict)
    print('Filling in missing status data...')
    df_for_year_count = cust_func.fill_dataframe(master_stat_time_df,valid_status_dict,iucn_start_year,current_year)
    # count years spent in each category for all species
    years_in_each_category = cust_func.get_years_spent_in_each_category(df_for_year_count)
    # write the status change data to file
    final_years_count_array = np.array([list(years_in_each_category.keys()),list(years_in_each_category.values())]).T
    np.savetxt(os.path.join(outdir,'years_spent_in_each_category.txt'),final_years_count_array,fmt='%s\t%s')
    change_type_dict_array = np.array([list(change_type_dict.keys()),list(change_type_dict.values())]).T
    np.savetxt(os.path.join(outdir,'change_type_dict.txt'),change_type_dict_array,fmt='%s\t%s')   

    # sample transition rates for all types of changes_________________________
    if n_rep > 0:
        sim_reps = n_rep
    else:
        if gl_data_available:
            sim_reps = gl_matrix.shape[1]
        else:
            print('ERROR: Use --n_rep flag to set the number of simulation replicates.')
            quit()
    print('Preparing data for %i simulation replicates'%sim_reps)
    status_change_coutn_df = pd.DataFrame(data=np.zeros([6,6]).astype(int),index = ['LC','NT','VU','EN','CR','DD'],columns=['LC','NT','VU','EN','CR','DD'])
    for status_change in change_type_dict.keys():
        states = status_change.split('->')
        original_state = states[0]
        new_state = states[1]
        count = change_type_dict[status_change]
        status_change_coutn_df.loc[original_state,new_state] = count
    sampled_rates_df = pd.DataFrame(columns = ['status_change']+ ['rate_%i'%i for i in np.arange(0,sim_reps)])
    for status_a in status_change_coutn_df.columns:
        column = status_change_coutn_df[status_a]
        for status_b in column.index:
            if not status_a == status_b:
                count = column[status_b]
                total_time = years_in_each_category[status_a]
                rates = sample_rate(count, total_time, n_samples = sim_reps, range_factor = rate_sampling_range, n_bins = rate_bins)
                sampled_rates_df = sampled_rates_df.append(pd.DataFrame(data=np.matrix(['%s->%s'%(status_a,status_b)]+list(rates)),columns = ['status_change']+ ['rate_%i'%i for i in np.arange(0,sim_reps)]),ignore_index=True)
    sampled_rates_df[['rate_%i'%i for i in np.arange(0,sim_reps)]] = sampled_rates_df[['rate_%i'%i for i in np.arange(0,sim_reps)]].apply(pd.to_numeric)
    sampled_rates_df.to_csv(os.path.join(outdir,'sampled_status_change_rates.txt'),sep='\t',index=False,float_format='%.8f')

    # get current status for all species we want to simulate___________________
    # get list of all species that we don't have IUCN data for already
    if status_list:
        current_status_list = status_list_data
    else:
        missing_species = np.array([species for species in species_list if not species in most_recent_status_dict.keys()])
        if len(missing_species) > 0:
            missing_species_file = os.path.join(outdir,'missing_species_list.txt')
            np.savetxt(missing_species_file,missing_species,fmt='%s')
            # extract the current status for those missing species
            print('Extracting current status for missing species...')
            iucn_cmd = ['Rscript','./iucn_sim/r_scripts/get_current_iucn_status_missing_species.r', missing_species_file, iucn_key, iucn_outdir]
            #iucn_error_file = os.path.join(iucn_outdir,'get_current_iucn_status_missing_species_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            seqtk = subprocess.Popen(iucn_cmd)
            seqtk.wait()
            missing_species_status_file = os.path.join(iucn_outdir,'current_status_missing_species.txt')
            if os.path.exists(missing_species_status_file):
                missing_species_status_data = pd.read_csv(missing_species_status_file,sep='\t',header=None)
            # clean up temporary files
            if os.path.exists(missing_species_status_file):
                os.remove(missing_species_status_file)
            if os.path.exists(missing_species_file):
                os.remove(missing_species_file)
        current_status_list = []
        for species in species_list:
            if species in most_recent_status_dict.keys():
                current_status = most_recent_status_dict[species]
            else:
                current_status = missing_species_status_data[missing_species_status_data[0]==species][1].values[0]
                if current_status not in ['LC','NT','VU','EN','CR','DD']:
                    print('Species %s was not found in IUCN database. Current status is therefore set to DD (data deficient).'%species)
                    current_status = 'DD'
            current_status_list.append(current_status)
    final_df_current_status = pd.DataFrame(np.array([species_list,current_status_list]).T,columns=['species','current_status'])
    final_df_current_status.to_csv(os.path.join(iucn_outdir,'current_status_all_species.txt'),sep='\t',index=False)


    # calculate EN and CR extinction risk using GL information_________________
    # calculate yearly extinction risks for categories EN and CR
    if gl_data_available:
        en_risks = np.array([p_e_year(np.minimum(np.maximum([20]*len(gl_array),5*gl_array),100),0.2) for gl_array in gl_matrix])
    else:
        print('Warning: No generation length (GL) data found. Extinction risks for status EN and CR are calculated without using GL data.')
        en_risks = np.array([[p_e_year(20,0.2)]*sim_reps]*len(species_list))
    if en_risks.shape[1] == 1:
        en_risks = np.array([en_risks for i in range(sim_reps)])[:,:,0].T        
    en_risks_df = pd.DataFrame(np.zeros((len(species_list),sim_reps+1)))
    en_risks_df.columns = ['species']+ ['en_extinction_risk_yearly_%i'%i for i in np.arange(0,sim_reps)]
    en_risks_df.species = species_list
    en_risks_df.iloc[:,1:] = en_risks
    en_risks_df.to_csv(os.path.join(outdir,'en_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')

    if gl_data_available:
        cr_risks = np.array([p_e_year(np.minimum(np.maximum([10]*len(gl_array),3*gl_array),100),0.5) for gl_array in gl_matrix])
    else:
        cr_risks = np.array([[p_e_year(10,0.5)]*sim_reps]*len(species_list))
    if cr_risks.shape[1] == 1:
        cr_risks = np.array([cr_risks for i in range(sim_reps)])[:,:,0].T
    cr_risks_df = pd.DataFrame(np.zeros((len(species_list),sim_reps+1)))
    cr_risks_df.columns = ['species']+ ['cr_extinction_risk_yearly_%i'%i for i in np.arange(0,sim_reps)]
    cr_risks_df.species = species_list
    cr_risks_df.iloc[:,1:] = cr_risks
    cr_risks_df.to_csv(os.path.join(outdir,'cr_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')






# code scraps__________________________________________________________________    
#    # count years spent in each category for all species, use df without DD values for this
#    # replace DD with NaN
#    master_stat_time_df_no_dd = master_stat_time_df.copy()
#    master_stat_time_df_no_dd.replace('DD', np.nan,inplace=True)
#    # extract the valid status series for the no_dd df
#    valid_status_dict_no_dd,most_recent_status_dict_no_dd,status_series_no_dd,taxon_series_no_dd = cust_func.extract_valid_statuses(master_stat_time_df_no_dd)
#    # count how often each status change occurs
#    change_type_dict_no_dd = cust_func.count_status_changes(master_stat_time_df_no_dd,valid_status_dict_no_dd)
#    # fill in missing fields in dataframe
#    print('Filling in missing status data...')
#    df_for_year_count = cust_func.fill_dataframe(master_stat_time_df_no_dd,valid_status_dict_no_dd,iucn_start_year,current_year)
#    # count years spent in each category
#    years_in_each_category = cust_func.get_years_spent_in_each_category(df_for_year_count)
#    change_type_dict_no_dd_array = np.array([list(change_type_dict_no_dd.keys()),list(change_type_dict_no_dd.values())]).T
#    np.savetxt(os.path.join(outdir,'change_type_dict_no_dd.txt'),change_type_dict_no_dd_array,fmt='%s\t%s')   


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download IUCN data and estimate status transition rates

Created on Mon Oct 28 14:43:44 2019
@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
#import matplotlib.pyplot as plt
np.random.seed(1234)
import subprocess
import os
import datetime
import iucn_sim.functions as cust_func
from urllib.request import urlopen
from io import StringIO

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
        help="How many different transition-rate and extinction risk estimates to produce for simulations (default == number of provided GL value columns per species under 'input_data' flag)."
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
        '--n_gen',
        default=100000,
        help="Number of generations for MCMC for transition rate estimation (default=100000)."
    )
    parser.add_argument(
        '--burnin',
        default=1000,
        help="Burn-in for MCMC for transition rate estimation (default=1000)."
    )
    
    
def main(args):   
    
#    import argparse
#    p = argparse.ArgumentParser()
#    args = p.parse_args()    
#    args.input_data = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/gl_data_carnivora_some_missing_gl.txt'
#    args.reference_group = 'Mammalia'
#    args.reference_rank = 'class'
#    args.n_rep = 0
#    args.iucn_key = '01524b67f4972521acd1ded2d8b3858e7fedc7da5fd75b8bb2c5456ea18b01ba'
#    args.outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/carnivora_output_test/transition_rates'
#    args.allow_precompiled_iucn_data = 1
#    args.n_gen = 100000
#    args.burnin = 1000
#    args.status_list = 0    
    
    # get user input
    input_data = args.input_data
    taxon_reference_group = args.reference_group
    reference_rank = args.reference_rank
    n_rep = int(args.n_rep)
    iucn_key = args.iucn_key
    outdir = args.outdir
    allow_precompiled_iucn_data = int(args.allow_precompiled_iucn_data)
    n_gen = int(args.n_gen)
    burnin = int(args.burnin)
    status_list = args.status_list
    
    # create the r-scripts to be used later on:
    cust_func.write_r_scripts(outdir)
    
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
    
    if status_list:
        status_list_data = list(pd.read_csv(status_list,sep='\t',header=None).iloc[:,0].values)

        if len(species_list) != len(status_list_data):
            print('Error: Length of provided status list does not match length of species list provided as --input_data!')
            quit()

    # get IUCN history_________________________________________________________
    iucn_outdir = os.path.join(outdir,'iucn_data')    
    if not os.path.exists(iucn_outdir):
        os.makedirs(iucn_outdir)
    precompiled_taxon_groups = []
    precompiled_taxon_group_files = []
    if allow_precompiled_iucn_data:
        for taxon_group in taxon_reference_groups:
            try:
                # look for precompiled files online    
                url = 'https://raw.githubusercontent.com/tobiashofmann88/iucn_extinction_simulator/master/data/precompiled/iucn_history/%s_iucn_history.txt'%taxon_group.upper()
                urlpath =urlopen(url)
                string = urlpath.read().decode('utf-8')        
                string_input = StringIO(string)
                ref_group_data = pd.read_csv(string_input, sep="\t")
                hist_outfile = os.path.join(iucn_outdir,'iucn_history_%s.txt'%str.lower(taxon_group))
                ref_group_data.to_csv(hist_outfile,sep='\t',index=False)
                precompiled_taxon_groups.append(str.lower(taxon_group))
                precompiled_taxon_group_files.append(hist_outfile)
            except:
                pass

    iucn_history_files = []
    for i,taxon_group in enumerate(taxon_reference_groups):
        if str.lower(taxon_group) in precompiled_taxon_groups and allow_precompiled_iucn_data:
            print('Using precompiled IUCN history data for %s.'%taxon_group)                    
            iucn_history_files.append([file for file in precompiled_taxon_group_files if os.path.basename(file).endswith(str.lower(taxon_group)+'.txt')][0])
        else:
            print('Fetching IUCN history using rredlist')
            rank = reference_ranks[i]
            iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_iucn_status_data_and_species_list.r'), str.upper(taxon_group), str.lower(rank), iucn_key, iucn_outdir]
            if not iucn_key:
                quit('***IUCN-KEY ERROR:*** Need to download IUCN history for specified reference group. Please provide a valid IUCN key (using the --iucn_key flag). Alternatively choose a precompiled reference group (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/).'%(taxon_group))
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
    status_change_coutn_df = pd.DataFrame(data=np.zeros([6,6]).astype(int),index = ['LC','NT','VU','EN','CR','DD'],columns=['LC','NT','VU','EN','CR','DD'])
    for status_change in change_type_dict.keys():
        states = status_change.split('->')
        original_state = states[0]
        new_state = states[1]
        count = change_type_dict[status_change]
        status_change_coutn_df.loc[original_state,new_state] = count
    status_change_coutn_df.to_csv(os.path.join(outdir,'status_change_counts.txt'),sep='\t',index=True)
    print('Counted the following transition occurrences in IUCN history of reference group:')
    print(status_change_coutn_df)
    print('Estimating rates for %i simulation replicates'%sim_reps)
    sampled_rates_df = pd.DataFrame(columns = ['status_change']+ ['rate_%i'%i for i in np.arange(0,sim_reps)])
    for status_a in status_change_coutn_df.columns:
        row = status_change_coutn_df.loc[status_a]
        for status_b in row.index.values:
            if not status_a == status_b:
                count = row[status_b]
                total_time = years_in_each_category[status_a]
                rates = sample_rate_mcmc(count, total_time, n_samples = sim_reps, n_gen = n_gen, burnin = burnin)
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
            print('Extracting current status for species that were not found in reference group...')
            iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_current_iucn_status_missing_species.r'), missing_species_file, iucn_key, iucn_outdir]
            #iucn_error_file = os.path.join(iucn_outdir,'get_current_iucn_status_missing_species_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            if not iucn_key:
                quit('***IUCN-KEY ERROR:*** Trying to download current status for species %s from IUCN. Please provide a valid IUCN key (using the --iucn_key flag) to access IUCN data. Alternatively provide the current IUCN status for all your input species using the --status_list flag.'%(str(list(missing_species))))
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
                    print('Species %s was not found in IUCN database. Current status is therefore set to NE (Not Evaluated).'%species)
                    current_status = 'NE'
            current_status_list.append(current_status)
    final_df_current_status = pd.DataFrame(np.array([species_list,current_status_list]).T,columns=['species','current_status'])
    final_df_current_status.to_csv(os.path.join(iucn_outdir,'current_status_all_species.txt'),sep='\t',index=False)


    # calculate EN and CR extinction risk using GL information_________________
    # calculate yearly extinction risks for categories EN and CR
    if gl_data_available:
        en_risks = []
        for gl_array in gl_matrix:
            #replace all nan values with the standard en extinction risk
            en_risks_species = p_e_year(np.minimum(np.maximum([20]*len(gl_array),5*gl_array),100),0.2)
            n_nan = len(en_risks_species[en_risks_species!=en_risks_species])
            en_risks_species[en_risks_species!=en_risks_species] = [p_e_year(20,0.2)]*n_nan
            en_risks.append(en_risks_species)
        en_risks = np.array(en_risks)
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
        cr_risks = []
        for gl_array in gl_matrix:
            #replace all nan values with the standard en extinction risk
            cr_risks_species = p_e_year(np.minimum(np.maximum([10]*len(gl_array),3*gl_array),100),0.5)
            n_nan = len(cr_risks_species[cr_risks_species!=cr_risks_species])
            cr_risks_species[cr_risks_species!=cr_risks_species] = [p_e_year(10,0.5)]*n_nan
            cr_risks.append(cr_risks_species)
        cr_risks = np.array(cr_risks)
    else:
        cr_risks = np.array([[p_e_year(10,0.5)]*sim_reps]*len(species_list))
    if cr_risks.shape[1] == 1:
        cr_risks = np.array([cr_risks for i in range(sim_reps)])[:,:,0].T
    cr_risks_df = pd.DataFrame(np.zeros((len(species_list),sim_reps+1)))
    cr_risks_df.columns = ['species']+ ['cr_extinction_risk_yearly_%i'%i for i in np.arange(0,sim_reps)]
    cr_risks_df.species = species_list
    cr_risks_df.iloc[:,1:] = cr_risks
    cr_risks_df.to_csv(os.path.join(outdir,'cr_extinction_risks_all_species.txt'),sep='\t',index=False, float_format='%.12f')




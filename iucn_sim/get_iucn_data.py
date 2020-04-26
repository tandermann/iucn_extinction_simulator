#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download IUCN data for future simulations

Created on Thu Apr  9 12:05:39 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import subprocess
import os
import iucn_sim.functions as cust_func
from urllib.request import urlopen
from io import StringIO
import datetime



def add_arguments(parser):
    parser.add_argument(
        '--reference_group',
        default=0,
        metavar='taxon-group',
        help="Name of taxonomic group (or list of groups) to be used for calculating status transition rates (e.g. 'Mammalia' or 'Rodentia,Chiroptera'). Alternatively provide path to text file containing a list of species names, compatible with IUCN taxonomy (>1000 species recommended). If none provided, the target species list ('--target_species_list') will be used for calculating transition rates. Tip: Use precompiled group for faster processing and in case you don't have an IUCN key (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/ or request specific groups to be added: tobias.andermann@bioenv.gu.se)"
    )
    parser.add_argument(
        '--reference_rank',
        default=0,
        metavar='rank',
        help="Provide the taxonomic rank of the provided reference group(s). E.g. in case of 'Mammalia', provide 'class' for this flag, in case of 'Rodentia,Chiroptera' provide 'order,order'. Has to be at least 'Family' or above. This flag is not needed if species list is provided as reference_group or if reference group is already pre-compiled."
    )
    parser.add_argument(
        '--target_species_list',
        default=1,
        metavar='<path>',
        help="File containing the list of species that you want to simulate future extinctions for. In case you have generation length (GL) data available, provide the file containing the GL data for each species here (including the species names). This function will output one central data file for downstream processing that contains the current status information as well as the GL data (if available) for each species. You can provide multiple GL values for each species, e.g. several randomely sampled values from the GL uncertainty interval of a given species. Set this flag to 0 if you want to supress downloading of current status information, e.g. if you already have current status information for your species (may be necessary if you don't have a valid IUCN key). Set to 1 if you want to use the same taxa that are present in the reference group. See https://github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/ for examples of the format of GL data input files and the format of the output file conataining current status information.",
    )
    parser.add_argument(
        '--outdir',
        required=True,
        metavar='<path>',
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--iucn_key',
        default=0,
        metavar='<IUCN-key>',
        help="Provide your IUCN API key (see https://apiv3.iucnredlist.org/api/v3/token) for downloading IUCN history of your provided reference group. Not required if using precompiled reference group and a manually compiled current status list (to be used in the 'transition_rates' function). Also not required if all species in your target_species_list are present in the precompiled reference_group)."
    )
    parser.add_argument(
        '--no_online_sync',
        action='store_true',
        help='Turn off the online-search for precompiled IUCN history files for your reference group.',
        default=False
    )


def main(args):
    # get user input
    taxon_reference_group = args.reference_group
    reference_rank = args.reference_rank
    species_list_file = args.target_species_list
    outdir = args.outdir
    iucn_key = args.iucn_key
    no_precompiled_iucn_data = args.no_online_sync
    
    # create the r-scripts to be used later on:
    cust_func.write_r_scripts(outdir)
    
    if taxon_reference_group:
        taxon_reference_groups = taxon_reference_group.split(',')
    else:
        print('No reference group provided. Use the --reference_group to provide the name of a taxonomic group, e.g. Mammalia')
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
    if not no_precompiled_iucn_data:
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
        if str.lower(taxon_group) in precompiled_taxon_groups and not no_precompiled_iucn_data:
            print('Using precompiled IUCN history data for %s.'%taxon_group)                    
            iucn_history_files.append([file for file in precompiled_taxon_group_files if os.path.basename(file).startswith(str.upper(taxon_group))][0])
        else:
            print('Fetching IUCN history using rredlist')
            rank = reference_ranks[i]
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


    # get current IUCN status all target species_______________________________
    if species_list_file:
        # process the IUCN history data
        iucn_start_year = 2001
        current_year = datetime.datetime.now().year  
        master_stat_time_df = pd.DataFrame(columns=['species']+list(np.arange(iucn_start_year,current_year+1).astype(str)))
        statuses_through_time = pd.read_csv(iucn_history_file, delimiter = '\t')
        target_columns = [column for column in master_stat_time_df.columns if column in statuses_through_time.columns]
        master_stat_time_df[target_columns] = statuses_through_time[target_columns]
    
        # check if we have sufficient number of species for rate estimation
        if len(master_stat_time_df) < 1000:
            print('\n\n%s'%('#'*50),'\nWarning: Only %i species in reference dataset. This may not be sufficient for proper estimation of status transition rates. It is recommended to choose a larger reference group encompassing >1000 species at the least!\n\nContinuing processing IUCN history data of reference group ...\n%s\n\n'%(len(master_stat_time_df),'#'*50))
        # treat EW as EX
        master_stat_time_df.replace('EW', 'EX',inplace=True)
        
        # clean and sort df
        master_stat_time_df = master_stat_time_df.sort_values(by='species')
        master_stat_time_df = master_stat_time_df.drop_duplicates()
        master_stat_time_df.index = np.arange(len(master_stat_time_df))   
        
        # set most recent status to NE for species without any IUCN status information
        na_row_indeces = np.where(master_stat_time_df.iloc[:,1:].T.isnull().all().values)
        for index in na_row_indeces:
            master_stat_time_df.iloc[index,-1] = 'NE'
        
        # extract most recent valid status for each taxon from reference group
        valid_status_dict_refgroup,most_recent_status_dict_refgroup,status_series_refgroup,taxon_series_refgroup = cust_func.extract_valid_statuses(master_stat_time_df)

        # see if we want to only stick with the reference taxa
        try:
            species_list_file = int(species_list_file)
            gl_data = []
            joined_df = pd.DataFrame([taxon_series_refgroup,status_series_refgroup]).T
        except:
            # get list of species we want to simulate
            species_list_data = pd.read_csv(species_list_file,sep='\t',header=None)
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
            except:
                gl_data = []
            
            # see which species are not in reference group and hence need current IUCN status extraction
            remaining_species_wo_iucn_status = sorted(list(set(species_list)-set(taxon_series_refgroup)))
            
            fraction_overlap = np.round(1-(len(remaining_species_wo_iucn_status)/len(species_list)),3)        
            if fraction_overlap <= 0.5:
                print('\n\n%s'%('#'*50),'\nWarning: Only %.3f of target species found in reference group. Is this intended? You may want to reconsider chosing a reference group encompassing more of your target species.\n\nContinuing with downloading current status information for all remaining species...\n%s\n\n'%(fraction_overlap,'#'*50))
            else:
                print('\n\nA fraction of %.3f of the specified target species is present in reference group.\n\n'%fraction_overlap)
                
            tmp_outdir = os.path.join(outdir,'other_files')
            if not os.path.exists(tmp_outdir):
                os.makedirs(tmp_outdir)
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
            target_reference_taxa = list(set(species_list)-set(current_status_missing_taxa[0].values))
            status_series_refgroup = np.array(status_series_refgroup).astype(str)
            taxon_series_refgroup = np.array(taxon_series_refgroup).astype(str)
            status_reference_taxa = [status_series_refgroup[taxon_series_refgroup == i][0] for i in target_reference_taxa]
            current_status_reference_taxa = pd.DataFrame(data = np.array([target_reference_taxa,status_reference_taxa]).T)
            joined_df = pd.concat([current_status_missing_taxa,current_status_reference_taxa],ignore_index=True).sort_values(by=[0])

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




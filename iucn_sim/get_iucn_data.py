#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 12:05:39 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
#import matplotlib.pyplot as plt
#np.random.seed(1234)
import subprocess
import os
import datetime
import iucn_sim.functions as cust_func
from urllib.request import urlopen
from io import StringIO



def add_arguments(parser):
    parser.add_argument(
        '--reference_group',
        default=0,
        metavar='Mammalia',
        help="Name of taxonomic group (or list of groups) to be used for calculating status transition rates (e.g. 'Mammalia' or 'Rodentia,Chiroptera'). Alternatively provide path to text file containing a list of species names, compatible with IUCN taxonomy (>1000 species recommended). If none provided, the input species list with GL data will be used for calculating transition rates. Tip: Use precompiled group for significantly faster processing (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/)"
    )
    parser.add_argument(
        '--reference_rank',
        default=0,
        metavar='class',
        help="Provide the taxonomic rank of the provided reference group(s). E.g. in case of 'Mammalia', provide 'class' for this flag, in case of 'Rodentia,Chiroptera' provide 'order,order'. Has to be at least 'Family' or above. This flag is not needed if species list is provided as reference_group or if reference group is already pre-compiled."
    )
    parser.add_argument(
        '--iucn_key',
        default=0,
        metavar='<IUCN-key>',
        help="Provide your IUCN API key (see https://apiv3.iucnredlist.org/api/v3/token) for downloading IUCN history of your provided reference group. Not required if using precompiled reference group."
    )
    parser.add_argument(
        '--outdir',
        required=True,
        metavar='<path>',
        help="Provide path to outdir where results will be saved."
    )
    parser.add_argument(
        '--allow_precompiled_iucn_data',
        default=1,
        metavar='0/1',
        help="Set this flag to 0 if you want to avoid using precompiled IUCN history data. By default (1) this data is used if available for your specified reference organism group."
    )


    
    import argparse
    p = argparse.ArgumentParser()
    args = p.parse_args()    
    args.reference_group = 'Aves'
    args.reference_rank = 'class'
    args.iucn_key = '01524b67f4972521acd1ded2d8b3858e7fedc7da5fd75b8bb2c5456ea18b01ba'
    args.outdir = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/birds_output_test/iucn_data'
    args.allow_precompiled_iucn_data = 1
     


def main(args):
    # get user input
    input_data = args.input_data
    taxon_reference_group = args.reference_group
    reference_rank = args.reference_rank
    iucn_key = args.iucn_key
    outdir = args.outdir
    allow_precompiled_iucn_data = int(args.allow_precompiled_iucn_data)
    
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

    # get gl data to calculate en and cr extinction risk for all species
    gl_data = pd.read_csv(input_data,sep='\t',header=None)
    species_list = gl_data.iloc[:,0].values
    # replace underscores in species name in case they are present
    species_list = np.array([i.replace('_',' ') for i in species_list])
    # Check if all species names are binomial
    for species in species_list:
        if len(species.split(' ')) != 2:
            print('ERROR','*'*50,'\nABORTED: All provided species names provided under --input_data flag must be binomial! Found non binomial name:\n%s\n'%species,'*'*55)
            quit()
    gl_data_available = False
    if gl_data.shape[1] > 1:
        gl_matrix = gl_data.iloc[:,1:].values
        gl_data_available = True


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
                hist_outfile = os.path.join(iucn_outdir,os.path.basename(url))
                ref_group_data.to_csv(hist_outfile,sep='\t',index=False)
                precompiled_taxon_groups.append(str.lower(taxon_group))
                precompiled_taxon_group_files.append(hist_outfile)
            except:
                pass

    iucn_history_files = []
    for i,taxon_group in enumerate(taxon_reference_groups):
        if str.lower(taxon_group) in precompiled_taxon_groups and allow_precompiled_iucn_data:
            print('Using precompiled IUCN history data for %s.'%taxon_group)                    
            iucn_history_files.append([file for file in precompiled_taxon_group_files if os.path.basename(file).startswith(str.upper(taxon_group))][0])
        else:
            print('Fetching IUCN history using rredlist')
            rank = reference_ranks[i]
            iucn_cmd = ['Rscript',os.path.join(outdir,'rscripts/get_iucn_status_data_and_species_list.r'), str.upper(taxon_group), str.lower(rank), iucn_key, iucn_outdir]
            if not iucn_key:
                print('ERROR','*'*50,'\nIUCN-KEY ERROR: Need to download IUCN history for specified reference group: %s. Please provide a valid IUCN key (using the --iucn_key flag). Alternatively choose a precompiled reference group (see available groups at github.com/tobiashofmann88/iucn_extinction_simulator/data/precompiled/iucn_history/)'%(taxon_group),'*'*55)
                quit()
            #iucn_error_file = os.path.join(iucn_outdir,'get_iucn_status_data_and_species_list_error_file.txt')
            #with open(iucn_error_file, 'w') as err:
            run_iucn_cmd = subprocess.Popen(iucn_cmd)
            run_iucn_cmd.wait()
            iucn_history_files.append(os.path.join(iucn_outdir,'%s_iucn_history.txt'%str.upper(taxon_group)))





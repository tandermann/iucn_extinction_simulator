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
import argparse

p = argparse.ArgumentParser()
p.add_argument('-input_data', help='Path to generation length (GL) data: first column taxon list, followed by n columns of GL values', default=None, required=True)
args = p.parse_args()


args.input_data = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/gl_data_all_mammals.txt'

# get user input
input_data = args.input_data
taxon_reference_group = 'Mammalia'
n_years = 81
n_rep = 0

data = pd.read_csv(input_data,sep='\t')
species_list = data.iloc[:,0].values
gl_matrix = data.iloc[:,1:].values
if n_rep:
    sim_reps = n_rep
else:
    sim_reps = gl_matrix.shape[1]
    
print('Preparing data for %i simulation replicates'%sim_reps)

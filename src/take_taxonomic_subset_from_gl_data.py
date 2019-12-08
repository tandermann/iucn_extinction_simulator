#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 11:33:40 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)

infile = '/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/gl_data_all_mammals.txt'
all_mammal_data = pd.read_csv(infile,sep='\t',header=None)

phylacine_file = '/Users/tobias/GitHub/iucn_predictions/data/raw/taxonomy/mammals/Trait_data.csv'
phylacine = pd.read_csv(phylacine_file)
phylacine['species_names'] = [' '.join(i.split('_')) for i in phylacine['Binomial.1.2']]
order_info = phylacine['Order.1.2'].values
species_order_dict = dict(zip(phylacine.species_names,order_info))
target_order = 'Carnivora'

subset_all_gl_dta = all_mammal_data[all_mammal_data.iloc[:,0].isin(species_order_dict.keys())].copy()


carnivora_gl_subset = subset_all_gl_dta[[True if species_order_dict[i] ==  target_order else False for i in subset_all_gl_dta.iloc[:,0].values]]
carnivora_gl_subset.to_csv('/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/gl_data_carnivora.txt',header=False,sep='\t',index=False)


carnivora_gl_subset_no_gl_data = carnivora_gl_subset.iloc[:,0]
carnivora_gl_subset_no_gl_data.to_csv('/Users/tobias/GitHub/iucn_extinction_simulator/data/example_data/gl_data_carnivora_no_gl.txt',header=False,sep='\t',index=False)
len(carnivora_gl_subset.shape)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 21:23:42 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""


import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)

# AVES_________________________________________________________________________
in_base = '/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_'
gl_mean = pd.read_csv(in_base+'gl_mean.txt',sep='\t',header=None)
gl_std = pd.read_csv(in_base+'gl_std.txt',sep='\t',header=None)
species_list = pd.read_csv(in_base+'species_list.txt',sep='\t',header=None)[0].values
final_gl_array = np.zeros([100000,len(species_list)])
for i in np.arange(100):
    print(i)
    gl_mean_array = gl_mean[i].values
    gl_std_array = gl_std[i].values
    final_gl_array[i*1000:(i+1)*1000] = np.array([np.random.normal(gl_mean_array,gl_std_array) for i in np.arange(1000)])
random_index_sample = np.random.choice(np.arange(0,final_gl_array.shape[0]),100)
selected_gl_values = final_gl_array[random_index_sample,:].T
final_gl_df = pd.DataFrame(species_list)
final_gl_df.columns = ['species']
for i in np.arange(0,100):
    final_gl_df[i] = np.exp(selected_gl_values[:,i])/365.2422
final_gl_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_final_gl_values.txt',sep='\t',index=False,header=False)



# MAMMALIA_____________________________________________________________________
in_base = '/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/mammals/phylogenetic_imputation/mammals_'
gl_mean = pd.read_csv(in_base+'gl_mean.txt',sep='\t',header=None)
gl_std = pd.read_csv(in_base+'gl_std.txt',sep='\t',header=None)
species_list = pd.read_csv(in_base+'species_list.txt',sep='\t',header=None)[0].values
final_gl_array = np.zeros([100000,len(species_list)])
for i in np.arange(100):
    print(i)
    gl_mean_array = gl_mean[i].values
    gl_std_array = gl_std[i].values
    final_gl_array[i*1000:(i+1)*1000] = np.array([np.random.normal(gl_mean_array,gl_std_array) for i in np.arange(1000)])
random_index_sample = np.random.choice(np.arange(0,final_gl_array.shape[0]),100)
selected_gl_values = final_gl_array[random_index_sample,:].T
final_gl_df = pd.DataFrame(species_list)
final_gl_df.columns = ['species']
for i in np.arange(0,100):
    final_gl_df[i] = np.exp(selected_gl_values[:,i])/365.2422
final_gl_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/mammals/phylogenetic_imputation/mammals_final_gl_values.txt',sep='\t',index=False,header=False)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 17:08:10 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
from ete3 import Tree

# AVES_________________________________________________________________________
# get GL data
gl_data = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/generation_length/birds/birds_bm_gl.txt',sep='\t')

# array of species list, gl, and bm (logtransformed and proper nan values)
species_list = gl_data.species.values
gl_values_log = np.log(gl_data.gl_values_days.values.astype(float))
bm_values_log = np.log(gl_data.bm_values.values.astype(float))
# replace nan values with proper np.nan
gl_values_log[gl_values_log!=gl_values_log] = np.nan
bm_values_log[bm_values_log!=bm_values_log] = np.nan
# jetz species list
jetz_tree = '/Users/tobias/GitHub/iucn_predictions/data/raw/phylogenies/birds/example_tree.tre'
t=Tree(jetz_tree)
species_list_jetz = t.get_leaf_names()
species_list_jetz = [i.replace('_',' ') for i in species_list_jetz]
# get ferran sayol's translation info and create dict
sayol_tax_match_file = '/Users/tobias/GitHub/iucn_predictions/data/raw/taxonomy/birds/TabMatch_HBWBirdLife_Jetz.xlsx'
sayol_tax_match_df = pd.read_excel(sayol_tax_match_file,sheet_name='HBW-BirdLife3_to_Jetz')
iucn_hbw_names = sayol_tax_match_df['species.hbw3'].values
jetz_names = sayol_tax_match_df['species.jetz'].values
iucn_jetz_translation_dict = dict(zip(iucn_hbw_names,jetz_names))
# see how many species match with and without taxonomic correction
matches_raw = [i for i in species_list if i in species_list_jetz]
matches_corrected = []
all_data = []
for i in species_list:
    if i in iucn_jetz_translation_dict.keys():
        if iucn_jetz_translation_dict[i] in species_list_jetz:
            species_name_jetz = iucn_jetz_translation_dict[i]
            bm = bm_values_log[species_list==i]
            gl = gl_values_log[species_list==i]
            matches_corrected.append(i)
            all_data.append([species_name_jetz.replace(' ','_'),bm[0],gl[0]])
    else:
        if i in species_list_jetz:
            species_name_jetz = i
            bm = bm_values_log[species_list==i]
            gl = gl_values_log[species_list==i]
            matches_corrected.append(i)
            all_data.append([species_name_jetz.replace(' ','_'),bm[0],gl[0]])
# before
len(matches_raw)/len(species_list)
# after
len(matches_corrected)/len(species_list)

# create dataframe with all final data, compatible with newest IUCN mammal species list
final_complete_df = pd.DataFrame(data = np.array(all_data),columns=['species','log_bm','log_gl'])
#plt.scatter(final_complete_df.log_bm.values.astype(float),final_complete_df.log_gl.values.astype(float))
# write to file
final_complete_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/birds/trait_data_birds_in_phylogeny.txt',sep='\t',index=False)






# MAMMALIA_____________________________________________________________________
# get GL data
gl_data = pd.read_excel('/Users/tobias/GitHub/iucn_predictions/data/raw/generation_length/mammals/Generation Lenght for Mammals.xlsx')
# set no info values to NaN
gl_data.Calculated_GL_d[gl_data.Calculated_GL_d == 'no information'] = np.nan
gl_values = gl_data.Calculated_GL_d.values.astype(float)
# select those species with real data
subset_with_gl_data = gl_data[gl_values>0].copy()
subset_with_gl_data.Calculated_GL_d = np.log(subset_with_gl_data.Calculated_GL_d.values.astype(float))
# extract species list and GL values
species_list_gl = subset_with_gl_data.Scientific_name.values
gl_data = subset_with_gl_data.Calculated_GL_d.values
# create df
gl_df = pd.DataFrame(data = np.array([species_list_gl,gl_data]).T,columns = ['species','gl'])


# get body mass data
phylacine_traits = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/taxonomy/mammals/Trait_data.csv')
phylacine_species_list = phylacine_traits['Binomial.1.2'].values
phylacine_bm_list = np.log(phylacine_traits['Mass.g'].values)
# set body mass to nan for all species with imputed data in phylacine
imputed_species = phylacine_traits[phylacine_traits['Mass.Method'] == 'Imputed']
phylacine_bm_list[np.isin(phylacine_species_list,imputed_species)] = np.nan

bm_dict = dict(zip(phylacine_species_list,phylacine_bm_list))


# get all available BM and GL info for all species in phylacine
all_data = []
for species in phylacine_species_list:
    if species in bm_dict.keys():
        bm = bm_dict[species]
    else:
        bm = np.nan
    taxon = species.replace('_',' ')
    if taxon in gl_df.species.values:
        gl = gl_df[gl_df.species==taxon].gl.values[0]
    else:
        gl = np.nan
    all_data.append([species,bm,gl])

# create dataframe with all final data, compatible with newest IUCN mammal species list
final_complete_df = pd.DataFrame(data = np.array(all_data),columns=['species','log_bm','log_gl'])

# write to file
final_complete_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/mammals/trait_data_mammals_in_phylogeny.txt',sep='\t',index=False)


plt.plot(final_complete_df.log_bm.values.astype(float),final_complete_df.log_gl.values.astype(float),'ro')







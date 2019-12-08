#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 09:40:37 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)


# AVES_________________________________________________________________________
# get newest IUCN species list
iucn_bird_species_list = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/iucn/species_list_birds.txt',sep='\t',header=None)[0].values

# get array of 100 GL values for each species
gl_data = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/birds/phylogenetic_imputation/birds_final_gl_values.txt',sep='\t',header=None)

# match taxonomy between GL data (jetz taxonomy) with IUCN taxonomy
sayol_tax_match_file = '/Users/tobias/GitHub/iucn_predictions/data/raw/taxonomy/birds/TabMatch_HBWBirdLife_Jetz.xlsx'
sayol_tax_match_df = pd.read_excel(sayol_tax_match_file,sheet_name='HBW-BirdLife3_to_Jetz')
iucn_hbw_names = sayol_tax_match_df['species.hbw3'].values
jetz_names = sayol_tax_match_df['species.jetz'].values
jetz_iucn_translation_dict = dict(zip(jetz_names,iucn_hbw_names))
gl_data['iucn_names'] = [jetz_iucn_translation_dict[name.replace('_',' ')] if name.replace('_',' ') in jetz_iucn_translation_dict.keys() else name.replace('_',' ') for name in gl_data[0].values]
# remove df rows with same iucn name in multiple rows (only keep one)
gl_data = gl_data.iloc[gl_data['iucn_names'].drop_duplicates().index.values,:]
gl_data.index = np.arange(len(gl_data))
# add genus information to df
gl_data['genus'] = [i.split(' ')[0] for i in gl_data['iucn_names'].values]
# get mean GL for each genus
genus_gl_mean = {}
for genus,data in gl_data.groupby(['genus']):
    #print(genus)
    gl_data_genus = data.iloc[:,1:101]
    mean_100_rep = np.mean(gl_data_genus).values
    genus_gl_mean[genus] = mean_100_rep
# get list of taxa for which we don't have gl values yet
final_gl_values_all_species_df = pd.DataFrame(np.zeros((len(iucn_bird_species_list),101)))
final_gl_values_all_species_df.columns = ['species']+ ['gl_years_%i'%i for i in np.arange(0,100)]
final_gl_values_all_species_df.species = iucn_bird_species_list
final_gl_array = np.zeros((len(iucn_bird_species_list),100))
for i,species in enumerate(sorted(iucn_bird_species_list)):
    print(species)
    # if we have gl data for this species, get the 100 values
    if species in gl_data.iucn_names.values:
        gl_values = gl_data[gl_data.iucn_names==species].values[0][1:101].astype(float)
    # if not, get the genus/family/order mean
    else:
        genus = species.split(' ')[0]
        if genus == 'Caloramphus':
            genus = 'Calorhamphus'
        elif genus == 'Chrysocorythus':
            genus = 'Serinus'
        elif genus == 'Cryptopipo':
            genus = 'Xenopipo'
        elif genus == 'Horizocerus':
            genus = 'Tropicranus'
        elif genus == 'Ophrysia':
            genus = 'Rollulus'
        elif genus == 'Poliocrania':
            genus = 'Myrmeciza'
        elif genus == 'Rhamphocharis':
            genus = 'Melanocharis'
        elif genus == 'Rhodonessa':
            genus = 'Anas'
        elif genus == 'Systellura':
            genus = 'Caprimulgus'
        elif genus == 'Thapsinillas':
            genus = 'Alophoixus'
        elif genus == 'Trachylaemus':
            genus = 'Trachyphonus'
        elif genus == 'Tunchiornis':
            genus = 'Hylophilus'
        if genus in genus_gl_mean.keys():
            gl_values = genus_gl_mean[genus]
        else:
            break
    final_gl_array[i,:] = gl_values

# fill data into df
final_gl_values_all_species_df.iloc[:,1:] = final_gl_array
# sort by species name
final_gl_values_all_species_df = final_gl_values_all_species_df.sort_values(by=['species'])
# write the final df to file
final_gl_values_all_species_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/birds/gl_data/gl_data_all_birds.txt',sep='\t',index=False)




# MAMMALIA_____________________________________________________________________
# get newest IUCN species list
iucn_mammal_species_list = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/iucn/species_list_mammals.txt',sep='\t',header=None)[0].values
# get taxonomy information
phylacine_data = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/taxonomy/mammals/Trait_data.csv',sep=',')
# create genus-family and genus-order dict
genus_fam_dict = {}
genus_order_dict = {}
for i in phylacine_data.iterrows():
    genus = i[1]['Genus.1.2']
    family = i[1]['Family.1.2']
    order = i[1]['Order.1.2']
    genus_fam_dict[genus] = family
    genus_order_dict[genus] = order
# get array of 100 GL values for each species
gl_data = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/mammals/phylogenetic_imputation/mammals_final_gl_values.txt',sep='\t',header=None)
# get list of species for which we have GL data
species_with_gl_values = [i.replace('_',' ') for i in gl_data[0].values]
# add genus and family information to df
gl_data['genus'] = [i.split(' ')[0] for i in species_with_gl_values]
gl_data['family'] = [genus_fam_dict[i] for i in gl_data['genus'].values]
gl_data['order'] = [genus_order_dict[i] for i in gl_data['genus'].values]
# get mean GL for each genus
genus_gl_mean = {}
for genus,data in gl_data.groupby(['genus']):
    #print(genus)
    gl_data_genus = data.iloc[:,1:101]
    mean_100_rep = np.mean(gl_data_genus).values
    genus_gl_mean[genus] = mean_100_rep
# get mean GL for each family
family_gl_mean = {}
for family,data in gl_data.groupby(['family']):
    #print(genus)
    gl_data_family = data.iloc[:,1:101]
    mean_100_rep = np.mean(gl_data_family).values
    family_gl_mean[family] = mean_100_rep
# get mean GL for each order
order_gl_mean = {}
for order,data in gl_data.groupby(['order']):
    #print(genus)
    gl_data_order = data.iloc[:,1:101]
    mean_100_rep = np.mean(gl_data_order).values
    order_gl_mean[order] = mean_100_rep
# get list of taxa for which we don't have gl values yet
final_gl_values_all_species_df = pd.DataFrame(np.zeros((len(iucn_mammal_species_list),101)))
final_gl_values_all_species_df.columns = ['species']+ ['gl_years_%i'%i for i in np.arange(0,100)]
final_gl_values_all_species_df.species = iucn_mammal_species_list
final_gl_array = np.zeros((len(iucn_mammal_species_list),100))
for i,species in enumerate(iucn_mammal_species_list):
    # if we have gl data for this species, get the 100 values
    if species in [i.replace('_',' ') for i in gl_data[0].values]:
        gl_values = gl_data[gl_data[0]==species.replace(' ','_')].values[0][1:101].astype(float)
    # if not, get the genus/family/order mean
    else:
        genus = species.split(' ')[0]
        if genus == 'Micaelamys':
            genus = 'Aethomys'
        if genus == 'Plecturocebus':
            genus = 'Callicebus'
        if genus == 'Brassomys':
            genus = 'Melomys'
        if genus == 'Gyldenstolpia':
            genus = 'Kunsia'
        if genus == 'Gardnerycteris':
            genus = 'Mimon'
        if genus == 'Nannospalax':
            genus = 'Spalax'
        if genus == 'Parastrellus' or 'Perimyotis':    
            genus = 'Pipistrellus'
        if genus in genus_gl_mean.keys():
            gl_values = genus_gl_mean[genus]
        else:
            print(species)
            break
#            family = genus_fam_dict[genus]
#            if family in family_gl_mean.keys():
#                gl_values = family_gl_mean[family]
#            else:
#                if genus in order_gl_mean.keys():
#                    gl_values = order_gl_mean[genus]
#                else:
#                    print(species)
#                    break   
    final_gl_array[i,:] = gl_values
# fill data into df
final_gl_values_all_species_df.iloc[:,1:] = final_gl_array
# sort by species name
final_gl_values_all_species_df = final_gl_values_all_species_df.sort_values(by=['species'])
# write the final df to file
final_gl_values_all_species_df.to_csv('/Users/tobias/GitHub/iucn_predictions/data/processed/trait_data/mammals/gl_data/gl_data_all_mammals.txt',sep='\t',index=False)
    










#set([i.split(' ')[0] for i in iucn_mammal_species_list])-set([i.split(' ')[0] for i in species_with_gl_values])

# if genus is present in gl data, take genus mean:



# if body mass available to genus/family/order regression

# else take genus/family/order mean


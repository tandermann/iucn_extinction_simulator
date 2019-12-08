#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 10:13:44 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)
import urllib3
from bs4 import BeautifulSoup
import requests
import re
from selenium import webdriver
import time

taxon_ids = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/iucn/mammals_taxon_ids.txt',sep='\t',header=None)[0].values
species_names = pd.read_csv('/Users/tobias/GitHub/iucn_predictions/data/raw/iucn/species_list_mammals.txt',sep='\t',header=None)[0].values
species_id_dict = dict(zip(species_names,taxon_ids))

species_url_dict = {}
for species in species_id_dict.keys():
    taxon_id = species_id_dict[species]
    
    driver=webdriver.Chrome('/Users/tobias/bin/chromedriver')
    url = 'https://www.iucnredlist.org/search/grid?query=%s&searchType=species'%species
    driver.get(url)   
    driver.find_elements_by_class_name('layout-card--split__minor')
    html=driver.page_source

    target = 'href="/species/%i/'%taxon_id
    left,sep,right = html.partition(target)
    time.sleep(5)
    if sep:
        random_int = right.split('"')[0]
        target_url = 'https://www.iucnredlist.org/species/%s/%s'%(str(taxon_id),str(random_int))
        species_url_dict[species] = target_url
    else:
        print('species not found:',species)
        break


for species in species_url_dict.keys():
#    driver.get('https://www.iucnredlist.org/species/3129/46364616')
    driver=webdriver.Chrome('/Users/tobias/bin/chromedriver')
    url = species_url_dict[species]
    driver.get(url)   
    driver.find_elements_by_class_name('layout-card--split__minor')
    html=driver.page_source
    left,sep,right = html.partition('<h3 class="heading">Generation length (years)</h3><p class="card__data card__data--std card__data--accent">')
    if sep:
        generation_time = right.split('</p>')[0]
        print(species, generation_time)


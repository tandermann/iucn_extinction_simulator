#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 20:59:28 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)


    # replace status with index
    final_change_types, final_years_count = cust_func.replace_iucn_status_with_int(change_type_dict_no_dd,years_in_each_category)
    
    # generate q-matrices__________________________________________________________
    qmatrix_dict_list = []
    for i in np.arange(100):
        print(i)
        en_risks_rep = en_risks.T[i]
        cr_risks_rep = cr_risks.T[i]
        q_matrix_dict = {}
        for j,species in enumerate(species_list):
            en_risk = en_risks_rep[j]
            cr_risk = cr_risks_rep[j]
            status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,en_risk,cr_risk]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category
            q_matrix = cust_func.qmatrix(final_change_types,final_years_count, status_specific_p_e)
            q_matrix_dict[species] = q_matrix
        qmatrix_dict_list.append(q_matrix_dict)
    
    # also calculate regular q-matrix, not accounting for GL
    status_specific_p_e = np.array([0.000000155728,0.000041551152,0.001053050310,0.011095167095,0.066967008463]) # These values are the category specific probabilities of extinction per year calculated from IUCN definition of each category
    q_matrix = cust_func.qmatrix(final_change_types,final_years_count, status_specific_p_e)

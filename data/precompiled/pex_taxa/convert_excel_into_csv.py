#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 14:22:49 2020

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import pandas as pd

file = './2020_1_RL_Stats_Table_9.xlsx'
data_df = pd.read_excel(file,header=0)
data_df.to_csv(file.replace('.xlsx','.txt'),sep='\t',index=False)

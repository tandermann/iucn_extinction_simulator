#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:30:52 2019

@author: Tobias Andermann (tobias.andermann@bioenv.gu.se)
"""

import numpy as np
np.set_printoptions(suppress=True)
import pandas as pd
import matplotlib.pyplot as plt
np.random.seed(1234)



def sample_rate(count, tot_time, n_samples = 1, range_factor = 100, n_bins = 10000):
    def get_loglik(count, dT, rate):
        return np.log(rate)*count - dT*rate
    mle = count/tot_time
    if count == 0:
        minRange = 1*np.exp(-10)
        maxRange = range_factor/tot_time
    else:
        minRange = mle/range_factor
        maxRange = mle*range_factor
    rates = np.linspace(minRange, maxRange,n_bins)
    lik = get_loglik(count,tot_time,rates)
    lik = lik - np.max(lik)
    sample_rates = np.random.choice(rates,size=n_samples,p=np.exp(lik)/np.sum(np.exp(lik)),replace=1)
    #plt.plot(np.log(rates),lik)
    #np.log(mle)
    return sample_rates



max_t = 100
x = np.array([33,12,45,100,100,54,100,33,44,6])
n_bins = 10000
n_samples=100

q_samples = np.linspace(0, 1, n_bins)
all_sampled_rates = []
for x_value in x:
    if x_value < max_t:
        #extinct        
        p = q_samples*np.exp(-q_samples*x_value)
    else:
        # not extinct
        p = np.exp(-q_samples*x_value)
    p = p/sum(p)
    sample_rates = np.random.choice(q_samples,size=n_samples,p=p,replace=1)
    all_sampled_rates.append(sample_rates)
all_rates = [item for sublist in all_sampled_rates for item in sublist]
np.mean(all_rates)





# 10000 q-s from uniform between 0 and 1
# calculate posterior for all q-s and each given x-value
# resample 100 q-values with probability=posterior --> needs to be normalized, summing up to 1, replace=True






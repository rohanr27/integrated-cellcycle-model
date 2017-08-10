#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 10:58:51 2017

@author: rohanroy
"""


from pandas import *
from cobra.flux_analysis.sampling import OptGPSampler

sampler1 = OptGPSampler(dnamodel1, processes=2, thinning=200)
sampler2 = OptGPSampler(dnamodel2, processes=2, thinning=200)

dnamodel1samples = sampler1.sample(4500)
dnamodel2samples = sampler2.sample(4500)

#%%

for reaction in dnamodel1.reactions:
    reaction.averageflux=dnamodel1samples[reaction.id].mean()
    reaction.stdflux=dnamodel1samples[reaction.id].std()
    reaction.medianflux=dnamodel1samples[reaction.id].median()
    
for reaction in dnamodel2.reactions:
    reaction.averageflux = dnamodel2samples[reaction.id].mean()
    reaction.stdflux=dnamodel2samples[reaction.id].std()
    reaction.medianflux=dnamodel2samples[reaction.id].median()
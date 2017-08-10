# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:51:58 2017

@author: anwara
"""


import cobra
from pylab import *
import numpy

model = cobra.io.read_sbml_model("yeast_7.00_cobra.xml")
#dnamodel=model.copy()
sol=model.optimize().objective_value

#%% DNA replication Model 1: Create New Reaction for the Synthesis of DNA and maximize alongside growth

from cobra import Model, Reaction, Metabolite


all_biomass = iter(model.metabolites.s_0450.reactions)
r = [' ']*(len(model.metabolites.s_0450.reactions))
for i in range(len(model.metabolites.s_0450.reactions)):
    r[i] = str(next(all_biomass))
    
biomass1 = iter(model.reactions.get_by_id(r[0]).metabolites)
reac1 = model.reactions.get_by_id(r[0]).reaction
r1met= [' ']*(len(model.reactions.get_by_id(r[0]).metabolites))
for x in range(len(model.reactions.get_by_id(r[0]).metabolites)):
    a = model.metabolites.get_by_id(str(next(biomass1)))
    r1met[x] = [a.id, a.name]
    
biomass2 = iter(model.reactions.get_by_id(r[1]).metabolites)
reac2 = model.reactions.get_by_id(r[1]).reaction
r2met= [' ']*(len(model.reactions.get_by_id(r[1]).metabolites))
for x in range(len(model.reactions.get_by_id(r[1]).metabolites)):
    a = model.metabolites.get_by_id(str(next(biomass2)))
    r2met[x] = [a.id, a.name]
    
biomass3 = iter(model.reactions.get_by_id(r[2]).metabolites)
reac3 = model.reactions.get_by_id(r[2]).reaction
r3met= [' ']*(len(model.reactions.get_by_id(r[2]).metabolites))
for x in range(len(model.reactions.get_by_id(r[2]).metabolites)):
    a = model.metabolites.get_by_id(str(next(biomass3)))
    r3met[x] = [a.id, a.name]
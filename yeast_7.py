#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 16:58:01 2017

@author: rohanroy
"""

import cobra
from pylab import *

model = cobra.io.read_sbml_model("yeast_7.00_cobra.xml")
#dnamodel=model.copy()
sol=model.optimize().objective_value

names_met=[]
names_rxn=[]
for metabolite in model.metabolites:
    names_met.append(metabolite.name)

for reaction in model.reactions:
    names_rxn.append(reaction.name)


#%% DNA replication Model 1: Create New Reaction for the Synthesis of DNA and maximize alongside growth

from cobra import Model, Reaction, Metabolite
dnamodel=model.copy()
# reaction for the synthesis of DNA: dATP + dTTp + dGTP + dCTP --> DNA + PP
reaction = Reaction('DNAsynth')
reaction.name = 'DNA synthesis'
reaction.lower_bound=0
reaction.upper_bound=1000
s_0586=dnamodel.metabolites.s_0586 #dATP
s_0650=dnamodel.metabolites.s_0650 #dTTP
s_0617=dnamodel.metabolites.s_0617 #dGTP
s_0590=dnamodel.metabolites.s_0590 #dCTP
s_0434=dnamodel.metabolites.s_0434 #ATP
s_0425=dnamodel.metabolites.s_0425 #AMP
s_0637=dnamodel.metabolites.s_0637 #diphosphate
DNA = Metabolite('DNA',name='DNA')
reaction.add_metabolites({s_0586: -0.0036, s_0650:-0.0036, s_0617:-0.0024, s_0590:-0.0024, DNA: 0.012, s_0637: 0.012})
dnamodel.add_reactions([reaction])

# exchange reaction for DNA: DNA -->
reaction1=Reaction('DNA')
reaction1.name='DNA replication'
reaction1.lower_bound=0
reaction1.upper_bound=1000
reaction1.add_metabolites({DNA:-1})
dnamodel.add_reactions([reaction1])

# transport reaction for muclear pyrophosphate: PP [nucleus] <--> PP [cytoplasm]
reaction2=Reaction('r_5000')
reaction2.name='Diphosphate transport, nucleus-cytoplasm'
reaction2.lower_bound = -1000
reaction2.upper_bound = 1000
reaction2.add_metabolites({model.metabolites.s_0637:-1,model.metabolites.s_0633:1})
dnamodel.add_reactions([reaction2]) 

# remove dNMP from growth expression: dAMP (s_0584), dCMP (s_0589), dGMP (s_0615), dTMP (s_0649)
growth=dnamodel.reactions.r_2133
growth.subtract_metabolites({dnamodel.metabolites.s_0584:-0.0036, dnamodel.metabolites.s_0649:-0.0036, dnamodel.metabolites.s_0589: -0.0024, dnamodel.metabolites.s_0615: -0.0024})

# change coefficients of objective function to reflect need for DNA synthesis
c = linspace(0.19, 0.24, 21)

for i in c:
    dnamodel.objective = {dnamodel.reactions.r_2111:1, dnamodel.reactions.DNA:i}
    print dnamodel.objective.expression
    dnamodel.optimize()
    print dnamodel.summary()
    print dnamodel.reactions.r_0976.flux

#%% DNA replication Model 2: Change Coefficients of existing biomass objective function

# Levels of all four dNTPS increase 2.5-4 fold as cells enter S phase
change = 4
dnamodel=model.copy()
growth = dnamodel.reactions.r_2133

# add dAMP (s_0584), dCMP (s_0589), dGMP (s_0615), dTMP (s_0649)
growth.add_metabolites({dnamodel.metabolites.s_0584:-0.0036*change,
                        dnamodel.metabolites.s_0649:-0.0036*change,
                        dnamodel.metabolites.s_0589: -0.0024*change, 
                        dnamodel.metabolites.s_0615: -0.0024*change})
dnamodel.optimize()

#%% DNA replication Model 3: Remove dNMP from biomass objective function and replace with dNTP

dnamodel=model.copy()
growth=dnamodel.reactions.r_2133

#Subtract dNMP
growth.subtract_metabolites({dnamodel.metabolites.s_0584:-0.0036, dnamodel.metabolites.s_0649:-0.0036, dnamodel.metabolites.s_0589: -0.0024, dnamodel.metabolites.s_0615: -0.0024})

# add dATP (s_0586), dCTP (s_0590), dGTP (s_0617), dTTP (s_0650) to the biomass pseudoreaction
growth.add_metabolites({model.metabolites.s_0586:-0.0036, model.metabolites.s_0590:-0.0024, model.metabolites.s_0617:-0.0024, model.metabolites.s_0650:-0.0036})
#Add additional phosphate [cytoplasm] (s_1322)
growth.add_metabolites({dnamodel.metabolites.s_1322: 0.012*2})

dnamodel.optimize()

#%% DNA replication Model 4: Minimizing difference between dNTP consumption and dNTP production


from cobra import Model, Reaction, Metabolite

dnamodel=model.copy()

#dNTP consumption from S phase
vc=8

#reactions that produce dNTP:
#r_0970 (dATP)
#r_0796 (dATP)
#r_0971 (dCTP)
#r_0797 (dCTP)
#r_0972 (dGTP)
#r_0798 (dGTP)
#r_1750 (dTTP) ??
#r_0799 (dTTP)


import cplex
import optlang
from optlang import cplex_interface

dnamodel.solver = optlang.cplex_interface

same_flux_purine = dnamodel.problem.Constraint(
        dnamodel.reactions.r_0970.flux_expression
        +dnamodel.reactions.r_0796.flux_expression
        -dnamodel.reactions.r_0799.flux_expression,lb=0,ub=0)
dnamodel.add_cons_vars(same_flux_purine)

same_flux_pyrimidine = dnamodel.problem.Constraint(
        dnamodel.reactions.r_0971.flux_expression
        +dnamodel.reactions.r_0797.flux_expression
        -dnamodel.reactions.r_0972.flux_expression
        -dnamodel.reactions.r_0798.flux_expression,
        lb=0,ub=0)
dnamodel.add_cons_vars(same_flux_pyrimidine)

#dna_composition = dnamodel.problem.Constraint(
#        (dnamodel.reactions.r_0970.flux_expression+dnamodel.reactions.r_0796.flux_expression)
#        /(dnamodel.reactions.r_0972.flux_expression+dnamodel.reactions.r_0798.flux_expression), lb=1.632, ub=1.632)
#dnamodel.add_cons_vars(dna_composition)


obj= dnamodel.problem.Objective(
    (dnamodel.reactions.r_0970.flux_expression+dnamodel.reactions.r_0796.flux_expression
    +dnamodel.reactions.r_0971.flux_expression+dnamodel.reactions.r_0797.flux_expression
    +dnamodel.reactions.r_0972.flux_expression+dnamodel.reactions.r_0798.flux_expression
    +dnamodel.reactions.r_0799.flux_expression - vc)**2, direction = 'min')



#dnamodel.solver = 'glpk'
#obj= dnamodel.problem.Objective(
#    dnamodel.reactions.r_0970.flux_expression+dnamodel.reactions.r_0796.flux_expression
#    +dnamodel.reactions.r_0971.flux_expression+dnamodel.reactions.r_0797.flux_expression
#    +dnamodel.reactions.r_0972.flux_expression+dnamodel.reactions.r_0798.flux_expression
#    +dnamodel.reactions.r_0799.flux_expression - vc, direction = 'min')

dnamodel.objective = obj






#%% Functions


# write reaction in terms of species name instead of IDs
def rxnwithnames(model,reaction, line = False):
    if line is True:
        rxn=reaction.reaction
        name=''
        for i in range (0, len(rxn)-5):
            if rxn[i]+rxn[i+1] == 's_':
                name = name + model.metabolites.get_by_id(rxn[i:i+6]).name
            name = name + rxn[i]
            if len(name)>6:
                if name[len(name)-6] + name[len(name)-5] == 's_':
                    name = name[0:len(name)-6]    
        name = name[0:len(name)-1]    
        return name
    if line is False:
        y={}
        for x in reaction.metabolites:
            y[x.name]= reaction.metabolites[x]
        return y
        
# yeast 6 biomass psuedoreaction

def biomass6(model):
    for key in model.reactions.get_by_id("r_2133").metabolites:
        print key.name, model.reactions.get_by_id("r_2133").metabolites[key] 
        


# yeast 5 biomass pseudoreaction
def biomass5(model):
    for key in model.reactions.get_by_id("r_2110").metabolites:
        print key.name, model.reactions.get_by_id("r_2110").metabolites[key]

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
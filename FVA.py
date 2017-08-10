#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:11:26 2017

@author: rohanroy
"""

#Flux Variability Analysis

import cobra
import optlang
from optlang import gurobi_interface
from cobra import Model, Reaction, Metabolite

model = cobra.io.read_sbml_model("yeast_7.00_cobra.xml")
dnamodel=model.copy()

#dNTP consumption from S phase
#vc=8

#reactions that produce dNTP:
#r_0970 (dATP)
#r_0796 (dATP)
#r_0529 (dATP)
#r_0971 (dCTP)
#r_0797 (dCTP)
#r_0972 (dGTP)
#r_0798 (dGTP)
#r_1750 (dTTP) 
#r_0799 (dTTP)

dntp = [dnamodel.reactions.r_0970, dnamodel.reactions.r_0796, dnamodel.reactions.r_0529,
        dnamodel.reactions.r_0971, dnamodel.reactions.r_0797, dnamodel.reactions.r_0972,
        dnamodel.reactions.r_0798, dnamodel.reactions.r_0799]

glycolysis = [dnamodel.reactions.r_0467, dnamodel.reactions.r_0892, dnamodel.reactions.r_0962,
              dnamodel.reactions.r_0366, dnamodel.reactions.r_0713]

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


dnamodel.solver = optlang.cplex_interface
vc1 = 0.7879604847584999
fbasolution = -0.13904901487262922

'''
fba = dnamodel.problem.Constraint((dnamodel.reactions.DNAsynth.flux_expression - vc1)**2-dnamodel.reactions.r_2111.flux_expression, lb=fbasolution,ub=fbasolution)
dnamodel.add_cons_vars(fba)

minflux = {}
maxflux = {}
for reaction in dnamodel.reactions:
    objmin = dnamodel.problem.Objective(reaction.flux_expression, direction = 'min')
    dnamodel.objective = objmin
    dnamodel.optimize()
    minflux[reaction.id]=dnamodel.optimize().objective_value
    
    objmax = dnamodel.problem.Objective(reaction.flux_expression, direction = 'max')
    dnamodel.objetive = objmax
    dnamodel.optimize()
    maxflux[reaction.id]=dnamodel.optimize().objective_value
'''
#%% dnamodel 1
dnamodel1 = dnamodel.copy()
biomassfba = 0.13905039909452585
vsfba = 0.7867839551600597

biomass = dnamodel1.problem.Constraint(dnamodel1.reactions.r_2111.flux_expression, lb = biomassfba, ub = biomassfba)
dnamodel1.add_cons_vars(biomass)

dnasynth = dnamodel1.problem.Constraint(dnamodel1.reactions.DNAsynth.flux_expression, lb=vsfba, ub=vsfba)
dnamodel1.add_cons_vars(dnasynth)


minflux1 = {}
maxflux1 = {}
for reaction in dnamodel1.reactions:
    objmin = dnamodel1.problem.Objective(reaction.flux_expression, direction = 'min')
    dnamodel1.objective = objmin
    dnamodel1.optimize()
    minflux1[reaction.id]=dnamodel1.optimize().objective_value
    
    objmax = dnamodel1.problem.Objective(reaction.flux_expression, direction = 'max')
    dnamodel1.objective = objmax
    dnamodel1.optimize()
    maxflux1[reaction.id]=dnamodel1.optimize().objective_value

diff1 = {}
for x in maxflux1:
    diff1[x]=maxflux1[x]-minflux1[x]


#%% dnamodel 2
dnamodel2 = dnamodel.copy()
biomassfba = 0.14090174871323874
vsfba = 1.4090174831616546e-05

biomass = dnamodel2.problem.Constraint(dnamodel2.reactions.r_2111.flux_expression, lb = biomassfba, ub = biomassfba)
dnamodel2.add_cons_vars(biomass)

dnasynth = dnamodel2.problem.Constraint(dnamodel2.reactions.DNAsynth.flux_expression, lb=vsfba, ub=vsfba)
dnamodel2.add_cons_vars(dnasynth)


minflux2 = {}
maxflux2 = {}
for reaction in dnamodel2.reactions:
    objmin = dnamodel2.problem.Objective(reaction.flux_expression, direction = 'min')
    dnamodel2.objective = objmin
    dnamodel2.optimize()
    minflux2[reaction.id]=dnamodel2.optimize().objective_value
    
    objmax = dnamodel2.problem.Objective(reaction.flux_expression, direction = 'max')
    dnamodel2.objective = objmax
    dnamodel2.optimize()
    maxflux2[reaction.id]=dnamodel2.optimize().objective_value

diff2 = {}
for x in maxflux2:
    diff2[x]=maxflux2[x]-minflux2[x]




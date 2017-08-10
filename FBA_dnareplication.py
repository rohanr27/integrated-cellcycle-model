#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 10:11:08 2017

@author: rohanroy
"""


import cobra
from pylab import *

model = cobra.io.read_sbml_model("yeast_7.00_cobra.xml")
sol=model.optimize().objective_value

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

#save all reaction fluxes from optimized model to text file
def fluxfile(dnamodel,name):
    file_ = open(name,"w")
    for reactions in dnamodel.reactions:
        if reactions.flux != 0:
            file_.write('\n' + reactions.id + ' ' + reactions.name + '\n' 
                        + rxnwithnames(dnamodel, reactions, line= True)
                        + '\n' + str(reactions.averageflux) + '\n')        


from cobra import Model, Reaction, Metabolite

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

#%%

import cplex
import optlang
from optlang import cplex_interface
dnamodel.solver = optlang.cplex_interface


genome = 1.21e7         # length of S. Cerevisiae genome (bp)
lengthS = 17.0/60.0    # experimental length of S phase (hour)
cellmass = 15.0         # dry mass of cell (pg)
stoich = dnamodel.reactions.DNAsynth.metabolites[dnamodel.metabolites.DNA]  #stoichiometric coefficient of DNA (sum of dNTP) in DNAsynth reaction
k = 2*1e3*1e12/6.022e23

#approximate rate of DNA consumption
vc1 = (k*genome)/(lengthS*cellmass*stoich)

dnamodel1=dnamodel.copy()
obj1= dnamodel1.problem.Objective(
            (dnamodel1.reactions.DNAsynth.flux_expression - vc1)**2-dnamodel1.reactions.r_2111.flux_expression, direction = 'min')
dnamodel1.objective = obj1
dnamodel1.optimize()       
#fluxfile(dnamodel1,'vc_'+str(vc1)+'_fbamodel')

#predicted length of S phase 
#lengthSpred = (k*genome)/(cellmass*dnamodel1.reactions.DNAsynth.flux*stoich)
#diff=100*(abs(lengthSpred-lengthS))/lengthS


#%%

#dnamodel2, growth (without dna consumption) objective function. dna consumption rate set to zero
vc2 = 0
dnamodel2=dnamodel.copy() 
obj2= dnamodel2.problem.Objective(
            (dnamodel2.reactions.DNAsynth.flux_expression-vc2)**2-dnamodel2.reactions.r_2111.flux_expression, direction = 'min')
dnamodel2.objective = obj2
dnamodel2.optimize()       

#fluxfile(dnamodel2,'vc_'+str(vc2)+'_fbamodel')
      

#%%
'''
#sort dnamodel1 by difference (dnamodel1 reaction flux - dnamodel2 reaction flux) for each reaction
x = len(dnamodel1.reactions)
for i in range(0,x):
    dnamodel1.reactions[i].diff = dnamodel1.reactions[i].flux-dnamodel2.reactions[i].flux

dnamodel1sort=sorted(dnamodel1.reactions, key=lambda reaction:reaction.diff, reverse=True)


file_ = open('flux_differences_sort','w')
file_.write('DNA Replication Model Flux - Growth Model Flux \n')
for i in dnamodel1sort:
    file_.write('\n' + i.id + ' ' + i.name + '\n'
                + rxnwithnames(dnamodel1, i,line=True)+ '\n'
                + str(i.diff)+'\n')
'''




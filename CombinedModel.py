#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:06:07 2017

@author: rohanroy
"""
from pylab import *
import cobra
import roadrunner




#Combined FBA-Cell Cycle Model for a single Saccharomyces cerevisiae Cell

#Number of cells to simulate
ncells = 3000

#%% Cell Cycle Simulation

#Assign Initial Conditions of stochastic cell cycle model from a single cell cycle simulation of the deterministic model
tyson = roadrunner.RoadRunner("test.xml")
varIC = roadrunner.RoadRunner("test.xml")

tyson.model['init(m)']=0.45
tyson.model['mu']=0.008
tyson.timeCourseSelections=['time','Cdc20a','Cdc20t','Cdh1','CKIt','m','CycBt','IEP','SK']
result = tyson.simulate(0,94,95)

#Initial Condition is the last row of values of simulation
Initial_Conditions = result[94,:]
for i in range(1, len(Initial_Conditions)):
    varIC.model['init('+tyson.timeCourseSelections[i]+')']=Initial_Conditions[i]

 
varIC.timeCourseSelections=['time','CycB','Cdc20a','Cdh1','m','CKIt','SK']
#Growth Rate of Cells
varIC.model['mu']=0.008


G1length=[]
S_G2_Mlength=[]
broke=0
samplemass = []
for i in range(0,ncells):
    varIC.reset()
    for species in varIC.getIndependentFloatingSpeciesIds():
        if species == 'm':
            division = np.random.normal(0.5,0.05)
            varIC.model['init(m)']= division*0.90
        else:
            varIC.model['init('+species+')']=np.random.normal(varIC.model['init('+species+')'],varIC.model['init('+species+')']*0.01)         
             #Variable Initial Conditions drawn from Normal distribution, mu = reported value, sigma = mu*0.01

    #Vary the cell mass @ division by changing the time-lag associated with IEP (CycB -> activation of cdc20)
    varIC.model['k9'] = np.random.normal(0.1,0.01)
    result = varIC.simulate(0,150,151)
    
    #Differentiate Cdh1 concentration
    times = result[:,0]
    dCdh1dt = diff(result[:,3])/diff(times)
    
    samplemass.append(result[np.random.randint(0,120),4])

    startS =0
    endM =0
    count =0
    while startS == 0:
        if dCdh1dt[count]<-0.20:
            startS = count /60.0
        count += 1
    count =0
    
    while endM ==0 and count < len(dCdh1dt)-1:
        if dCdh1dt[count]>0.40:
            endM = count/60.0
        count +=1
    
    if startS>0 and endM >0 and endM > startS:
        G1length.append(startS)
        S_G2_Mlength.append(endM - startS)
    else:
        broke+=1

#%% Building FBA model

from cobra import Model, Reaction, Metabolite
model = cobra.io.read_sbml_model("yeast_7.00_cobra.xml")
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

import cplex
import optlang
from optlang import cplex_interface
dnamodel.solver = optlang.cplex_interface


genome = 1.21e7         # length of S. Cerevisiae genome (bp)
lengthS = 17.0/60.0    # experimental length of S phase (hour)
cellmass = 50.0         # mass of cell (pg)
stoich = dnamodel.reactions.DNAsynth.metabolites[dnamodel.metabolites.DNA]  #stoichiometric coefficient of DNA (sum of dNTP) in DNAsynth reaction
k = 2*1e3*1e12/6.022e23/0.30

#approximate rate of DNA consumption
vc1 = (k*genome)/(lengthS*cellmass*stoich)

dnamodel1=dnamodel.copy()
obj1= dnamodel1.problem.Objective(
            (dnamodel1.reactions.DNAsynth.flux_expression - vc1)**2-dnamodel1.reactions.r_2111.flux_expression, direction = 'min')
dnamodel1.objective = obj1
    

vc2 = 0
dnamodel2=dnamodel.copy() 
obj2= dnamodel2.problem.Objective(
            (dnamodel2.reactions.DNAsynth.flux_expression-vc2)**2-dnamodel2.reactions.r_2111.flux_expression, direction = 'min')
dnamodel2.objective = obj2       

#%% Flux Sampling

from pandas import *
from cobra.flux_analysis.sampling import OptGPSampler

sampler1 = OptGPSampler(dnamodel1, processes=2, thinning=200)
sampler2 = OptGPSampler(dnamodel2, processes=2, thinning=200)

dnamodel1samples = sampler1.sample(ncells)
dnamodel2samples = sampler2.sample(ncells)

for reaction in dnamodel1.reactions:
    reaction.averageflux=dnamodel1samples[reaction.id].mean()
    reaction.stdflux=dnamodel1samples[reaction.id].std()
       
for reaction in dnamodel2.reactions:
    reaction.averageflux = dnamodel2samples[reaction.id].mean()
    reaction.stdflux=dnamodel2samples[reaction.id].std()


#%% ratio   


#ratio of dna to protein is relatively constant

datp=[0]*(len(G1length))
arg=[0]*(len(G1length))
ratio=[0]*(len(G1length))

for l in range(0,len(G1length)):
    
    lengthSpred = (k*genome)/(cellmass*dnamodel1samples['DNAsynth'][l]*stoich)


    datp[l] = 0.0036*(dnamodel2samples['DNAsynth'][l]*G1length[l]
                     +dnamodel1samples['DNAsynth'][l]*lengthSpred
                     +dnamodel2samples['DNAsynth'][l]*(S_G2_Mlength[l]-lengthSpred))
    
    arg[l] = 0.1607*(dnamodel2samples['r_2133'][l]*G1length[l]
                    +dnamodel1samples['r_2133'][l]*lengthSpred
                    +dnamodel2samples['r_2133'][l]*(S_G2_Mlength[l]-lengthSpred))
    
    ratio[l] = arg[l]/datp[l]




#%% Figures

# G1 Length Distribution
figure()
hist(G1length,20)
xlabel('time(hr)')
ylabel('count')
title('G1 Duration')
xlim(0.3,1.8) 

# S/G2/M Length Distribution
figure()
hist(S_G2_Mlength,20)
xlabel('time(hr)')
ylabel('count')
title('S_G2_M Duration')
xlim(0.3,1.8)              
            
# G1 length vs. S/G2/M length scatterplot
figure()
scatter(G1length,S_G2_Mlength)            
xlim(0.4,1.8)            
ylim(0.4,1.8)
xlabel('G1 Duration (h)')
ylabel('S/G2/M Duration (h)')

#Cell Mass Distribution
s = []*len(samplemass)
for i in samplemass:
    s.append(i*100)
figure()
bins = bincount(s, minlength= 101)
prob = bins/ float(bins.sum())
perc = []*len(prob)
for x in prob:
    perc.append(x*100)
bar(arange(len(bins)),perc,width=1)
xlim([0,120])
xlabel('cell mass')
ylabel('Number %')    
    
#Ratio of arginine (protein) to datp (nucleotide)
figure()
hist(ratio,20)        
plt.axvline(x=average(ratio),color='r',linewidth='2.0')            
#plt.axvline(x=0.1607/0.0036,color='g',linewidth='2.0')            
xlabel('arg/datp')
ylabel('count')            
            
















#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:08:17 2017

@author: rohanroy
"""

#Variable Cell Cycle Model

import roadrunner
import cobra
from pylab import *


'''
from libsbml import *
reader = SBMLReader()
document = reader.readSBML("BIOMD0000000195.xml")
document.getNumErrors()
tyson=document.getModel()
division = tyson.createParameter()
division.setId('division')
division.setConstant(False)
division.setValue(0.5)
celldivision=tyson.getEvent(0)
new_mathml=parseFormula('m * division')
celldivision.getEventAssignment(0).setMath(new_mathml)
writeSBML(document, 'test.xml')
'''


#Calibration of Stochastic Cell Cycle Model

#Assign Initial Conditions of stochastic cell cycle model from the end points of a single cell cycle simulation of the deterministic model

tyson = roadrunner.RoadRunner("BIOMD0000000195.xml")
varIC = roadrunner.RoadRunner("BIOMD0000000195.xml")

#Initialize starting mass to 0.45 
tyson.model['init(m)']=0.45

#Initialize growth rate to 0.008
tyson.model['mu']=0.008

#Select all 8 dependendent species to track
tyson.timeCourseSelections=['time','Cdc20a','Cdc20t','Cdh1','CKIt','m','CycBt','IEP','SK']

#Simulate ODE model for 94 minutes (exactly one cell cycle)
result = tyson.simulate(0,300,300)
times = result[:,0]

#Plot ODE results

for i in range(1, len(tyson.timeCourseSelections)):
    series=result[:,i]
    name = tyson.timeCourseSelections[i]
    plot(times, series, label = str(name))
    #legend(loc='best')


#Initial Condition is the last row of values of simulation
Initial_Conditions = result[94,:]
for i in range(1, len(Initial_Conditions)):
    tyson.model['init('+tyson.timeCourseSelections[i]+')']=Initial_Conditions[i]


#%%
#Running the Stochastic Cell Cycle Model

sim=1000            #number of simulations
varIC.timeCourseSelections=['time','CycB','Cdc20a','Cdh1','m','CKIt','SK']
G1length=[]         #length of G1 Phase
S_G2_Mlength=[]     #length of combined S/G2/M phse
broke=0             #Broke simulations have uncontrolled growth without division
samplemass = []     #'Relative Mass' of cell randomly sampled at any point in the cell cycle
massatS = []        #'Relative Mass' of cell at the beginning of S phase
begmass = []        #'Relative Mass' of cell at beginning of cell cycle
massatG2 = []       #'Relative Mass' of cell at the beginning of G2 phase

for i in range(0,sim):
    varIC.reset()
    for species in varIC.getIndependentFloatingSpeciesIds():
        # Each cell has a different size at birth due to assymetric division
        if species == 'm':
            division = np.random.normal(0.5,0.05)
            varIC.model['init(m)']= division*0.90
        #Variable Initial Conditions drawn from Normal distribution, mu = reported value, sigma = mu*0.01
        else:
            varIC.model['init('+species+')']=np.random.normal(tyson.model['init('+species+')'],tyson.model['init('+species+')']*0.01)         
             
             
    #Vary the cell mass @ division by changing the time-lag associated with IEP (CycB -> activation of cdc20)
    varIC.model['k9'] = np.random.normal(0.1,0.01)
    
    #Default Growth Rate
    varIC.model['mu']=0.008
    result = varIC.simulate(0,250,251)
    times = result[:,0]
    
    '''
    #Plot of stochastic simulation
    figure()
    for j in range(1, len(varIC.timeCourseSelections)):
        series=result[:,j]
        name = varIC.timeCourseSelections[j]
        plot(times, series, label = str(name))
        #legend(loc='best')
        xlabel('Time (min)')
        ylabel('concentration (1/l)')
        title('Stochastic Cell Cycle')
    '''
    
    #Derivate of Cdh1 used to calculate length of G1 phase
    dCdt = diff(result[:,3])/diff(times)

    '''
    #Plot of dCdh1/dt versus t    
    figure()
    plot(dCdt)
    xlabel('time')
    ylabel('d[Cdh1]/dt')
    '''
    
    #Randomly sample 'relative mass' of cell
    samplemass.append(result[np.random.randint(0,80),4])
    
    #Heuristic to determine G1 Length from Cdh1 derivative
    startS = 0
    endM =0
    count =0
    
    while startS == 0:
        if dCdt[count]<-0.10:
            startS = count
        count += 1
    count =0
    
    while endM ==0 and count < len(dCdt)-1:
        if dCdt[count]>0.40:
            endM = count
        count +=1
    
    if startS>0 and endM >0 and endM > startS:
        G1length.append(startS)
        massatS.append((result[startS,4]/0.45)*15*log(2))
        begmass.append((result[0,4]/0.45)*15*log(2))
        S_G2_Mlength.append(endM - startS)
        
    else:
        broke+=1
           
#Histogram of Length of G1               
figure()
hist(G1length,20)
xlabel('time(min)')
ylabel('count')
title('G1 Duration')
xlim(10,100)  

#Histogram of Length of S/G2/M combined
figure()
hist(S_G2_Mlength,20)
xlabel('time(min)')
ylabel('count')
title('S_G2_M Duration')
xlim(10,100)              

lengths=[]
#run fba model first
#length of S phase determined by variations in cell mass at beginning of cell cycle
for x in massatS:
    k = 2*1e3*1e12/6.022e23
    genome = 1.21e7         # length of S. Cerevcplex_interfaceisiae genome (bp)
    stoich = dnamodel.reactions.DNAsynth.metabolites[dnamodel.metabolites.DNA]  #stoichiometric coefficient of DNA (sum of dNTP) in DNAsynth reaction
    lengths.append(60*(k*genome)/(x*dnamodel1.reactions.DNAsynth.flux*stoich))


G2_Mlength=[]
for x in range(0, len(lengths)):
    G2_Mlength.append(S_G2_Mlength[x]-lengths[x])
    pleh = int(round(lengths[x]+G1length[x]))
    
    
    massatG2.append((result[pleh,4]/0.45)*15*log(2))



import seaborn
figure()
seaborn.distplot(G1length, hist=False, color = 'r', label = 'G1')
seaborn.distplot(G2_Mlength, hist = False, color = 'g', label = 'G2/M')
seaborn.distplot(lengths, hist = False, color = 'b', label =  'S')
legend(loc= 'best')
ylim([0,0.15])
xlabel('Time (min)')


s = []*len(samplemass)
for i in samplemass:
    s.append((i/0.45)*50*log(2))
    
    
figure()
seaborn.distplot(s, hist=False, color = 'b')
xlabel('Cell Size (pg)')
title('Cell Size Distribution for an Asynchronous Population')


#ratio of dna to protein is relatively constant

datp=[0]*(len(G1length))
arg=[0]*(len(G1length))
ratio=[0]*(len(G1length))


for l in range(0,len(G1length)):            
    datp[l] = 0.0036*(dnamodel2.reactions.DNAsynth.flux*G1length[l]*begmass[l]
                     +dnamodel1.reactions.DNAsynth.flux*lengths[l]*massatS[l]
                     +dnamodel2.reactions.DNAsynth.flux*G2_Mlength[l]*massatG2[l])
    
    arg[l] = 0.1607*(dnamodel2.reactions.r_2133.flux*G1length[l]*begmass[l]
                    +dnamodel1.reactions.r_2133.flux*lengths[l]*massatS[l]
                    +dnamodel2.reactions.r_2133.flux*G2_Mlength[l]*massatG2[l])
    
    
    ratio[l] = arg[l]/datp[l]

figure()
hist(ratio,20, color = 'green')        
plt.axvline(x=average(ratio),color='r',linewidth='2.0')            
plt.axvline(x=0.1607/0.0036,color='b',linewidth='2.0')            
xlabel('[protein]/[dNTP]')
ylabel('count')            
            
       
            
            
            
            
            
            
            
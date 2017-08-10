#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:08:16 2017

@author: rohanroy
"""

import roadrunner
from pylab import *

chen = roadrunner.RoadRunner("BIOMD0000000056.xml")
varIC = roadrunner.RoadRunner("BIOMD0000000056.xml")

chen.timeCourseSelections = ['time',
 'MASS',
 'CDC14',
 'BUD',
 'CDC20',
 'CDH1',
 'IEP',
 'CLB2',
 'CDC6',
 'SIC1',
 'CLB5',
 'NET1',
 'RENTP',
 'CDC6P',
 'SIC1P',
 'CDC20i',
 'PDS1',
 'SWI5',
 'NET1P',
 'F2P',
 'C5P',
 'C2P',
 'F5P',
 'C2',
 'F5',
 'CLN2',
 'CDC15',
 'PPX',
 'SPN',
 'TEM1GTP',
 'CDH1i',
 'ESP1',
 'SWI5P',
 'RENT',
 'C5',
 'F2',
 'ORI',]

#Simulate ODE model for 94 minutes (exactly one cell cycle)
result = chen.simulate(0,196,197)
times = result[:,0]



for i in range(1, 8):
    series=result[:,i]
    name = chen.timeCourseSelections[i]
    plot(times, series, label = str(name))
    legend(loc='best')


#Initial Condition is the last row of values of simulation
Initial_Conditions = result[196,:]

for i in range(1, len(Initial_Conditions)):
    varIC.model['init('+chen.timeCourseSelections[i]+')']=Initial_Conditions[i]



G1length=[]         #length of G1 Phase
S_G2_Mlength=[]     #length of combined S/G2/M phase
broke=0             #Broke simulations have uncontrolled growth without division
samplemass = []     #'Relative Mass' of cell randomly sampled at any point in the cell cycle
massatS = []        #'Relative Mass' of cell at the beginning of S phase
begmass = []        #'Relative Mass' of cell at beginning of cell cycle
massatG2 = []       #'Relative Mass' of cell at the beginning of G2 phase

sim=1000            #number of simulations
for i in range(0,sim):
    varIC.reset()
    for species in varIC.getIndependentFloatingSpeciesIds():
        # Each cell has a different size at birth due to assymetric division
        if species == 'MASS':
            division = np.random.normal(0.5,0.05)
            varIC.model['init(MASS)']= division*2.48
        #Variable Initial Conditions drawn from Normal distribution, mu = reported value, sigma = mu*0.01
        else:
            varIC.model['init('+species+')']=np.random.normal(chen.model['init('+species+')'],chen.model['init('+species+')']*0.30)                    
        
        varIC.timeCourseSelections = ['time',
                                 'MASS',
                                 'CDC14',
                                 'BUD',
                                 'CDC20',
                                 'CDH1',
                                 'CLB2',
                                 'CDC6',
                                 'SIC1',
                                 'ORI',
                                 'CLB5',
                                 'NET1',
                                 'RENTP',
                                 'CDC6P',
                                 'SIC1P',
                                 'CDC20i',
                                 'PDS1',
                                 'SWI5',
                                 'NET1P',
                                 'F2P',
                                 'C5P',
                                 'C2P',
                                 'F5P',
                                 'C2',
                                 'F5',
                                 'CLN2',
                                 'CDC15',
                                 'PPX',
                                 'SPN',
                                 'TEM1GTP',
                                 'CDH1i',
                                 'ESP1',
                                 'SWI5P',
                                 'RENT',
                                 'C5',
                                 'F2',
                                 'IEP',]
    result = varIC.simulate(0,250,251)
    times = result[:,0]
    
    '''
    #Plot of stochastic simulation
    figure()
    for j in range(1, 8):
        series=result[:,j]
        name = varIC.timeCourseSelections[j]
        plot(times, series, label = str(name))
        #legend(loc='best')
        xlabel('Time (min)')
        ylabel('concentration')
        title('Stochastic Cell Cycle')
       

    #ORI Plot
    figure()
    plot(times,result[:,9], label = 'ORI')
    legend(loc='best')
    xlabel('Time (min)')
    '''
    
    samplemass.append(result[np.random.randint(0,90),1])
    
    #Heuristic to determine G1 Length from ORI
    startS = 0
    endM =0
    count =0
    
    while startS == 0:
        if result[count,9]>1.0:
            startS = count
        count += 1

    
    while endM ==0:
        if result[count,9]<1.0:
            endM = count
        count +=1
    
    if startS>0 and endM >0 and endM > startS:
        G1length.append(startS)
        massatS.append((result[startS,1]/1.14)*15*log(2))
        begmass.append((result[0,1]/1.14)*15*log(2))
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
    
    
    massatG2.append((result[pleh,1]/1.14)*15*log(2))



import seaborn
figure()
seaborn.distplot(G1length, hist=False, color = 'r', label = 'G1')
seaborn.distplot(G2_Mlength, hist = False, color = 'g', label = 'G2/M')
seaborn.distplot(lengths, hist = False, color = 'b', label =  'S')
legend(loc= 'best')
xlabel('Time (min)')

   
s = []*len(samplemass)
for i in samplemass:
    s.append((i/1.14)*50*log(2))
    
    
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
            
              
    

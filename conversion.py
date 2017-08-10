#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:41:31 2017

@author: rohanroy
"""

import json

       
conversion = {     
'r_0958': 'PC',
'r_0984': 'RPE', 
'r_0466': 'G6PDH2r', 
'r_0887': 'PFK_3',
'r_1048': 'TALA', 
'r_0961': 'PDHm',
'r_0990': 'FBA3', 
'r_0886': 'PFK',
'r_0962': 'PYK', 
'r_0714': 'MDH',
'r_0568': 'PPA', 
'r_0569': 'PPAm', 
'r_0454': 'FRDm',
'r_0173': 'ALDD2y',
'r_0884': 'PPCK',
'r_0982': 'RPI', 
'r_0449': 'FBP',
'r_0534': 'HEX1', 
'r_0356': 'DPGM',
'r_0959': 'PYRDC', 
'r_1797': 'FRUK',
'r_0091': 'PGL',
'r_0366': 'ENO',
'r_0893': 'PGM',
'r_0889': 'GND', 
'r_0892': 'PGK', 
'r_0300': 'CSm', 
'r_0112': 'ACS',
'r_0113': 'ACSm', 
'r_0486': 'GAPD', 
'r_0467': 'PGI', 
'r_0450': 'FBA', 
'r_1054': 'TPI', 
'r_0568': 'PPA',
'r_0568': 'PPA',
'r_1166': 'GLCt1',
'r_1049': 'TKT1', 
'r_1050': 'TKT2',
'r_0322': 'FBA2',
'r_2034': 'PYRt2m',
'r_1254': 'PYRt2',
'r_1633': 'ACALDt',
'r_0163': 'ALCD2x_copy1',
'r_1762': 'ETOHt',
'r_1763': 'ETOHtm',
'r_0165': 'ALCD2irm',
'r_0174': 'ALDD2xm',
'r_0175': 'ALDD2ym',
'r_1106': 'ACt2r',
'r_1245': 'PIt2m',
'r_1128': 'CITtcm',
'r_1686': 'CITt2r',
'r_2131': 'ICDHym',
'r_0658': 'ICDHxm',
'r_0832': 'AKGDam',
'r_1696': 'CO2tm',
'r_1697': 'CO2t',
'r_0505': 'GCC2cm_copy1',
'r_0831': 'AKGDbm',
'r_1022': 'SUCOASm',
'r_1264': 'SUCCtm',
'r_1265': 'SUCFUMtm',
'r_1796': 'FUMt2r',
'r_1021': 'SUCD2_u6m',
'r_0454': 'FRDm',
'r_0451': 'FUMm',
'r_0713': 'MDHm',
'r_1226': 'MALtm',
'r_0718': 'ME1m',
'r_1901': 'MALt2r',
'r_0714': 'MDH',
'r_1239': 'OAAt2m',
}

for reaction in dnamodel1.reactions:
    reaction.biggid = None
    for rid in conversion:
        if reaction.id == rid:
            reaction.biggid = conversion[rid]
    

for reaction in dnamodel2.reactions:
    reaction.biggid = None
    for rid in conversion:
        if reaction.id == rid:
            reaction.biggid = conversion[rid]
        

dnamodel1_escher = {}
for reactions in dnamodel1.reactions:
    if reactions.biggid is not None:
        dnamodel1_escher[reactions.biggid]=reactions.averageflux

dnamodel2_escher = {}
for reactions in dnamodel2.reactions:
    if reactions.biggid is not None:
        dnamodel2_escher[reactions.biggid]=reactions.averageflux
        
comparison_escher = [dnamodel2_escher,dnamodel1_escher]

ratio = {}
for rxn in dnamodel1_escher:
   if dnamodel2_escher[rxn]==0 and dnamodel1_escher[rxn]==0:
        ratio[rxn] = 1
   elif dnamodel2_escher[rxn]==0:
       ratio[rxn] = 1000
   else:
       ratio[rxn] = dnamodel1_escher[rxn]/dnamodel2_escher[rxn]
   
        
import json

with open('dnamodel1_escher.json','w') as f:
    json.dump(dnamodel1_escher, f)

with open('dnamodel2_escher.json', 'w') as f:
    json.dump(dnamodel2_escher, f)

with open('comparison_escher.json', 'w') as f:
    json.dump(comparison_escher, f)
    
with open('ratio_escher.json', 'w') as f:
    json.dump(ratio, f)    
    

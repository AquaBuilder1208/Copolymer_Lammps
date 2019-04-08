# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 13:38:41 2018

@author: HP ElifeBook
"""

import numpy as np
import numpy.random as nr
#import pandas as pd
import Init 

# This python is used to generated random number, only for -CF3.

#box,atomlist=Init.Split_file('PVDF-model.cif')

def Count_C(atomlist):
    count=0
    for atom in atomlist.value:count=count+atom.flag
    return count

def Randomlize(atomlist):
    rate=nr.randint(38,42) #HFP rate
    total_C=Count_C(atomlist)
    total_HVP=int(total_C/2*rate/100+0.5)
    print('We Would Need to Randomlize '+str(total_HVP)+' indexs')
    other_index,choose_index=Return_index_Other(atomlist),[]
    while len(choose_index)<total_HVP:
        new_index=nr.choice(other_index,1)
        monoindex=atomlist.search_index(new_index).Monomer_print_index()
        for index in monoindex:
            if index in other_index:other_index.remove(index)  
        choose_index.append(new_index)
    return np.array(choose_index)       
                
def Return_index_Other(atomlist):
    index=[]
    for atom in atomlist.value:
        if atom.flag==0: index.append(atom.index)
    return index
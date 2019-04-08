# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:27:48 2018

@author: Gaussian
"""

import numpy as np
import numpy.random as nr
import pandas as pd
#from scipy.optimize import root,fsolve
import Init,Random,Rebuild

class Dummy(Init.Atom):
    def __init__(self,x,y,z,name,index):
        super(Dummy,self).__init__(x,y,z,name,index)
        self.appeared,self.flag=0,3

#box,atomlist=Init.Split_file('PVDF-model.cif')
#oldatomlist=atomlist
#print('A Good Thing is that We Successfully Initialize the File!')
#build_index=Random.Randomlize(atomlist)
#atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
#if not Good_random:
#    while Good_random!=True:
#        print('Starting Regeneration, this would be time-consuming')
#        box,atomlist=Init.Split_file('PVDF-model.cif')
#        #atomlist=oldatomlist
#        build_index=Random.Randomlize(atomlist)
#        atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
#print('Luckily We have a Good Random Number! Though Still '+ str(count)+' Warning was Made. Forget it.')

#max_index=Rebuild.Max_index(atomlist)

def Calc_mass(atomlist):
    mass=0
    for atom in atomlist.value:
        if atom.name[0]=='C' or atom.flag==1:mass=mass+12
        elif atom.name[0]=='F':mass=mass+18.998
        elif atom.name[0]=='H':mass=mass+1.008
        elif atom.name[0]=='O':mass=mass+15.9994 
    return mass

def Calc_dummy(atomlist):
    mass=Calc_mass(atomlist)
    num_dummy=int(mass*0.15/278)+5
    print('We Would Need to Randomlize '+str(num_dummy)+' Plasticizers')
    rate=num_dummy*243/mass
    print('That would be approx. '+str(int(rate*100))+'% in weight')
    return num_dummy

def Random_dummy(atomlist,box):
    max_index=Rebuild.Max_index(atomlist)
    other_index=Random.Return_index_Other(atomlist)
    num_dummy=Calc_dummy(atomlist)
    choose_index,choose_list=[],[]
    while len(choose_index)<num_dummy:
        new_index=nr.choice(other_index,1)
        monoindex=atomlist.search_index(new_index).Monomer_print_index()
        for index in monoindex:
            if index in other_index:other_index.remove(index)  
        choose_index.append(new_index)
        choose_list.append(atomlist.search_index(new_index))
    choose_list=Init.Atomlist(choose_list)
    dummylist=[]
    for atom in choose_list.value:
        bx,by,bz=box.x,box.y,box.z
        carbon=atom.neighbour[0]
        rx,ry,rz=atom.x*bx,atom.y*by,atom.z*bz
        Rx,Ry,Rz=carbon.x*bx,carbon.y*by,carbon.z*bz
        vect1=np.array([Rx-rx,Ry-ry,Rz-rz])
        vect1=vect1/np.linalg.norm(vect1)*5#Estimated Radii
        qx,qy,qz=(np.array([Rx,Ry,Rz])+vect1).tolist()
        dummy=Dummy(qx/bx,qy/by,qz/bz,'X'+str(max_index+1)+'x',max_index+1)
        max_index=max_index+1
        atomlist.value.append(dummy)
        dummylist.append(dummy)
    dummylist=Init.Atomlist(dummylist)
    print('Successful Creat Dummy Atoms')
    return atomlist,dummylist
    

        
        
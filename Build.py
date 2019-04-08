# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 20:06:54 2018

@author: Gaussian
"""

import numpy as np
#import numpy.random as nr
import pandas as pd
#from scipy.optimize import root,fsolve
import Init,Random,Rebuild,Dummy

# This python outputs Data & POSCAR, the later one is only used for test.

#def build_data(inputfile='PVDF-model.cif',filename='Data'):
def build_data(atomlist,box,filename='Data'):
    #box,atomlist=Init.Split_file(inputfile)
    #oldatomlist=atomlist
    #print('A Good Thing is that We Successfully Initialize the File!')
    #build_index=Random.Randomlize(atomlist)
    #atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
    #if not Good_random:
    #    while Good_random!=True:
    #        print('Starting Regeneration, this would be time-consuming')
    #        box,atomlist=Init.Split_file(inputfile)
    #        #atomlist=oldatomlist
    #        build_index=Random.Randomlize(atomlist)
    #        atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
    #print('Luckily We have a Good Random Number! Though Still '+ str(count)+' Warning was Made. Forget it.')
        #到这里为止 成功的加入了大量的CF3基团
    #atomlist,dummylist=Dummy.Random_dummy(atomlist,box)
    print('Build A File Of LAMMPS DATA for Main Part')            
    savingmsg = open(filename+'.txt', 'w')
    savingmsg.write('Auto Generated Data From Cif file\n\n')
    savingmsg.write(str(len(atomlist.value))+' '+'atoms\n')
#    savingmsg.write(str(len(atomlist.value)-len(dummylist.value))+' '+'atoms\n')
    savingmsg.write('4'+' '+'atom types\n')
    savingmsg.write('\n') 
    savingmsg.write('0 '+str(box.ux)+' xlo xhi\n') 
    savingmsg.write('0 '+str(box.uy)+' ylo yhi\n')
    savingmsg.write('0 '+str(box.uz)+' zlo zhi\n')
    savingmsg.write('\nMasses\n\n') 
    savingmsg.write('1 12.011\n') 
    savingmsg.write('2 18.998\n') 
    savingmsg.write('3 1.008\n')
    savingmsg.write('4 16\n')
    savingmsg.write('\nAtoms # charge\n\n')
    for atom in atomlist.value:
        if atom.name[0]=='C' or atom.flag==1: type1=1
        elif atom.name[0]=='F': type1=2
        elif atom.name[0]=='H': type1=3
        elif atom.name[0]=='O': type1=4
        savingmsg.write(str(atom.index)+' '+str(type1)+' '+str(0.0)+' '+str(atom.x*box.x)+' '+str(atom.y*box.y)+' '+str(atom.z*box.z)+'\n')
    savingmsg.close()
#    build_data_dummy(dummylist,box,filename)
    build_poscar(box,atomlist,filename)
    
def build_data_dummy(dummylist,box,filename):
    print('Build A File Of LAMMPS DATA For Dummy Atoms')            
    savingmsg = open(filename+'_dummy.txt', 'w')
    savingmsg.write('Auto Generated Data From Cif file\n\n')
    savingmsg.write(str(len(dummylist.value))+' '+'atoms\n')
    savingmsg.write('1'+' '+'atom types\n')
    savingmsg.write('\n') 
    savingmsg.write('0 '+str(box.ux)+' xlo xhi\n') 
    savingmsg.write('0 '+str(box.uy)+' ylo yhi\n')
    savingmsg.write('0 '+str(box.uz)+' zlo zhi\n')
    savingmsg.write('\nMasses\n\n') 
    savingmsg.write('1 243.00\n') 
    savingmsg.write('\nAtoms # charge\n\n')
    for atom in dummylist.value:
        #if atom.name[0]=='C' or atom.flag==1: type1=1
        #elif atom.name[0]=='F': type1=2
        #elif atom.name[0]=='H': type1=3
        type1=1 
        savingmsg.write(str(atom.index)+' '+str(type1)+' '+str(0.0)+' '+str(atom.x*box.x)+' '+str(atom.y*box.y)+' '+str(atom.z*box.z)+'\n')
    savingmsg.close()

def build_poscar(atomlist,box,filename='Data'):
    # 我知道这看起来很傻逼
    #box,atomlist=Init.Split_file('PVDF-model.cif')
    #build_index=Random.Randomlize(atomlist)
    #atomlist=Rebuild_H_to_CF3(atomlist,build_index,box)
    print('Build A File Of POSCAR')
    savingmsg = open(filename+'.vasp', 'w')
    savingmsg.write('Auto Generated Poscar From Cif file\n1.0\n')
    struct=np.zeros([3,3])
    struct[0,0],struct[1,1],struct[2,2]=box.ux,box.uy,box.uz 
    savingmsg.close()
    pd.DataFrame(struct).to_csv(filename+'.vasp',sep=' ',header=False,index=False,mode='a')
    savingmsg = open(filename+'.vasp', 'a')
    savingmsg.write('F H C O\n')
    count_C,count_F,count_H,count_O=0,0,0,0
    for atom in atomlist.value:
        if atom.name[0]=='C' or atom.flag==1: count_C=count_C+1
        elif atom.name[0]=='F': count_F=count_F+1
        elif atom.name[0]=='H': count_H=count_H+1
        elif atom.name[0]=='O': count_O=count_O+1
    savingmsg.write(str(count_F)+' '+str(count_H)+' '+str(count_C)+' '+str(count_O)+'\n')
    savingmsg.write('Direct\n')
    for atom in atomlist.value:
        if atom.name[0]=='F' and atom.flag==0: savingmsg.write(str(atom.x)+' '+str(atom.y)+' '+str(atom.z)+'\n')
    for atom in atomlist.value:
        if atom.name[0]=='H' and atom.flag==0: savingmsg.write(str(atom.x)+' '+str(atom.y)+' '+str(atom.z)+'\n')
    for atom in atomlist.value:
        if atom.name[0]=='C' or atom.flag==1: savingmsg.write(str(atom.x)+' '+str(atom.y)+' '+str(atom.z)+'\n')
    for atom in atomlist.value:
        if atom.name[0]=='O' and atom.flag==0: savingmsg.write(str(atom.x)+' '+str(atom.y)+' '+str(atom.z)+'\n')
    savingmsg.close()
    

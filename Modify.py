# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 19:25:56 2018

@author: HP ElifeBook
"""

import numpy as np
#import numpy.random as nr
import pandas as pd
#from scipy.optimize import root,fsolve
import Init,Random,Rebuild,Dummy,Build
import random

class Box_locate(Init.Box):
    def __init__(self,x,y,z):
        super(Box_locate,self).__init__(x,y,z)
        atomlist=[]
        self.atomlist=Init.Atomlist(atomlist)
    def locate(self,lx,ly,lz):self.lx,self.ly,self.lz=float(lx),float(ly),float(lz)
    def make_box(self,ux,uy,uz,mv):
        _,atomlist=Build_one()
        self.atomlist=atomlist
        for atom in atomlist.value:
            atom.x,atom.y,atom.z=(atom.x*self.x+self.lx)/ux,(atom.y*self.y+self.ly-mv[0])/uy,(atom.z*self.z+self.lz-mv[1])/uz
    def add_atom(self,atomlist,ux,uy,uz):
        dx=np.random.random(1)[0]
        self.dx=dx
        #atom=Init.Atom(0.5,0.5,0.5,'O',0)
        #atom.flag=0
        #atomlist.value=[atom]
        for atom in atomlist.value:
            atom.x=((atom.x)*15+self.lx)/ux
            atom.y=(atom.y*15+self.ly)/uy
            atom.z=(atom.z*15+self.lz)/uz
        #for atom in atomlist.value:
        #    atom.x,atom.y,atom.z=(atom.x*self.x+self.lx)/ux,(atom.y*self.y+self.ly)/uy,(atom.z*self.z+self.lz)/uz
        self.atomlist=atomlist
        self.flag=1
    def addmore_atom(self,atomlist,ux,uy,uz):
        #dx=np.random.random(1)[0]/3
        #self.dx=self.dx+dx
        for atom in atomlist.value:
            atom.x=((atom.x+self.dx)*15+15*self.flag+self.lx)/ux
            atom.y=(atom.y*15+self.ly)/uy
            atom.z=(atom.z*15+self.lz)/uz
        self.atomlist.value=self.atomlist.value+atomlist.value
        self.flag=self.flag+1
        
class Modbox(object):
    def __init__(self,x,y,z,num=1,dm=0):
         self.x,self.y,self.z=float(x),float(y),float(z)
         self.num=int(num)
         self.ux,self.uy,self.uz=x,num*y+dm*num,num*z+dm*num
         self.dm=dm
         self.otherboxlist=[]
    def divide_box(self,mv):
        boxlist=[]
        for j in range(self.num):
           for k in range(self.num):
               box_locate=Box_locate(self.x,self.y,self.z)
               box_locate.locate(0,j*self.y+self.dm*j,k*self.z+self.dm*k)#changed
               box_locate.make_box(self.ux,self.uy,self.uz,mv)
               boxlist.append(box_locate)
        print('\n######## Finish Initializing Supercell ########\n')
        emptyboxlist=[]
        ir=int(self.ux/(self.dm+0.1))
        for i in range(ir):
            for j in range(self.num):
                for k in range(self.num):
                    box_locate1=Box_locate(self.x,self.dm,self.z)
                    box_locate1.locate(i*self.dm,(j+1)*self.y+j*self.dm,k*self.z+k*self.dm)
                    box_locate2=Box_locate(self.x,self.y,self.dm)
                    box_locate2.locate(i*self.dm,j*self.y+j*self.dm,(k+1)*self.z+k*self.dm)
                    box_locate3=Box_locate(self.x,self.dm,self.dm)
                    box_locate3.locate(i*self.dm,(j+1)*self.y+(j+1)*self.dm,(k+1)*self.z+(1+k)*self.dm)
                    emptyboxlist=emptyboxlist+[box_locate1,box_locate2,box_locate3]
        self.boxlist=boxlist
        self.emptyboxlist=emptyboxlist
    def add_atomlist(self):
        atomlist_tot1=[]
        atomlist_tot2=[]
        for box in self.boxlist:
            atomlist_tot1=atomlist_tot1+box.atomlist.value
        for box in self.otherboxlist:
            atomlist_tot2=atomlist_tot2+box.atomlist.value
        self.atomlist1=Init.Atomlist(atomlist_tot1)
        self.atomlist1.reindex()
        self.atomlist2=Init.Atomlist(atomlist_tot2)
        self.atomlist2.reindex()        
            
    
class Superbox(Modbox):
    def __init__(self,x,y,z):
        super(Superbox,self).__init__(x,y,z)
        self.atomlist=Init.Atomlist([])
        self.flag=0
    def add_modbox(self,modbox):
        for atom in modbox.atomlist.value:
            atom.x=(atom.x*modbox.ux+self.flag*modbox.ux)/self.x
        self.atomlist.value=self.atomlist.value+modbox.atomlist.value
        self.flag=self.flag+1

def Build_one(inputfile='PVDF-model.cif'):
    box,atomlist=Init.Split_file(inputfile)
    #oldatomlist=atomlist
    #print('A Good Thing is that We Successfully Initialize the File!')
    build_index=Random.Randomlize(atomlist)
    atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
    if not Good_random:
        while Good_random!=True:
            print('Starting Regeneration, this would be time-consuming')
            box,atomlist=Init.Split_file(inputfile)
            #atomlist=oldatomlist
            build_index=Random.Randomlize(atomlist)
            atomlist,Good_random,count=Rebuild.Rebuild_H_to_CF3(atomlist,build_index,box)
    print('Luckily We have a Good Random Number! Though Still '+ str(count)+' Warning was Made. Forget it.')
    return box,atomlist

def Add_plastic(modbox,numplast):
    #print(len(modbox.emptyboxlist+modbox.emptyboxlist))
    boxchoice=random.sample(modbox.emptyboxlist,numplast)
    #boxchoice=modbox.emptyboxlist+modbox.emptyboxlist
    for box in boxchoice:
        index=np.random.randint(7)+1
        _,atomlist=Init.Split_file('DBP-'+str(index)+'.cif')
        if len(box.atomlist.value)==0:
            box.add_atom(atomlist,modbox.ux,modbox.uy,modbox.uz)
        else:
            #print('REUSE!')
            box.addmore_atom(atomlist,modbox.ux,modbox.uy,modbox.uz)
    modbox.otherboxlist=boxchoice

def Start(inputfile='PVDF-model.cif',size=2000,dm=15,mv=[0,0]):
    size=size/2
    box,atomlist=Init.Split_file(inputfile)
    num=int(pow(size/len(atomlist.value),1/2)+1)
    print('Now Generating a Supercell of ',1,'*',num,'*',num,'\n')
    print('######## Now Initializing Supercell ########\n')
    modbox=Modbox(box.x,box.y,box.z,num,dm)
    modbox.divide_box(mv)
    modbox.add_atomlist()
    #print('Though required size was ',size,', in total we make ',len(modbox.atomlist.value),' atoms\n')
    print('######## Now Generating Plasticizer ########\n')
    numplast=Dummy.Calc_dummy(modbox.atomlist1)
    Add_plastic(modbox,numplast)
    modbox.add_atomlist()
    print('Though required size was ',size*2,', in total we make ',len(modbox.atomlist1.value+modbox.atomlist2.value),' atoms, including ',len(modbox.atomlist2.value),' atoms from plasticizers\n')
    mass1=Dummy.Calc_mass(modbox.atomlist1)
    mass2=Dummy.Calc_mass(modbox.atomlist2)
    print('The actual Mass Ratio is ',mass2/(mass1+mass2)*100,'\n')
    Build.build_poscar(modbox.atomlist1,modbox,filename='Polymer')
    Build.build_poscar(modbox.atomlist2,modbox,filename='Plastic')
    Build.build_poscar(Init.Atomlist(modbox.atomlist1.value+modbox.atomlist2.value),modbox,filename='Combine')
    #modbox.uy,modbox.uz=modbox.uy+modbox.num*10,modbox.uz+modbox.num*10
    #return modbox

def get_key(value,modboxdict):return list (modboxdict.keys()) [list (modboxdict.values()).index (value)]

def Continue(num):
    modboxlist=[]
    for i in range(num):
        modbox=Start()
        modboxlist.append(modbox)
    numlist=np.arange(num).tolist()
    modboxdict=dict(zip(numlist,modboxlist))
    superbox=Superbox(modboxlist[0].ux*num,modboxlist[0].uy,modboxlist[0].uz)
    for modbox in modboxlist:
        superbox.add_modbox(modbox)
    Build.build_data(superbox.atomlist,superbox)
    
    
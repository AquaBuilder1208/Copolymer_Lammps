# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 20:41:21 2018

@author: Gaussian
"""

import numpy as np
import pandas as pd

class Box(object):
    def __init__(self,x,y,z):
        self.x,self.y,self.z=float(x),float(y),float(z)
        
        
# Class Box is the simulation box. Sorry for that it must be ortho at this time
#   But you may change it as you want. Just remember to redefine the build part

class Atom(object):
    def __init__(self,x,y,z,name,index):
        self.x,self.y,self.z,self.name=float(x),float(y),float(z),str(name)
        self.neighbour,self.monomer=[],[]
        self.index=int(index)
    def Reindex(self,index):self.index=index
    def Appear(self):self.appeared=1 # record that a carbon has appeared in a monomer.
    def Append_monomer(self,atom):self.monomer.append(atom)
    def Append_neighbour(self,atom):self.neighbour.append(atom)
    def Monomer_print(self):
        namelist=[]
        for monomer in self.monomer:namelist.append(monomer.name)
        return np.unique(namelist)
    def Monomer_print_index(self):
        namelist=[]
        for monomer in self.monomer:namelist.append(monomer.index)
        return np.unique(namelist)
    def Clean_monomer(self,atomlist):
        temp=[]
        for name in self.Monomer_print():temp.append(atomlist.search_name(str(name)))
        self.monomer=temp
 
# Class Atom is a typical atom defined. With keyword x,y,z,index,monomer,neighbour.      
# The keyword Flag in class Atom is use to determine whether it is a Carbon or not.
# The keyword Monomer records a typical MONOMER defined. C2H2F2.
# The keyword Neighbour records atoms it bonded to
        
class Carbon(Atom):
    def __init__(self,x,y,z,name,index):
        super(Carbon,self).__init__(x,y,z,name,index)
        self.appeared,self.flag=0,1
        
class Other(Atom):
    def __init__(self,x,y,z,name,index):
        super(Other,self).__init__(x,y,z,name,index)
        self.appeared,self.flag=0,0

# Class Atomlist is a list of atoms, with methods to search for certain atoms using ID or name.
#   A subclass named Chain is used to record atoms of a full chain. 
        
class Atomlist(object):
    def __init__(self,value):self.value=value
    def search_index(self,index):
        for atom in self.value:
            if atom.index==index:return atom   
    def search_name(self,name):
        for atom in self.value: 
            if atom.name==name:return atom  
    def reindex(self):
        for i in range(len(self.value)):
            self.value[i].index=i+1
        
def Strip_line(lines):
    new_lines=[]
    for subline in lines:
        subline=subline.split('\n')
        subline=subline[0].split()
        new_lines.append(subline)
    return new_lines

def Split_file(file):
    with open (file,'r') as f:
        line=f.readlines()
        info_lines=pd.DataFrame(Strip_line(line[9:15]))
        box=Read_infolines(info_lines)
        atomlist=Loop_for_coord(line[24:])
        return box,atomlist

def Read_infolines(info_lines):
    xyz,angle=np.array(info_lines[1]).astype('float')[0:3],np.array(info_lines[1]).astype('float')[3:6]
    #if sum( angle-90 == 0 ) != 3:
    #    print('Error, not a Ortho box, System Abort')
    #    return -1
    return Box(xyz[0],xyz[1],xyz[2])

def Loop_for_coord(line):
    coord_lines,bond_lines,flag,atomlist=[],[],0,[]
    for subline in line:
        if ('loop_' in subline):flag = 1
        if flag == 0:coord_lines.append(subline.split())
        else:bond_lines.append(subline.split())
    coord,bond=numer_coord(pd.DataFrame(coord_lines).iloc[:,0:5]),np.array(pd.DataFrame(bond_lines[6:]).iloc[:,0:3])
    for i in range(int(coord.size/6)):
        print('Now Processing ',i+1,' of ',int(coord.size/6),' atoms',end="\r")
        if str(coord.iloc[i,0][0])=='C':
            atomlist.append(Carbon(coord.iloc[i,2],coord.iloc[i,3],coord.iloc[i,4],coord.iloc[i,0],coord.iloc[i,5]))
        else:
            atomlist.append(Other(coord.iloc[i,2],coord.iloc[i,3],coord.iloc[i,4],coord.iloc[i,0],coord.iloc[i,5]))
    atomlist=Atomlist(atomlist)
    print('###################################')
    return purify_monomer(changing_bond(atomlist,bond))
           
def numer_coord(coord):
    #Note the index is given at the beginning and shoud NEVER changed.
    index=range(int(coord.size/5))
    new_coord=pd.concat([coord,pd.DataFrame(np.array(index)+1)],axis=1)
    new_coord.columns=(['name','type','x','y','z','index'])
    return new_coord
                    
def changing_bond(atomlist,bond): 
    print('Bond info detected, totally',len(bond),'bonds')
    count=0
    for bond_info in bond:
        print('Now processing', count+1,' of ', len(bond),' bond infos',end="\r")
        count=count+1
        atom1_name,atom2_name=bond_info[0],bond_info[1]
        atom1,atom2=atomlist.search_name(atom1_name),atomlist.search_name(atom2_name)               
        if (atom1.flag == 1 and atom2.flag == 1):
            #此时两个都是Carbon类
            if (atom1.appeared == 0 and atom2.appeared == 0):
                #（此时两者均未出现，未被计入其他单体），不可忽略本信息
                atom1.Append_monomer(atom2)
                atom2.Append_monomer(atom1)
                atom1.Appear()
                atom2.Appear()
            else:
                # 此时至少其中一者已出现过，信息无意义
                atom1.Append_neighbour(atom2)
                atom2.Append_neighbour(atom1)
                continue
        if (atom1.flag == 0 and atom2.flag == 0):1
            #此时信息属于有误,两个异原子不应该相互成键。
            #print('Wrong Bonding Information Detected')
            #raise ValueError
        else:
            #此时两者有一方是碳，只需要互相加入Monomer即可
            atom1.Append_monomer(atom2)
            atom2.Append_monomer(atom1)
        atom1.Append_neighbour(atom2)
        atom2.Append_neighbour(atom1)
    print('###################################')
    return atomlist
            
def purify_monomer(atomlist):  
    #This Function is used Twice at one time, this is because the low efficiency of my code. Sorry. 
    print('Cleaning Procedure 1 in running:')
    count=0
    for atom in atomlist.value:
        count=count+1
        print('Now Cleaning ', count+1 ,' of ',len(atomlist.value),' atoms',end="\r")
        for name in atom.Monomer_print():
            new_atom=atomlist.search_name(str(name))
            atom.monomer=atom.monomer+new_atom.monomer
        atom.Clean_monomer(atomlist)
    print('###################################')
    print('Cleaning Procedure 2 in running:')
    count=0
    for atom in atomlist.value:
        count=count+1
        print('Now Cleaning ', count+1 ,' of ',len(atomlist.value),' atoms',end="\r")
        for name in atom.Monomer_print():
            new_atom=atomlist.search_name(str(name))
            atom.monomer=atom.monomer+new_atom.monomer
        atom.Clean_monomer(atomlist)
    print('###################################')
    return atomlist

def pairinglist(atomlist):
        #if sum(self.poscar1.natom)!=sum(self.poscar2.natom): raise ArithmeticError
        #for atomtype in self.poscar1.atomdict:
        #    if atomtype not in self.poscar2.atomdict: raise AttributeError
        #    if self.poscar1.atomdict[atomtype]!=self.poscar2.atomdict[atomtype]: raise ArithmeticError
        #print('Input file Examined.')
        monomerlist,totlist,otherlist=[],[],[]
        print('Exanming Procedure in running:')
        count=0
        for atom in atomlist.value:
            print('Now Examining ', count+1 ,' of ',len(atomlist.value),' atoms',end="\r")
            count=count+1
            newlist,flag=[],0
            for inname in totlist:
                if atom.name==inname: flag=1
            if flag==0:
                #atom=self.cif1.search_name(name)
                monomer=atom.monomer
                for matom in monomer:newlist.append(matom.name)
                for matom in monomer:
                    if matom.name in totlist:
                        for mono in monomerlist:
                            if matom.name in mono:mono=np.unique(mono+newlist)
                    else:                                
                        newlist=np.unique(newlist).tolist()
                        monomerlist.append(newlist)
                        totlist=totlist+newlist
        for atom in atomlist.value:
            name=atom.name
            if name not in totlist: otherlist.append(name)
        print('###################################')
        return monomerlist,totlist,otherlist
                
            
            
        
        
        
    

    
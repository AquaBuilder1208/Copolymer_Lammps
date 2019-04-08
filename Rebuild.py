# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 15:11:34 2018

@author: Gaussian
"""

import numpy as np
import numpy.random as nr
#import pandas as pd
#from scipy.optimize import root,fsolve
import Init

# In this python, you will see how stupid I'm. This is used to rebuild CF3 group.
# Only god knows what I wrote.
# A vertical vector is used to rebuild. That's fxxking easy.

def Delta_x(x1,x2):
    if abs(x1-x2)<=abs(abs(x1-x2)-1):return x1-x2
    elif abs(x1-1-x2)<=0.5:return x1-1-x2
    else:return x1+1-x2
    
def Calc_position_C(atom,carbon,box,rebond=1.53):
    cx,cy,cz=carbon.x,carbon.y,carbon.z
    ax,ay,az=atom.x,atom.y,atom.z
    x,y,z=box.x,box.y,box.z
    bond=np.sqrt((Delta_x(cx,ax)*x)**2+(Delta_x(cy,ay)*y)**2+(Delta_x(cz,az)*z)**2)
    return [cx+(ax-cx)/bond*rebond,cy+(ay-cy)/bond*rebond,cz+(az-cz)/bond*rebond]
    
def Max_index(atomlist):  
    indexlist=[]
    for atom in atomlist.value:indexlist.append(atom.index)
    return(np.max(indexlist))
    
def Find_S(atom,carbon,box):
    bx,by,bz=box.x,box.y,box.z
    rx,ry,rz=atom.x*bx,atom.y*by,atom.z*bz
    Rx,Ry,Rz=carbon.x*bx,carbon.y*by,carbon.z*bz
    # Compute first F with random
    vect1=np.array([Rx-rx,Ry-ry,Rz-rz])/1.73
    # e1=np.array([1,0,0]) generate too much problems.
    e1=np.array([nr.rand(1)[0],nr.rand(1)[0],nr.rand(1)[0]])
    e1=e1/np.linalg.norm(e1)
    vect2=Cross(vect1,e1)
    ibl=1.2 # imagninal bond length
    ix,iy,iz=Calc_position_C(atom,carbon,box,rebond=1.73)
    ilocation=np.array([ix*bx,iy*by,iz*bz])
    xf1,yf1,zf1=ilocation+ibl*vect2
    vect3=Cross(vect2,vect1)
    xf2,yf2,zf2=ilocation-1/2*ibl*vect2+np.sqrt(3)/2*ibl*vect3
    xf3,yf3,zf3=ilocation-1/2*ibl*vect2-np.sqrt(3)/2*ibl*vect3
    return [xf1,yf1,zf1,xf2,yf2,zf2,xf3,yf3,zf3]

def Cross(vect,e1):
    # A cross for 1D vector. FXXK NUMPY.
    vect2=np.array([vect[1]*e1[2]-vect[2]*e1[1],vect[2]*e1[0]-vect[0]*e1[2],vect[0]*e1[1]-vect[1]*e1[0]])
    return vect2/np.linalg.norm(vect2)

def Distance(atom,sub_f1,sub_f2,sub_f3,box):
    bx,by,bz=box.x,box.y,box.z
    xa,ya,za=atom.x,atom.y,atom.z
    sublist=[sub_f1,sub_f2,sub_f3]
    flag=1
    for sub in sublist:
        xf,yf,zf=sub.x,sub.y,sub.z
        dist=np.sqrt((Delta_x(xf,xa)*bx)**2+(Delta_x(yf,ya)*by)**2+(Delta_x(zf,za)*bz)**2)
        if dist <= 1.2:
            flag=0
    return flag
        

def Rebuild_H_to_CF3(atomlist,build_index,box):
    max_index=Max_index(atomlist)
    count=0
    Good_random=True
    for index in build_index:
        atom=atomlist.search_index(index)
        carbon=atom.neighbour[0]
        nx,ny,nz=Calc_position_C(atom,carbon,box)
        atom.flag,atom.x,atom.y,atom.z=1,nx,ny,nz
        atom.name='C'+atom.name[1:]+'s' #Yet it is still a Other Class
        xf1,yf1,zf1,xf2,yf2,zf2,xf3,yf3,zf3=Find_S(atom,carbon,box)
        i=0
        while i<20:
            flag=1
            xf1,yf1,zf1,xf2,yf2,zf2,xf3,yf3,zf3=Find_S(atom,carbon,box)
            sub_f1=Init.Other(xf1/box.x,yf1/box.y,zf1/box.z,'F'+str(max_index+1)+'s',max_index+1)
            sub_f1.neighbour.append(atom)
            sub_f2=Init.Other(xf2/box.x,yf2/box.y,zf2/box.z,'F'+str(max_index+2)+'s',max_index+2)
            sub_f2.neighbour.append(atom)
            sub_f3=Init.Other(xf3/box.x,yf3/box.y,zf3/box.z,'F'+str(max_index+3)+'s',max_index+3)
            sub_f3.neighbour.append(atom)
            Fslist=[]
            for atomo in atomlist.value:
                if atomo.name=='Fs': Fslist.append(atomo)
            for atomo in Fslist:
                if (Distance(atomo,sub_f1,sub_f2,sub_f3,box)==0):#当前距离太近
                    #print('Flag was set as zero')
                    flag=0
                    xf1,yf1,zf1,xf2,yf2,zf2,xf3,yf3,zf3=Find_S(atom,carbon,box)
                    sub_f1=Init.Other(xf1/box.x,yf1/box.y,zf1/box.z,'Fs',max_index+1)
                    sub_f2=Init.Other(xf2/box.x,yf2/box.y,zf2/box.z,'Fs',max_index+2)
                    sub_f3=Init.Other(xf3/box.x,yf3/box.y,zf3/box.z,'Fs',max_index+3)
            if flag==1:
                break
            i=i+1
        if i==20:
            print('################################################################################')
            print('#                                                                              #')
            print('#WARNING: A Too Near -CF3 Group was Detected and We can not Fix it by Rotating.#')
            print('#                                                                              #')
            print('################################################################################')
            count=count+1
            #build_index=Random.Randomlize(atomlist)
        #elif i>0:
            #print('A Warning about too near -CF3 Group was fixed by Changing Rotational Matrix')
        if count/len(build_index) >= 0.05:
            print('################################################################################')
            print('#                                                                              #')
            print('#        Too many unsolved Warning Detected, Regenerate a Initial Guess        #')
            print('#                                                                              #')
            print('################################################################################')
            Good_random=False
            break
        atomlist.value.append(sub_f1)
        atomlist.value.append(sub_f2)
        atomlist.value.append(sub_f3)
        max_index=max_index+3
    return atomlist,Good_random,count


    

    
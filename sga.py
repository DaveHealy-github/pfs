#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 09:44:10 2021

functions to replicate MATLAB of Allemendinger et al. 
Structural Geology Algorithms 

@author: davidhealy
"""
 
import numpy as np 

def deg2rad(angle):
    return angle * np.pi / 180. 

def rad2deg(angle):
    return angle * 180. / np.pi 

#   force angles to range 0, 2*pi 
def ZeroTwoPi(angle):
    #   returns angle in radians in range 0, 2*pi  
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    angle2 = angle 
    
    if angle2 < 0.:
        angle2 += 2. * np.pi 
    
    elif angle2 >= (2. * np.pi):
        angle2 -= 2. * np.pi 
        
    return angle2  

#   convert coordinates from Cartesian to spherical
def CartToSph(cn, ce, cd):
    #   converts direction cosines in NED convention to [trd,plg] tuple
    #   trd & plg in radians 
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    
    #   plunge 
    plg = np.arcsin(cd)
    
    #   trend 
    if cn == 0.:
        if ce < 0.:
            trd = 3. / 2. * np.pi   #   trend is west 
        else:
            trd = np.pi / 2.
            
    else:
        trd = np.arctan(ce/cn)
        if cn < 0.:
            trd += np.pi 
        trd = ZeroTwoPi(trd)    
    
    return trd, plg 

#   convert coordinates from spherical to Cartesian
def SphToCart(trd, plg, k):
    #   converts [trd,plg] in radians to direction cosines [cn, ce, cd]
    #   in NED convention 
    #   k is an integer flag:
    #       k = 0: the trend and plunge of a line 
    #       k = 1: the strike and dip of a plane in right hand rule 
    #   In this last case, the direction cosines of the pole to the plane 
    #   are returned
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    
    #   for a line 
    if k == 0:
        cd = np.sin(plg)
        ce = np.cos(plg) * np.sin(trd)
        cn = np.cos(plg) * np.cos(trd)
        
    #   for a plane 
    if k == 1:
        cd = np.cos(plg)
        ce = -np.sin(plg) * np.cos(trd)
        cn = np.sin(plg) * np.sin(trd)
        
    return cn, ce, cd 

#   direction cosines of right-handed Cartesian system from trend and plunge 
#   of X1 and trend of X3 
def DirCosAxes(tX1, pX1, tX3):
    #   Returns 3x3 matrix dC with direction cosines of X1 (row1), X2 (row2)
    #   and X3 (row3)
    #   All inputs in radians 
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    
    #   some constants 
    east = np.pi / 2. 
    west = 1.5 * np.pi
    #   tolerance for near zero values 
    tol = 1.e-6 
    
    dC = np.zeros([3,3])
    
    #   dir cos of X1
    dC[0,0], dC[0,1], dC[0,2] = SphToCart(tX1, pX1, 0)
    
    #   plunge of X3
    if np.abs(pX1) < tol:
        dt = np.abs(tX1-tX3) 
        if np.abs(dt - east) < tol or np.abs(dt - west) < tol:
            pX3 = 0. 
        else:
            pX3 = east 
            
    else:
        pX3 = np.arctan(-(dC[0,0]*np.cos(tX3) + dC[0,1]*np.sin(tX3))/dC[0,2])
        
    #   dir cos of X3 
    dC[2,0], dC[2,1], dC[2,2] = SphToCart(tX3, pX3, 0)
    
    #   dir cos of X2 as cross product of X3 & X1
    dC[1,0] = dC[2,1]*dC[0,2] - dC[2,2]*dC[0,1]
    dC[1,1] = dC[2,2]*dC[0,0] - dC[2,0]*dC[0,2]
    dC[1,2] = dC[2,0]*dC[0,1] - dC[2,1]*dC[0,0]
    
    #   convert X2 to unit vector 
    r = np.sqrt(dC[1,0]*dC[1,0] + dC[1,1]*dC[1,1] + dC[1,2]*dC[1,2])
    for i in range(0,3):
        dC[1,i] = dC[1,i]/r 
        
    return dC 

#   calculate principal stresses from any given stress tensor 
def PrincipalStress(stress, tX1, pX1, tX3):
    #   Given a stress tensor, Cauchy() returns tractions on given plane 
    #   All input angles in radians 
    #
    #   Returns:
    #   ptress: 3x3 matrix of mag, trend and plunge of max, int, and min stresses 
    #   dCp: 3x3 matrix of dir cos of princip stresses 
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    
    #   dir cos of X1, X2, X3 
    dC = DirCosAxes(tX1, pX1, tX3)
    
    pstress = np.zeros([3,3])
    
    #   get eigenvalues and eigenvectors of stress tensor 
    eVal, eVec = np.linalg.eig(stress)
    idx = eVal.argsort()[::-1]
    eVal = eVal[idx]
    eVec = eVec[:,idx]
    pstress[0,0] = eVal[2]  #   max prin stress
    pstress[1,0] = eVal[1]
    pstress[2,0] = eVal[0]  #   min prin stress 
    
    #   convert dir cos of principal stresses to NED coord system 
    tV = np.zeros([3,3])
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                tV[j,i] = dC[k,j] * eVec[k,i] + tV[j,i]
    
    dCp = np.zeros([3,3])
    
    #   trd & plg of max prin stress
    dCp[0,:] = [tV[0,2], tV[1,2], tV[2,2]]
    pstress[0,1], pstress[0,2] = CartToSph(tV[0,2], tV[1,2], tV[2,2])
    
    #   trd & plg of int prin stress
    dCp[1,:] = [tV[0,1], tV[1,1], tV[2,1]]
    pstress[1,1], pstress[1,2] = CartToSph(tV[0,1], tV[1,1], tV[2,1])
    
    #   trd & plg of min prin stress
    dCp[2,:] = [tV[0,0], tV[1,0], tV[2,0]]
    pstress[2,1], pstress[2,2] = CartToSph(tV[0,0], tV[1,0], tV[2,0])
    
    return pstress, dCp 

def ShearOnPlane(stress, tX1, pX1, tX3, strike, dip):
    #   Given a stress tensor, ShearOnPlane() returns tractions on given plane 
    #   All input angles in radians 
    #
    #   Returns:
    #   tractions: 3x3 matrix of mag, trend and plunge of normal & shear stresses 
    #
    #   MATLAB script written by Nestor Cardozo for the book Structural 
    #   Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. 
    #   Converted to Python by Dave Healy, 2021 
    
    tractions = np.zeros([3,3])
    dCTT = np.zeros([3,3])
    
    #   principal stresses and directions 
    pstress, dCp = PrincipalStress(stress, tX1, pX1, tX3)
    for i in range(0,3):
        stress[i,i] = pstress[i,0]
    
    #   dir cos of pole to plane 
    p = np.zeros([3,])
    p[0], p[1], p[2] = SphToCart(strike, dip, 1)      
    
    #   transform pole to plane to principal stress coordinates
    pT = np.zeros([3,])
    for i in range(0,3):
        for j in range(0,3):
            pT[i] = dCp[i,j] * p[j] + pT[i]
            
    #   calc tractions in principal stress coordinates
    T = np.zeros([3,])
    #   tractions using Cauchy's law
    for i in range(0, 3):
        for j in range(0, 3):
            T[i] = stress[i,j] * pT[j] + T[i]
            
    #   find B axis 
    B = np.zeros(3)
    B[0] = T[1]*pT[2] - T[2]*pT[1]
    B[1] = T[2]*pT[0] - T[0]*pT[2]
    B[2] = T[0]*pT[1] - T[1]*pT[0]
    
    #   find shear direction 
    S = np.zeros(3)
    S[0] = pT[1]*B[2] - pT[2]*B[1]
    S[1] = pT[2]*B[0] - pT[0]*B[2]
    S[2] = pT[0]*B[1] - pT[1]*B[0]
    
    #   convert T, B and S to unit vectors 
    rB = np.sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2])
    rS = np.sqrt(S[0]*S[0] + S[1]*S[1] + S[2]*S[2])
    for i in range(0,3):
        B[i] = B[i]/rB 
        S[i] = S[i]/rS 
        
    #   transformation matrix from prin stress to plane coords 
    a = np.zeros([3,3])
    a[0,:] = [pT[0], pT[1], pT[2]]
    a[1,:] = [B[0], B[1], B[2]]
    a[2,:] = [S[0], S[1], S[2]]
    
    #   stress ratio, R
    #R = (stress[1,1] - stress[0,0])/(stress[2,2] - stress[0,0])
    
    #   normal and shear tractions 
    for i in range(0,3):
        tractions[i,0] = (stress[0,0]*a[0,0]*a[i,0] + 
                          stress[1,1]*a[0,1]*a[i,1] +
                          stress[2,2]*a[0,2]*a[i,2])
    
    #   tractions in NED coord system 
    for i in range(0,3):
        for j in range(0,3):
            dCTT[0,i] = dCp[j,i]*pT[j] + dCTT[0,i]
            dCTT[1,i] = dCp[j,i]*B[j] + dCTT[1,i]
            dCTT[2,i] = dCp[j,i]*S[j] + dCTT[2,i]
            
    #   trend and plunge of traction on plane 
    tractions[0,1], tractions[0,2] = CartToSph(dCTT[0,0],dCTT[0,1],dCTT[0,2])
    tractions[1,1], tractions[1,2] = CartToSph(dCTT[1,0],dCTT[1,1],dCTT[1,2])
    tractions[2,1], tractions[2,2] = CartToSph(dCTT[2,0],dCTT[2,1],dCTT[2,2])

    return tractions 
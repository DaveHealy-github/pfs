#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 11:07:49 2021

pfs.py
Probability of Fault Slip = pfs 

Response Surface Methodology applied to fault & fracture stability 
- slip tendency, dilation tendency, fracture susceptibility 

Import library of functions 

@author: davidhealy
"""

#   standard stuff 
import numpy as np 
import matplotlib.pyplot as plt
import scipy.stats as stats
import sga

degree_sign = '\u00b0'

def deg2rad(angle):
    return angle * np.pi / 180. 

def rad2deg(angle):
    return angle * 180. / np.pi 

def plotMohr(s1, s2, s3, muRef, C0Ref):
    
    meanS1 = np.mean(s1)
    meanS2 = np.mean(s2)
    meanS3 = np.mean(s3)
    
    theta = np.arange(0.,360.,1.)   
    sin2theta = np.sin(2. * theta * np.pi/180.) 
    cos2theta = np.cos(2. * theta * np.pi/180.) 
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    for i in range(0, len(s1)-1):
        
        tau13Mohr = ((s1[i] - s3[i])/2.) * sin2theta  
        sigma13Mohr = (s1[i] + s3[i])/2. + ((s1[i] - s3[i])/2.) * cos2theta
        
        tau12Mohr = ((s1[i] - s2[i])/2.) * sin2theta
        sigma12Mohr = (s1[i] + s2[i])/2. + ((s1[i] - s2[i])/2.) * cos2theta  
        
        tau23Mohr = ((s2[i] - s3[i])/2.) * sin2theta  
        sigma23Mohr = (s2[i] + s3[i])/2. + ((s2[i] - s3[i])/2.) * cos2theta  
        
        ax.plot(sigma13Mohr, tau13Mohr, 'gainsboro', alpha=0.5)
        ax.plot(sigma12Mohr, tau12Mohr, 'gainsboro', alpha=0.5)
        ax.plot(sigma23Mohr, tau23Mohr, 'gainsboro', alpha=0.5)
    
    tau13Mean = ((meanS1 - meanS3)/2.) * sin2theta  
    sigma13Mean = (meanS1 + meanS3)/2. + ((meanS1 - meanS3)/2.) * cos2theta
    
    tau12Mean = ((meanS1 - meanS2)/2.) * sin2theta
    sigma12Mean = (meanS1 + meanS2)/2. + ((meanS1 - meanS2)/2.) * cos2theta  
    
    tau23Mean = ((meanS2 - meanS3)/2.) * sin2theta  
    sigma23Mean = (meanS2 + meanS3)/2. + ((meanS2 - meanS3)/2.) * cos2theta 
        
    ax.plot(sigma13Mean, tau13Mean, '-r', label=r'Mean $\sigma_{13}$')
    ax.plot(sigma12Mean, tau12Mean, '-g', label=r'Mean $\sigma_{12}$')
    ax.plot(sigma23Mean, tau23Mean, '-b', label=r'Mean $\sigma_{23}$')
    
    sMu = r'$\mu=$%1.2f' % muRef    
    ax.plot([0, s1.max()], [C0Ref, muRef*s1.max()+C0Ref], '--k', label=sMu, alpha=0.5)
    
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Mohr diagram for stress variations, n=%i' % len(s1))
    ax.set_xlim(0, np.max(s1)*1.1)
    ax.set_ylim(0, ((np.max(s1)-np.min(s3))/2.)*1.4) 
    ax.set_xlabel('Normal stress, MPa')
    ax.set_ylabel('Shear stress, MPa')
    ax.legend(loc='upper right')
    
    plt.savefig('pfsInputMohr.png', dpi=300)

    return 

def plotMohr2(s1, s2, s3, sn, tau):
    
    meanS1 = np.mean(s1)
    meanS2 = np.mean(s2)
    meanS3 = np.mean(s3)
    
    theta = np.arange(0.,360.,1.)   
    sin2theta = np.sin(2. * theta * np.pi/180.) 
    cos2theta = np.cos(2. * theta * np.pi/180.) 
    
    fig, ax = plt.subplots(figsize=(6,4))
    
    for i in range(0, len(s1)):
        
        tau13Mohr = ((s1[i] - s3[i])/2.) * sin2theta  
        sigma13Mohr = (s1[i] + s3[i])/2. + ((s1[i] - s3[i])/2.) * cos2theta
        
        tau12Mohr = ((s1[i] - s2[i])/2.) * sin2theta
        sigma12Mohr = (s1[i] + s2[i])/2. + ((s1[i] - s2[i])/2.) * cos2theta  
        
        tau23Mohr = ((s2[i] - s3[i])/2.) * sin2theta  
        sigma23Mohr = (s2[i] + s3[i])/2. + ((s2[i] - s3[i])/2.) * cos2theta  
        
        ax.plot(sigma13Mohr, tau13Mohr, 'gainsboro')
        ax.plot(sigma12Mohr, tau12Mohr, 'gainsboro')
        ax.plot(sigma23Mohr, tau23Mohr, 'gainsboro')
    
        ax.plot(sn, tau, 'ok')
        
    tau13Mean = ((meanS1 - meanS3)/2.) * sin2theta  
    sigma13Mean = (meanS1 + meanS3)/2. + ((meanS1 - meanS3)/2.) * cos2theta
    
    tau12Mean = ((meanS1 - meanS2)/2.) * sin2theta
    sigma12Mean = (meanS1 + meanS2)/2. + ((meanS1 - meanS2)/2.) * cos2theta  
    
    tau23Mean = ((meanS2 - meanS3)/2.) * sin2theta  
    sigma23Mean = (meanS2 + meanS3)/2. + ((meanS2 - meanS3)/2.) * cos2theta 
        
    ax.plot(sigma13Mean, tau13Mean, '-r', label='sigma1-3')
    ax.plot(sigma12Mean, tau12Mean, '-g', label='sigma1-2')
    ax.plot(sigma23Mean, tau23Mean, '-b', label='sigma2-3')
        
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Mohr diagram for stress variation, n=%i' % len(s1))
#    ax.set_xlim(0, np.max(s1)*1.05)
    ax.set_ylim(0, ((np.max(s1)-np.min(s3))/2.)*1.2) 
    ax.set_xlabel('Normal stress, MPa')
    ax.set_ylabel('Shear stress, MPa')
    
    plt.savefig('pfsInputMohr2.png', dpi=300)

    return 

def plotMohr3(s1, s2, s3, mu, modeMu, C0, modeC0):
    
    meanS1 = np.mean(s1)
    meanS2 = np.mean(s2)
    meanS3 = np.mean(s3)

    theta = np.arange(0.,360.,1.)   
    sin2theta = np.sin(2. * theta * np.pi/180.) 
    cos2theta = np.cos(2. * theta * np.pi/180.) 
    
    fig, ax = plt.subplots(figsize=(6,5))
    
    for i in range(0, len(s1)-1):
        
        tau13Mohr = ((s1[i] - s3[i])/2.) * sin2theta  
        sigma13Mohr = (s1[i] + s3[i])/2. + ((s1[i] - s3[i])/2.) * cos2theta
        
        tau12Mohr = ((s1[i] - s2[i])/2.) * sin2theta
        sigma12Mohr = (s1[i] + s2[i])/2. + ((s1[i] - s2[i])/2.) * cos2theta  
        
        tau23Mohr = ((s2[i] - s3[i])/2.) * sin2theta  
        sigma23Mohr = (s2[i] + s3[i])/2. + ((s2[i] - s3[i])/2.) * cos2theta  
        
        ax.plot(sigma13Mohr, tau13Mohr, '-r', alpha=0.1, lw=0.5)
        ax.plot(sigma12Mohr, tau12Mohr, '-g', alpha=0.1, lw=0.5)
        ax.plot(sigma23Mohr, tau23Mohr, '-b', alpha=0.1, lw=0.5)
    
        ax.plot([0, s1.max()*1.1], [C0[i], mu[i]*s1[i]*1.1+C0[i]], 'gainsboro', alpha=0.5)
        
    tau13Mean = ((meanS1 - meanS3)/2.) * sin2theta  
    sigma13Mean = (meanS1 + meanS3)/2. + ((meanS1 - meanS3)/2.) * cos2theta
    
    tau12Mean = ((meanS1 - meanS2)/2.) * sin2theta
    sigma12Mean = (meanS1 + meanS2)/2. + ((meanS1 - meanS2)/2.) * cos2theta  
    
    tau23Mean = ((meanS2 - meanS3)/2.) * sin2theta  
    sigma23Mean = (meanS2 + meanS3)/2. + ((meanS2 - meanS3)/2.) * cos2theta 
        
    ax.plot(sigma13Mean, tau13Mean, '-r', lw=4, label=r'$\sigma_{1}-\sigma_{3}$')
    ax.plot(sigma12Mean, tau12Mean, '-g', lw=4, label=r'$\sigma_{1}-\sigma_{2}$')
    ax.plot(sigma23Mean, tau23Mean, '-b', lw=4, label=r'$\sigma_{2}-\sigma_{3}$')
    
    sMu = r'$\mu=$%1.2f, $C_0$=%2.1f' % (modeMu, modeC0)    
    ax.plot([0, s1.max()*1.1], [modeC0, modeMu*s1.max()*1.1+modeC0], '-k', label=sMu, alpha=0.5)
    
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Mohr diagram for input variations, n=%i' % len(s1))
    ax.set_xlim(0, np.max(s1)*1.1)
    ax.set_ylim(0, ((np.max(s1)-np.min(s3))/2.)*1.7) 
    ax.set_xlabel('Normal stress, MPa')
    ax.set_ylabel('Shear stress, MPa')
    ax.legend(loc='upper left')
    
    plt.savefig('pfsInputMohr3.png', dpi=300)

    return 

def plotMohr4(s1, s2, s3, Pf, mu, modeMu, C0, modeC0):
    
    meanS1 = s1.mean()
    meanS2 = s2.mean()
    meanS3 = s3.mean()
    meanPf = Pf.mean()
    
    theta = np.arange(0.,360.,1.)   
    sin2theta = np.sin(2. * theta * np.pi/180.) 
    cos2theta = np.cos(2. * theta * np.pi/180.) 
    
    fig, ax = plt.subplots(figsize=(6,5))
    
    for i in range(0, len(s1)-1):
        
        tau13Mohr = ((s1[i] - s3[i])/2.) * sin2theta  
        sigma13Mohr = (s1[i] + s3[i])/2. + ((s1[i] - s3[i])/2.) * cos2theta
        
        tau12Mohr = ((s1[i] - s2[i])/2.) * sin2theta
        sigma12Mohr = (s1[i] + s2[i])/2. + ((s1[i] - s2[i])/2.) * cos2theta  
        
        tau23Mohr = ((s2[i] - s3[i])/2.) * sin2theta  
        sigma23Mohr = (s2[i] + s3[i])/2. + ((s2[i] - s3[i])/2.) * cos2theta  
        
        ax.plot(sigma13Mohr, tau13Mohr, 'gainsboro', alpha=0.1, lw=0.5)
        ax.plot(sigma12Mohr, tau12Mohr, 'gainsboro', alpha=0.1, lw=0.5)
        ax.plot(sigma23Mohr, tau23Mohr, 'gainsboro', alpha=0.1, lw=0.5)
    
        tau13Mohr = ((s1[i] - s3[i])/2.) * sin2theta  
        sigma13Mohr = (s1[i] + s3[i])/2. + ((s1[i] - s3[i])/2.) * cos2theta - Pf[i]
        
        tau12Mohr = ((s1[i] - s2[i])/2.) * sin2theta
        sigma12Mohr = (s1[i] + s2[i])/2. + ((s1[i] - s2[i])/2.) * cos2theta - Pf[i]  
        
        tau23Mohr = ((s2[i] - s3[i])/2.) * sin2theta  
        sigma23Mohr = (s2[i] + s3[i])/2. + ((s2[i] - s3[i])/2.) * cos2theta - Pf[i]  
        
        ax.plot(sigma13Mohr, tau13Mohr, '-r', alpha=0.1, lw=0.5)
        ax.plot(sigma12Mohr, tau12Mohr, '-g', alpha=0.1, lw=0.5)
        ax.plot(sigma23Mohr, tau23Mohr, '-b', alpha=0.1, lw=0.5)
        
        ax.plot([0, s1.max()*1.1], [C0[i], mu[i]*s1[i]*1.1+C0[i]], 'gainsboro', alpha=0.5)
        
    tau13Mean = ((meanS1 - meanS3)/2.) * sin2theta  
    sigma13Mean = (meanS1 + meanS3)/2. + ((meanS1 - meanS3)/2.) * cos2theta - meanPf
    
    tau12Mean = ((meanS1 - meanS2)/2.) * sin2theta
    sigma12Mean = (meanS1 + meanS2)/2. + ((meanS1 - meanS2)/2.) * cos2theta - meanPf 
    
    tau23Mean = ((meanS2 - meanS3)/2.) * sin2theta  
    sigma23Mean = (meanS2 + meanS3)/2. + ((meanS2 - meanS3)/2.) * cos2theta - meanPf
        
    ax.plot(sigma13Mean, tau13Mean, '-r', lw=4, label=r'$\sigma_{1}-\sigma_{3}$')
    ax.plot(sigma12Mean, tau12Mean, '-g', lw=4, label=r'$\sigma_{1}-\sigma_{2}$')
    ax.plot(sigma23Mean, tau23Mean, '-b', lw=4, label=r'$\sigma_{2}-\sigma_{3}$')
    
    sMu = r'$\mu=$%1.2f, $C_0$=%2.1f' % (modeMu, modeC0)    
    ax.plot([0, s1.max()*1.1], [modeC0, modeMu*s1.max()*1.1+modeC0], '-k', label=sMu, alpha=0.5)
    
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Mohr diagram for input variations, n=%i' % len(s1))
    ax.set_xlim(0, np.max(s1)*1.1)
    ax.set_ylim(0, ((np.max(s1)-np.min(s3))/2.)*1.7) 
    ax.set_xlabel('Effective normal stress, MPa')
    ax.set_ylabel('Shear stress, MPa')
    ax.legend(loc='upper left')
    
    plt.savefig('pfsInputMohr4.png', dpi=300)

    return 

#   return tau and sigmaN on specified plane 
def calcStressOnPlane(s1, s2, s3, sHaz, strike, dip):
    #   calculates shear and normal stress components on a plane 

    #   need to revise this to make it generic... 
    stress = np.zeros([3,3])
    stress[0,0] = s1 
    stress[1,1] = s3 
    stress[2,2] = s2 

    tX1 = deg2rad(sHaz) 
    pX1 = deg2rad(0.) 
    tX3 = deg2rad(sHaz)
    strike = deg2rad(strike)  
    dip = deg2rad(dip)
    
    tractions = sga.ShearOnPlane(stress, tX1, pX1, tX3, strike, dip)
    
    tau = tractions[2,0] 
    sigmaN = tractions[0,0] 
    
    return sigmaN, tau 

#   return tau and sigmaN on specified plane 
def calcAndersonianStressOnPlane(sV, sH, sh, sHaz, strike, dip):
    #   calculates shear and normal stress components on a plane 

    stress = np.zeros([3,3])
    stress[0,0] = sH        #   default x direction 
    stress[1,1] = sh        #   default y direction 
    stress[2,2] = sV        #   default z direction 

    tX1 = deg2rad(sHaz) 
    pX1 = deg2rad(0.) 
    tX3 = deg2rad(sHaz)
    strike = deg2rad(strike)  
    dip = deg2rad(dip)
    
    tractions = sga.ShearOnPlane(stress, tX1, pX1, tX3, strike, dip)
    
    tau = tractions[2,0]
    sigmaN = tractions[0,0] 
    
    return sigmaN, tau 

#   normal probability plot of residuals 
def plotResiduals(x, sLabel1, sLabel2):
    
    #   Calculate quantiles and least-square-fit curve
    (quantiles, values), (slope, intercept, r) = stats.probplot(x.reshape(-1), dist='norm')

    #   define ticks
    ticks_perc=[1, 5, 10, 20, 50, 80, 90, 95, 99]
    #   transfrom them from precentile to cumulative density
    ticks_quan=[stats.norm.ppf(i/100.) for i in ticks_perc]

    #   xlimits 
    xlim = np.max([np.abs(np.min(values)), np.max(values)])*1.05

    #   plot results
    fig, ax = plt.subplots(figsize=(4,4))
    
    if 's' in sLabel1:
        ax.plot(values, quantiles,'or')
    elif 'd' in sLabel1:
        ax.plot(values, quantiles,'oC0')
    else:
        ax.plot(values, quantiles,'og')
    ax.plot(quantiles * slope + intercept, quantiles, '--k')
    ax.set_xlim(-xlim, xlim)
    ax.set_yticks(ticks_quan)
    ax.set_yticklabels(ticks_perc)
    ax.set_xlabel('Residual ' + sLabel1)
    ax.set_ylabel('Probability, %')
    ax.set_box_aspect(1)
    ax.grid(True)
    ax.set_title('Normal probability plot - ' + sLabel1)
    plt.text(np.min(x), quantiles[-2], sLabel2, 
             bbox=dict(facecolor='white', edgecolor='none'))
    plt.tight_layout() 
    fn = 'pfs' + sLabel1 + '.png'
    plt.savefig(fn, dpi=300)

    return 

#   normal probability plot of residuals 
def plotResiduals2(x, sLabel1, sLabel2):
    
    #   Calculate quantiles and least-square-fit curve
    (quantiles, values), (slope, intercept, r) = stats.probplot(x.reshape(-1), dist='norm')

    #   define ticks
    ticks_perc=[1, 5, 10, 20, 50, 80, 90, 95, 99]
    #   transfrom them from precentile to cumulative density
    ticks_quan=[stats.norm.ppf(i/100.) for i in ticks_perc]

    #   xlimits 
    xlim = np.max([np.abs(np.min(values)), np.max(values)])*1.05
    xlim = 0.19
    
    #   plot results
    fig, ax = plt.subplots(figsize=(4,4))
    
    if 's' in sLabel1:
        ax.plot(values, quantiles,'or')
    elif 'd' in sLabel1:
        ax.plot(values, quantiles,'oC0')
    else:
        ax.plot(values, quantiles,'og')
    ax.plot(quantiles * slope + intercept, quantiles, '--k')
    ax.set_xlim(-xlim, xlim)
    ax.set_yticks(ticks_quan)
    ax.set_yticklabels(ticks_perc)
    ax.set_xlabel('Residual ' + sLabel1)
    ax.set_ylabel('Probability, %')
    ax.set_box_aspect(1)
    ax.grid(True)
    ax.set_title('Normal probability plot - ' + sLabel1)
    plt.text(-xlim*0.9, quantiles[-2], sLabel2, 
             bbox=dict(facecolor='white', edgecolor='none'))
    plt.tight_layout() 
    fn = 'pfs' + sLabel1 + '.png'
    plt.savefig(fn, dpi=300)

    return 

#   ANOVA & residuals for regression solution 
def calcANOVAresid(y, X, beta):
    
    nd, nq = X.shape
    ny = len(y)
    
    #   sums of squares 
    SST = np.dot(y.T, y) - pow(np.sum(y),2.) / ny ; 
    SSR = np.dot(np.dot(beta.T, X.T), y) - pow(np.sum(y),2.) / ny ;
    SSE = SST - SSR ; 
    
    q = nq - 1 ; 
    MSR = SSR / q ; 
    MSE = SSE / (ny - q - 1) ; 
    Fstat = MSR / MSE ; 
    
    Rsquared = SSR / SST ; 
    RsquaredAdj = 1 - (SSE / (ny - q - 1)) / (SST / (ny - 1)) ; 

    return Rsquared, RsquaredAdj, Fstat, MSR, MSE

def calctStats(X, beta, q, MSE):

    #   W = (X'X)^-1
    W = np.linalg.inv(np.dot(X.T,X)) 
    
    rX, cX = X.shape 
    N = rX 
    df = N - q - 1 
    alpha = 0.05 
    tcrit = stats.t.ppf(1 - alpha/2, df)
    
    print('t-statistics, t_critical (alpha = %1.2f) = %3.4f' % (alpha, tcrit))
    t = np.zeros([q,])
    for j in range(1,q+1):
        t[j-1] = beta[j] / np.sqrt(MSE * W[j,j])
        if np.abs(t[j-1]) > tcrit:  
            print('Variable %i: Significant, |t|=%3.4f' % (j, np.abs(t[j-1])))
        else:
            print('Variable %i: *** NOT significant, |t|=%3.4f' % (j, np.abs(t[j-1])))
        
    return 

#   main effects for specified variables 
def plotMainEffectsTsTd(yTs, yTd, X, statX):
    
    dummy1 = statX[:,0].T
    S1Ts = (np.mean(yTs[np.ix_(X[:,1]==-1)]), 
            np.mean(yTs[np.ix_(X[:,1]==0)]),
            np.mean(yTs[np.ix_(X[:,1]==1)]))
    S1Td = (np.mean(yTd[np.ix_(X[:,1]==-1)]), 
            np.mean(yTd[np.ix_(X[:,1]==0)]),
            np.mean(yTd[np.ix_(X[:,1]==1)]))

    dummy2 = statX[:,1].T
    S2Ts = (np.mean(yTs[np.ix_(X[:,2]==-1)]), 
            np.mean(yTs[np.ix_(X[:,2]==0)]),
            np.mean(yTs[np.ix_(X[:,2]==1)]))
    S2Td = (np.mean(yTd[np.ix_(X[:,2]==-1)]), 
            np.mean(yTd[np.ix_(X[:,2]==0)]),
            np.mean(yTd[np.ix_(X[:,2]==1)]))

    dummy3 = statX[:,2].T
    S3Ts = (np.mean(yTs[np.ix_(X[:,3]==-1)]), 
            np.mean(yTs[np.ix_(X[:,3]==0)]),
            np.mean(yTs[np.ix_(X[:,3]==1)]))
    S3Td = (np.mean(yTd[np.ix_(X[:,3]==-1)]), 
            np.mean(yTd[np.ix_(X[:,3]==0)]),
            np.mean(yTd[np.ix_(X[:,3]==1)]))

    dummy4 = statX[:,3].T
    SHazTs = (np.mean(yTs[np.ix_(X[:,4]==-1)]), 
              np.mean(yTs[np.ix_(X[:,4]==0)]),
              np.mean(yTs[np.ix_(X[:,4]==1)]))
    SHazTd = (np.mean(yTd[np.ix_(X[:,4]==-1)]), 
              np.mean(yTd[np.ix_(X[:,4]==0)]),
              np.mean(yTd[np.ix_(X[:,4]==1)]))

    dummy5 = statX[:,4].T
    strikeTs = (np.mean(yTs[np.ix_(X[:,5]==-1)]), 
                np.mean(yTs[np.ix_(X[:,5]==0)]),
                np.mean(yTs[np.ix_(X[:,5]==1)]))
    strikeTd = (np.mean(yTd[np.ix_(X[:,5]==-1)]), 
                np.mean(yTd[np.ix_(X[:,5]==0)]),
                np.mean(yTd[np.ix_(X[:,5]==1)]))

    dummy6 = statX[:,5].T
    dipTs = (np.mean(yTs[np.ix_(X[:,6]==-1)]), 
             np.mean(yTs[np.ix_(X[:,6]==0)]),
             np.mean(yTs[np.ix_(X[:,6]==1)]))
    dipTd = (np.mean(yTd[np.ix_(X[:,6]==-1)]), 
             np.mean(yTd[np.ix_(X[:,6]==0)]),
             np.mean(yTd[np.ix_(X[:,6]==1)]))

    minPlotTs = np.min([S1Ts, S2Ts, S3Ts, SHazTs, strikeTs, dipTs])*.9  
    maxPlotTs = np.max([S1Ts, S2Ts, S3Ts, SHazTs, strikeTs, dipTs])*1.1  
    minPlotTd = np.min([S1Td, S2Td, S3Td, SHazTd, strikeTd, dipTd])*.9  
    maxPlotTd = np.max([S1Td, S2Td, S3Td, SHazTd, strikeTd, dipTd])*1.1  
    
    #   one graph per variable for 'main effects' 
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,8))

    axs[0,0].plot(dummy1, S1Ts, 'rs-', label='Slip tendency')
    axs[0,0].set_xlabel(r'$\sigma_1$, MPa')
    axs[0,0].set_ylabel('Tendency')
    axs[0,0].set_ylim(minPlotTs, maxPlotTs)
    axs[0,0].grid(True)
    axs[0,0].legend()
    axs[0,0].set_title(r'Main effects: $\sigma_1$')
    
    axs[0,1].plot(dummy2, S2Ts, 'rs-', label='Slip tendency')
    axs[0,1].set_xlabel(r'$\sigma_2$, MPa')
    axs[0,1].set_ylabel('Tendency')
    axs[0,1].set_ylim(minPlotTs, maxPlotTs)
    axs[0,1].grid(True)
    axs[0,1].legend()
    axs[0,1].set_title(r'Main effects: $\sigma_2$')
    
    axs[0,2].plot(dummy3, S3Ts, 'rs-', label='Slip tendency')
    axs[0,2].set_xlabel(r'$\sigma_3$, MPa')
    axs[0,2].set_ylabel('Tendency')
    axs[0,2].set_ylim(minPlotTs, maxPlotTs)
    axs[0,2].grid(True)
    axs[0,2].legend()
    axs[0,2].set_title(r'Main effects: $\sigma_3$')
    
    axs[1,0].plot(dummy4, SHazTs, 'rs-', label='Slip tendency')
    axs[1,0].set_xlabel(r'$\sigma_H$ azimuth, degrees')
    axs[1,0].set_ylabel('Tendency')
    axs[1,0].set_ylim(minPlotTs, maxPlotTs)
    axs[1,0].grid(True)
    axs[1,0].legend()
    axs[1,0].set_title(r'Main effects: $\sigma_H$ azimuth')
    
    axs[1,1].plot(dummy5, strikeTs, 'rs-', label='Slip tendency')
    axs[1,1].set_xlabel('Fault strike, degrees')
    axs[1,1].set_ylabel('Tendency')
    axs[1,1].set_ylim(minPlotTs, maxPlotTs)
    axs[1,1].grid(True)
    axs[1,1].legend()
    axs[1,1].set_title('Main effects: fault strike')
    
    axs[1,2].plot(dummy6, dipTs, 'rs-', label='Slip tendency')
    axs[1,2].set_xlabel('Fault dip, degrees')
    axs[1,2].set_ylabel('Tendency')
    axs[1,2].set_ylim(minPlotTs, maxPlotTs)
    axs[1,2].grid(True)
    axs[1,2].legend()
    axs[1,2].set_title('Main effects: fault dip')
    
    fig.tight_layout() 
    plt.savefig('pfsMainEffects_Ts.png', dpi=300)

    #   one graph per variable for 'main effects' 
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(12,8))

    axs[0,0].plot(dummy1, S1Td, 's-', c='C0', label='Dilation tendency')
    axs[0,0].set_xlabel(r'$\sigma_1$, MPa')
    axs[0,0].set_ylabel('Tendency')
    axs[0,0].set_ylim(minPlotTd, maxPlotTd)
    axs[0,0].grid(True)
    axs[0,0].legend()
    axs[0,0].set_title(r'Main effects: $\sigma_1$')
    
    axs[0,1].plot(dummy2, S2Td, 's-', c='C0', label='Dilation tendency')
    axs[0,1].set_xlabel(r'$\sigma_2$, MPa')
    axs[0,1].set_ylabel('Tendency')
    axs[0,1].set_ylim(minPlotTd, maxPlotTd)
    axs[0,1].grid(True)
    axs[0,1].legend()
    axs[0,1].set_title(r'Main effects: $\sigma_2$')
    
    axs[0,2].plot(dummy3, S3Td, 's-', c='C0', label='Dilation tendency')
    axs[0,2].set_xlabel(r'$\sigma_3$, MPa')
    axs[0,2].set_ylabel('Tendency')
    axs[0,2].set_ylim(minPlotTd, maxPlotTd)
    axs[0,2].grid(True)
    axs[0,2].legend()
    axs[0,2].set_title(r'Main effects: $\sigma_3$')
    
    axs[1,0].plot(dummy4, SHazTd, 's-', c='C0', label='Dilation tendency')
    axs[1,0].set_xlabel(r'$\sigma_H$ azimuth, degrees')
    axs[1,0].set_ylabel('Tendency')
    axs[1,0].set_ylim(minPlotTd, maxPlotTd)
    axs[1,0].grid(True)
    axs[1,0].legend()
    axs[1,0].set_title(r'Main effects: $\sigma_H$ azimuth')
    
    axs[1,1].plot(dummy5, strikeTd, 's-', c='C0', label='Dilation tendency')
    axs[1,1].set_xlabel('Fault strike, degrees')
    axs[1,1].set_ylabel('Tendency')
    axs[1,1].set_ylim(minPlotTd, maxPlotTd)
    axs[1,1].grid(True)
    axs[1,1].legend()
    axs[1,1].set_title('Main effects: fault strike')
    
    axs[1,2].plot(dummy6, dipTd, 's-', c='C0', label='Dilation tendency')
    axs[1,2].set_xlabel('Fault dip, degrees')
    axs[1,2].set_ylabel('Tendency')
    axs[1,2].set_ylim(minPlotTd, maxPlotTd)
    axs[1,2].grid(True)
    axs[1,2].legend()
    axs[1,2].set_title('Main effects: fault dip')
    
    fig.tight_layout() 
    plt.savefig('pfsMainEffects_Td.png', dpi=300)

    return 

def plotMainEffectsSf(ySf, X, statX): 
    
    dummy7 = statX[:,0].T
    muSf = (np.mean(ySf[np.ix_(X[:,8]==-1)]), 
            np.mean(ySf[np.ix_(X[:,8]==0)]),
            np.mean(ySf[np.ix_(X[:,8]==1)]))

    dummy8 = statX[:,1].T
    C0Sf = (np.mean(ySf[np.ix_(X[:,9]==-1)]), 
            np.mean(ySf[np.ix_(X[:,9]==0)]),
            np.mean(ySf[np.ix_(X[:,9]==1)]))

    dummy9 = statX[:,2].T
    PfSf = (np.mean(ySf[np.ix_(X[:,4]==-1)]), 
            np.mean(ySf[np.ix_(X[:,4]==0)]),
            np.mean(ySf[np.ix_(X[:,4]==1)]))

    minPlotSf = np.min([muSf, C0Sf, PfSf])*.9  
    maxPlotSf = np.max([muSf, C0Sf, PfSf])*1.1  
    
    #   one graph per variable for 'main effects' 
    fig, axs = plt.subplots(figsize=(12,4))

    ax7 = plt.subplot(131)
    ax7.plot(dummy7, muSf, 'gs-', label='Fracture susceptibility')
    ax7.set_xlabel(r'$\mu$')
    ax7.set_ylabel('Susceptibility')
    ax7.set_ylim(minPlotSf, maxPlotSf)
    ax7.grid(True)
    ax7.legend()
    ax7.set_title(r'Main effects: $\mu$')
    
    ax8 = plt.subplot(132)
    ax8.plot(dummy8, C0Sf, 'gs-', label='Fracture susceptibility')
    ax8.set_xlabel(r'$C_0$, MPa')
    ax8.set_ylabel('Susceptibility')
    ax8.set_ylim(minPlotSf, maxPlotSf)
    ax8.grid(True)
    ax8.legend()
    ax8.set_title(r'Main effects: $C_0$')
    
    ax9 = plt.subplot(133)
    ax9.plot(dummy9, PfSf, 'gs-', label='Fracture susceptibility')
    ax9.set_xlabel(r'$P_f$, MPa')
    ax9.set_ylabel('Susceptibility')
    ax9.set_ylim(minPlotSf, maxPlotSf)
    ax9.grid(True)
    ax9.legend()
    ax9.set_title(r'Main effects: $P_f$')
    
    fig.tight_layout() 
    plt.savefig('pfsMainEffects_Sf.png', dpi=300)
    
    return 

def plotStressPolygon(sV, sH, sh):
    
    minV = sV.min()
    maxV = sV.max()
    meanV = sV.mean()
    maxH = sH.max()
    minh = sh.min()
    
    fig, ax = plt.subplots(figsize=(4,4))
    
    ax.plot([0, maxV],[0, maxV], '--k')
    ax.plot([minh, minh],[minh, maxH], '-r')
    ax.plot([minh, maxH],[maxH, maxH], '-b')
    ax.plot(meanV, meanV, 'oy')
    ax.plot([minV, maxV],[minV, maxV],'-y')
    ax.scatter(sh, sH, alpha=0.5)

    ax.set_xlim(0, maxV)
    ax.set_ylim(0, maxV)
    ax.set_xlabel(r'$\sigma_h$, MPa')
    ax.set_ylabel(r'$\sigma_H$, MPa')
    ax.set_aspect('equal')
    ax.grid(True)
    ax.set_title('Stress polygon')
    
    fig.tight_layout() 
    plt.savefig('pfsStressPolygon.png', dpi=300)
    
    return 

def getStats(x):
    #   get some basic statistical parameters/moments of a distribution 
    minx = np.min(x)
    maxx = np.max(x)
    meanx = np.mean(x)
    medianx = np.median(x)
    
    return minx, maxx, meanx, medianx 

def getCircStats(x):
    #   get some basic statistical parameters/moments of a circular distribution 
    minx = np.min(x)
    maxx = np.max(x)
    meanx = x.mean()
    medianx = np.median(x)
    
    return minx, maxx, meanx, medianx 

def getCircStats2(x):
    #   get some basic statistical parameters/moments of a circular distribution 
    minx = np.min(x)
    maxx = np.max(x)
    meanx = x.mean()
    medianx = np.median(x)
    stdx = np.std(x)
    
    return minx, maxx, meanx, medianx, stdx 

def getFrictionDist(N):
    #   function to generate a distribution of N friction coefficients 
    #   assuming skew normal distribution, controlled by shape parameter alpha
    bInput = True     
    while bInput:
        
        muMean = input('Enter approximate mean value for friction [0.65]: ')
        if len(muMean) == 0:
            muMean = 0.65 
        else:
            muMean = float(muMean) 
           
        muStdDev = input('Enter approximate standard deviation for friction [10% of mean]: ')
        if len(muStdDev) == 0:
            muStdDev = muMean * .1   
        else:
            muStdDev = float(muStdDev) 
           
        muAlpha = input('Enter shape value (float) for skewness (left<0, symmetric=0, right>0) [-2.5]: ')        
        if len(muAlpha) == 0:
            muAlpha = -2.5
        else:
            muAlpha = float(muAlpha) 
         
        #   generate a skew normal distribution     
        mu = stats.skewnorm.rvs(muAlpha, muMean, muStdDev, N)
    
        fig, ax = plt.subplots(figsize=(6,6))
        nmu, bmu, pmu = ax.hist(mu, 20)
        
        sortedbmu = bmu[nmu.argsort()]
        incrbmu = bmu[1] - bmu[0]
        muMode = sortedbmu[-1] + incrbmu/2.
        muMin, muMax, muMean, muMedian = getStats(mu) 
        print('\n*** Summary statistics for this distribution')
        print('Minimum: %1.3f' % muMin)
        print('Maximum: %1.3f' % muMax)
        print('Median:  %1.3f' % muMedian)
        print('Mode:    %1.3f' % muMode)
        print('Mean:    %1.3f' % muMean)
        
        ax.plot([muMode, muMode], [0, max(nmu)*1.1], '-r')
        ax.set_xlabel('Friction coefficient')
        ax.set_ylabel('Count')
        ax.set_xlim(0, max(mu)*1.1)
        ax.grid(True)
        ax.set_title('Friction coefficient distribution, n=%i' % N)
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            bInput = False 

    return mu

def getCohesionDist(N):
    #   function to generate a distribution of N cohesion values 
    #   assuming skew normal distribution, controlled by shape parameter alpha
    bInput = True     
    while bInput:
        
        c0Mean = input('Enter approximate mean value for cohesion (MPa) [30.]: ')
        if len(c0Mean) == 0:
            c0Mean = 30. 
        else:
            c0Mean = float(c0Mean) 
           
        c0StdDev = input('Enter approximate standard deviation for cohesion [10% of mean]: ')
        if len(c0StdDev) == 0:
            c0StdDev = c0Mean * .1   
        else:
            c0StdDev = float(c0StdDev) 
           
        c0Alpha = input('Enter shape value (float) for skewness (left<0, symmetric=0, right>0) [-2.5]: ')        
        if len(c0Alpha) == 0:
            c0Alpha = -2.5
        else:
            c0Alpha = float(c0Alpha) 
 
        #   generate a skew normal distribution     
        c0 = stats.skewnorm.rvs(c0Alpha, c0Mean, c0StdDev, N)
    
        fig, ax = plt.subplots(figsize=(6,6))
        nc0, bc0, pc0 = ax.hist(c0, 20)

        sortedbc0 = bc0[nc0.argsort()]
        incrbc0 = bc0[1] - bc0[0]
        c0Mode = sortedbc0[-1] + incrbc0/2.
        c0Min, c0Max, c0Mean, c0Median = getStats(c0) 
        print('\n*** Summary statistics for this distribution')
        print('Minimum: %3.2f MPa' % c0Min)
        print('Maximum: %3.2f MPa' % c0Max)
        print('Median:  %3.2f MPa' % c0Median)
        print('Mode:    %3.2f MPa' % c0Mode)
        print('Mean:    %3.2f MPa' % c0Mean)
        
        ax.plot([c0Mode, c0Mode], [0, max(nc0)*1.1], '-r')
        ax.set_xlabel('Cohesion, MPa')
        ax.set_ylabel('Count')
        ax.set_xlim(0, max(c0)*1.1)
        ax.grid(True)
        ax.set_title('Cohesion distribution, n=%i' % N)
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            bInput = False 

    return c0

def getStrikeDist(N):
    #   function to generate a distribution of N strike values 
    #   assuming von Mises (circular normal) distribution, 
    #   controlled by shape parameter kappa
    bInput = True     
    while bInput:
        
        strikeMean = input('Enter mean value for strike (in degrees): ')
        if len(strikeMean) == 0:
            continue  
        else:
            strikeMean = float(strikeMean) 
           
        strikeKappa = input('Enter shape value (kappa) for von Mises distribution (low=dispersed, high=clustered) [5.0]: ')        
        if len(strikeKappa) == 0:
            strikeKappa = 5.0
        else:
            strikeKappa = float(strikeKappa) 
 
        #   von Mises circular normal distribution 
        strike = stats.vonmises.rvs(strikeKappa, strikeMean*np.pi/180, 1/strikeKappa, N)
        strike = strike * 180/np.pi 
        strikeMin, strikeMax, strikeMean, strikeMedian = getCircStats(strike) 
        strike[np.ix_(strike<0)] += 180.
        strike[np.ix_(strike>180)] -= 180.
        
        fig, ax = plt.subplots(figsize=(6,6))
        nstrike, bstrike, pstrike = ax.hist(strike, 20)
        
        sortedbstrike = bstrike[nstrike.argsort()]
        incrbstrike = bstrike[1] - bstrike[0]
        strikeMode = sortedbstrike[-1] + incrbstrike/2.
        print('\n*** Summary statistics for this distribution')
        print('Minimum: %3.2f%s' % (strikeMin, degree_sign))
        print('Maximum: %3.2f%s' % (strikeMax, degree_sign))
        print('Median:  %3.2f%s' % (strikeMedian, degree_sign))
        print('Mode:    %3.2f%s' % (strikeMode, degree_sign))
        print('Mean:    %3.2f%s' % (strikeMean, degree_sign))
        
        ax.plot([strikeMean, strikeMean], [0, max(nstrike)*1.1], '-r')
        ax.set_xlabel('Strike, degrees')
        ax.set_ylabel('Count')
        ax.set_xlim(0., 180.)
        ax.grid(True)
        ax.set_title('Strike distribution, n=%i' % N)
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            bInput = False 
    
    return strike 

def getDipDisp(N):
    #   function to generate a distribution of N dip values 
    #   assuming von Mises (circular normal) distribution, 
    #   controlled by shape parameter kappa
    bInput = True     
    while bInput:
        
        dipMean = input('Enter mean value for dip (in degrees): ')
        if len(dipMean) == 0:
            continue  
        else:
            dipMean = float(dipMean) 
           
        dipKappa = input('Enter shape value (kappa) for von Mises distribution (low=dispersed, high=clustered) [5.0]: ')        
        if len(dipKappa) == 0:
            dipKappa = 5.0
        else:
            dipKappa = float(dipKappa) 
 
        #   von Mises circular normal distribution 
        dip = stats.vonmises.rvs(dipKappa, dipMean*np.pi/180, 1/dipKappa, N)
        dip = dip * 180/np.pi 
        dipMin, dipMax, dipMean, dipMedian = getCircStats(dip) 
        dip[np.ix_(dip<0)] = 0.
        dip[np.ix_(dip>90)] = 90.
        
        fig, ax = plt.subplots(figsize=(6,6))
        ndip, bdip, pdip = ax.hist(dip, 20)
        
        sortedbdip = bdip[ndip.argsort()]
        incrbdip = bdip[1] - bdip[0]
        dipMode = sortedbdip[-1] + incrbdip/2.
        print('\n*** Summary statistics for this distribution')
        print('Minimum: %3.2f%s' % (dipMin, degree_sign))
        print('Maximum: %3.2f%s' % (dipMax, degree_sign))
        print('Median:  %3.2f%s' % (dipMedian, degree_sign))
        print('Mode:    %3.2f%s' % (dipMode, degree_sign))
        print('Mean:    %3.2f%s' % (dipMean, degree_sign))
        
        ax.plot([dipMean, dipMean], [0, max(ndip)*1.1], '-r')
        ax.set_xlabel('Dip, degrees')
        ax.set_ylabel('Count')
        ax.set_xlim(0., 180.)
        ax.grid(True)
        ax.set_title('Dip distribution, n=%i' % N)
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            bInput = False 
        
    return dip 

def getSHazDist(N):
    #   function to generate a distribution of N SHmax azimuth values 
    #   assuming von Mises (circular normal) distribution, 
    #   controlled by shape parameter kappa
    bInput = True     
    while bInput:
        
        sHazMean = input('Enter mean value for sHaz (in degrees): ')
        if len(sHazMean) == 0:
            continue 
        else:
            sHazMean = float(sHazMean) 
           
        sHazKappa = input('Enter shape value (kappa) for von Mises distribution (low=dispersed, high=clustered) [5.0]: ')        
        if len(sHazKappa) == 0:
            sHazKappa = 5.0
        else:
            sHazKappa = float(sHazKappa) 
 
        #   von Mises circular normal distribution 
        sHaz = stats.vonmises.rvs(sHazKappa, sHazMean*np.pi/180, 1/sHazKappa, N)
        sHaz = sHaz * 180/np.pi 
        sHazMin, sHazMax, sHazMean, sHazMedian = getCircStats(sHaz) 
        sHaz[np.ix_(sHaz<0)] += 180.
        sHaz[np.ix_(sHaz>180)] -= 180.
        
        fig, ax = plt.subplots(figsize=(6,6))
        nsHaz, bsHaz, psHaz = ax.hist(sHaz, 20)
        
        sortedbsHaz = bsHaz[nsHaz.argsort()]
        incrbsHaz = bsHaz[1] - bsHaz[0]
        sHazMode = sortedbsHaz[-1] + incrbsHaz/2
        print('\n*** Summary statistics for this distribution')
        print('Minimum: %3.2f%s' % (sHazMin, degree_sign))
        print('Maximum: %3.2f%s' % (sHazMax, degree_sign))
        print('Median:  %3.2f%s' % (sHazMedian, degree_sign))
        print('Mode:    %3.2f%s' % (sHazMode, degree_sign))
        print('Mean:    %3.2f%s' % (sHazMean, degree_sign))
        
        ax.plot([sHazMean, sHazMean], [0, max(nsHaz)*1.1], '-r')
        ax.set_xlabel('sigmaH azimuth, degrees')
        ax.set_ylabel('Count')
        ax.set_xlim(0., 180.)
        ax.grid(True)
        ax.set_title('sigmaH azimuth distribution, n=%i' % N)
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            bInput = False    
    
    return sHaz 

def getStressDist(N):
    #   function to generate a distribution of N stress magnitudes 
    #   assuming (symmetric) normal distribution
    #   constrained to be s1 >= s2 >= s3, for all N
    #   returns stress, Nx3 array with principal stresses in columns 
    #       s1 in col 0, s2 in col 1, s3 in col 2 
    bInput = True     
    while bInput:
        
        s1Mean = input('Enter mean value for sigma 1 (most compressive, MPa): ')
        if len(s1Mean) == 0:
            continue 
        else:
            s1Mean = float(s1Mean) 
           
        s1StdDev = input('Enter standard deviation for sigma 1 [10% of mean]: ')
        if len(s1StdDev) == 0:
            s1StdDev = s1Mean * .1   
        else:
            s1StdDev = float(s1StdDev) 
           
        s2Mean = input('Enter mean value for sigma 2 (intermediate, MPa): ')
        if len(s2Mean) == 0:
            continue 
        else:
            s2Mean = float(s2Mean) 
           
        s2StdDev = input('Enter standard deviation for sigma 2 [10% of mean]: ')
        if len(s2StdDev) == 0:
            s2StdDev = s2Mean * .1   
        else:
            s2StdDev = float(s2StdDev) 
           
        s3Mean = input('Enter mean value for sigma 3 (least compressive, MPa): ')
        if len(s3Mean) == 0:
            continue 
        else:
            s3Mean = float(s3Mean) 
           
        s3StdDev = input('Enter standard deviation for sigma_3 [10% of mean]: ')
        if len(s3StdDev) == 0:
            s3StdDev = s3Mean * .1   
        else:
            s3StdDev = float(s3StdDev) 
          
        stress = np.zeros([N,3])    
        nStress = 0     
        while nStress < N: 
        
            #   get random normal picks for s1, s2, s3 
            s1 = stats.norm.rvs(s1Mean, s1StdDev, 1)
            s2 = stats.norm.rvs(s2Mean, s2StdDev, 1)
            s3 = stats.norm.rvs(s3Mean, s3StdDev, 1)
            
            if s1 < s2 or s2 < s3:
                print('Bad stress badger...')
                print(s1, s2, s3)
                continue 
            else:
                stress[nStress, :] = s1, s2, s3  
                nStress += 1 

        fig, ax = plt.subplots(figsize=(12,6))
        ns1, bs1, ps1 = ax.hist(stress[:,0], 20, label='sigma1')
        ns2, bs2, ps2 = ax.hist(stress[:,1], 20, label='sigma2')
        ns3, bs3, ps3 = ax.hist(stress[:,2], 20, label='sigma3')
        
        ax.set_xlabel('Stress, MPa')
        ax.set_ylabel('Count')
#        ax.set_xlim(0., 180.)
        ax.grid(True)
        ax.set_title('Stress distribution, n=%i' % N)
        ax.legend()
        fig.tight_layout() 
        plt.show()
        
        sFinished = input('Satisfied? (Y/N) [Y]: ')
        if 'Y' in sFinished.upper() or len(sFinished) == 0:
            print('Stresses assigned, N=%i' % nStress)
            bInput = False    
    
    return stress

def setStressDist(N, s1Mean, s1StdDev, s2Mean, s2StdDev, s3Mean, s3StdDev):
    #   function to generate a distribution of N stress magnitudes 
    #   assuming (symmetric) normal distribution
    #   constrained to be s1 >= s2 >= s3, for all N
    #   returns stress, Nx3 array with principal stresses in columns 
    #       s1 in col 0, s2 in col 1, s3 in col 2 
    
    stress = np.zeros([N,3])    
    nStress = 0     
    
    while nStress < N: 
    
        #   get random normal picks for s1, s2, s3 
        s1 = stats.norm.rvs(s1Mean, s1StdDev, 1)
        s2 = stats.norm.rvs(s2Mean, s2StdDev, 1)
        s3 = stats.norm.rvs(s3Mean, s3StdDev, 1)
        
        if s1 < s2 or s2 < s3:
#            print('Bad stress combo, skipping...')
#            print(s1, s2, s3)
            continue 
        else:
            stress[nStress, :] = s1, s2, s3  
            nStress += 1 
    
    return stress    

def setStressDist2(N, sVMean, sVStdDev, sHMean, sHStdDev, shMean, shStdDev):
    #   function to generate a distribution of N stress magnitudes 
    #   assuming (symmetric) normal distribution
    #   returns stress, Nx3 array with principal Andersonian stresses in columns 
    #       sV in col 0, sHmax in col 1, shmin in col 2 
    
    stress = np.zeros([N,3])    
    nStress = 0     
    
    while nStress < N: 
    
        #   get random normal picks for s1, s2, s3 
        sV = stats.norm.rvs(sVMean, sVStdDev, 1)
        sH = stats.norm.rvs(sHMean, sHStdDev, 1)
        sh = stats.norm.rvs(shMean, shStdDev, 1)
        
        stress[nStress, :] = sV, sH, sh  
        nStress += 1 
    
    return stress    

def setFaultDist(N, strikeMean, strikeKappa, dipMean, dipKappa):

    strike = np.zeros([N,])
    dip = np.zeros([N,])
    nFault = 0 
    
    while nFault < N:
        
        thisStrike = rad2deg(np.random.vonmises(deg2rad(strikeMean), strikeKappa, 1))
        thisDip = rad2deg(np.random.vonmises(deg2rad(dipMean), dipKappa, 1))

        if thisDip > 90.:
            thisDip = np.abs(180. - thisDip) 
            thisStrike += 180. 

        if thisDip < 0.: 
            thisDip = np.abs(thisDip)
            thisStrike += 180. 
   
        strike[nFault] = thisStrike 
        dip[nFault] = thisDip             
        nFault += 1 

    strikeTrue = np.copy(strike) 
    strike[np.ix_(strike<0.)] += 360.
    strike[np.ix_(strike>360.)] -= 360.
    
    return strike, strikeTrue, dip 

def plotTornadoSf(lows, highs, names, mode):

    widths = highs - lows     
    widths, lows, names = zip(*sorted(zip(widths, lows, names))[::-1])
        
    #   y position for each variable
    ys = range(len(names))[::-1]  # top to bottom
    
    fig, ax = plt.subplots(figsize=(4,4))
    #   plot the bars, one by one
    for y, low, w in zip(ys, lows, widths):
    
        #   each bar is a "broken" horizontal bar chart
        ax.broken_barh([(low, w)], (y - 0.4, 1.),
            facecolors=['white', 'white'],  
            edgecolors=['black', 'black'],
            linewidth=1)
    
        #   display the width as text, in the centre of the bar 
        #x = low + w / 2
        sw = '%3.2f MPa' % w
        ax.text(mode*1.2, y, sw, va='center', ha='center')
    
    #   draw a vertical line down the middle
    plt.axvline(mode, linestyle='--', color='black', lw=0.5)
    
    #   position the x-axis on the top, hide all the other spines (=axis lines)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xlabel('Fracture susceptibility, MPa')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')  
    
    #   make the y-axis display the variables
    plt.yticks(ys, names)    
    plt.ylim(-1, len(names))    

    fig.tight_layout() 
    plt.savefig('pfsTornadoSf.png', dpi=300)
    
    return 

def plotTornadoTsTd(lows, highs, names, mode, f):

    widths = highs - lows     
    widths, lows, names = zip(*sorted(zip(widths, lows, names))[::-1])
        
    #   y position for each variable
    ys = range(len(names))[::-1]  # top to bottom
    
    fig, ax = plt.subplots(figsize=(4,4))
    #   plot the bars, one by one
    for y, low, w in zip(ys, lows, widths):
    
        #   each bar is a "broken" horizontal bar chart
        ax.broken_barh([(low, w)], (y - 0.4, 1.),
            facecolors=['white', 'white'],  
            edgecolors=['black', 'black'],
            linewidth=1)
    
        #   display the width as text, in the centre of the bar 
        #x = low + w / 2
        sw = '%1.2f' % w
        ax.text(mode*1.1, y, sw, va='center', ha='center')
    
    #   draw a vertical line down the middle
    plt.axvline(mode, linestyle='--', color='black', lw=0.5)
    
    #   position the x-axis on the top, hide all the other spines (=axis lines)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if f == 1: 
        ax.set_xlabel('Slip tendency')
    else: 
        ax.set_xlabel('Dilation tendency')
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')  
    
    #   make the y-axis display the variables
    plt.yticks(ys, names)    
    plt.ylim(-1, len(names))    

    fig.tight_layout() 
    if f == 1: 
        plt.savefig('pfsTornadoTs.png', dpi=300)
    else:
        plt.savefig('pfsTornadoTd.png', dpi=300)
        
    return 
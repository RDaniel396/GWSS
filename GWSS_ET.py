#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 09:17:30 2020

@author: richard
"""
import numpy as np
from numpy import array
from numpy.random import choice 
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate
from scipy.stats import norm
import time
start=time.time()

c=299792458 #SI
G=6.674e-11 #SI
SolM=1.9891e30 #Solar mass SI
unit=SolM*G/(c**3) #Conversion of mass to seconds


"""
1. We take the cosmology from the output of class and interpolate the data 
    to create more points in the redshift we are concerned with, z_max<z<z_min.
    
2. Once we have the interpolated the data, we can use the distripution of 
    sources as seen on earth arXiv: 1608.08008, to create simulated data.
    
3. With the simulated data we create our own file "data" to use to determine
    the errors from insturuments and error with observing GW at higher redshift.
"""
##############################################################################
###      SECTION 1 - SIMULATING THE DATA
##Load the data from the output of class
data = np.loadtxt( "lcdm00_background.dat" )
## Determine the range of redhsift we are interested in
zmin=0.01
zmax=6
NoGW=200

z_step=4000  ##Higher steps, means more accurate simulations, but slower code 2000 works well 

## Pull the data from the output of class code:
z=data[:,0] ## redshift
t=data[:,1] ##proper time in terms of Gyrs
d_l=data[:,6] ## luminosity distance in Mpc
d_l=d_l*(3.08567758149137e22)/c       ##convert to seconds [Mpc]/[c]
d_c=data[:,5] ## comoving distance 
d_c=d_c*(3.08567758149137e22)/c       ##convert to seconds [Mpc]/[c]
H=data[:,3]  ## hubble constant
H=H*c/(3.08567758149137e22)           ##convert 1/s        [c]/[Mpc]


## We then interpolate the data we have picked as the implicit funtion is unkown
fH=interpolate.interp1d(z,H, kind='cubic')
fd_c=interpolate.interp1d(z,d_c, kind='cubic')
fd_l=interpolate.interp1d(z,d_l, kind='cubic')
ft=interpolate.interp1d(z,t, kind='cubic')


## We can then chose the redshift range we are interested in
## from this we create new data from the interpolation in the redshift range
Z=np.linspace(zmax,zmin,z_step).reshape(-1,1)
#Z=np.random.uniform(zmin,zmax,z_step).reshape(-1,1)

H_new=fH(Z)
d_c_new=fd_c(Z)
d_l_new=fd_l(Z)
t_new=ft(Z)

## Create a temporary array to hold the new interpolated data
D=[]
D=np.concatenate((Z, H_new, d_c_new, d_l_new, t_new),axis=1)


## From here we then need to use the distripution of sources measured, P:
## Following arXiv: 1608.08008, which uses approximation with a cutoff at z=5

i=0
R=[]
while i<= z_step-1:
    a=Z[i,0]
    if a<=1:
        r=1+2*a
    elif 1<a<=5:
        r=3*(5-a)/4
    else: r=0
    R.append(r)
    i=i+1
R=array([R],dtype=np.float32).T    
P=4*np.pi*(d_c_new**2)*R/(H_new*(1+Z))

## Alternitivley we can use the distribution from arXiv: 1201.3563 with some minor assumptions

# i=0
# rho=[]
# while i<= z_step-1:
#     a=Z[i,0]
#     if a<=1.04:
#         r=10**(-1.82)*(1+a)**(3.28)
#     elif 1<a<4.48:
#         r=10**(-0.724)*(1+a)**(-0.26)
#     else: r=10**(4.99)*(1+a)**(-8)
#     rho.append(r)
#     i=i+1
# rho=array([rho],dtype=np.float32).T    
# P=(4*np.pi/3)* d_c_new**3 * (np.log(t_new/0.2)*0.1/(1+Z))

P_dist=P/sum(P) ## Change this to a distripution of fractions/percentage


data=[]
def simulation(NoGW,data):
    
    ## We can then pick our sources of redshift weighted by the distribution 
    Z_new=choice(Z[:,0],p=P_dist[:,0],size=NoGW)
    
    ## At this point we now have 1000 simulated redshifts sources. 
    ## We then pick out the corresponding interpolated data, to be our simulation
    data=[]
    i=0
    while i<NoGW:
        d=D[D[:,0]==(Z_new[i])]
        data=np.concatenate((data,d[:,0],d[:,1],d[:,2],d[:,3],d[:,4]),axis=0)
        i=i+1
    
    ## All the units are the same as the output from class
    data=np.array(data).reshape((NoGW,5)) ## This now holds our simulated data
    ## c.1=z c.2=H c.3=comov dist c.4=lum dist c.5=time 
    
    
    
    
    ##############################################################################
    ###             SECTION 3 - ADDITIONAL SIMULATED DATA (Without the need for distribution)
    
    
    ## At this point I add other simulated measurements - but these do not need
    ## A distribution, they can be random or called upon from sources in arXiv
    
    ## Add in additional columns to the data to provide 
    ## theta, phi: relative location of the source in the sky and psi: polarisation angle
    
    ##Here I have used a random generator to provide locations and polarisation
    theta=np.random.rand(data[:,1].size,1)*2*np.pi 
    phi=np.random.rand(data[:,1].size,1)*np.pi
    psi=np.random.rand(data[:,1].size,1)*np.pi
    
    data=np.concatenate((data,theta,phi,psi),axis=1)
    ## c.1=z c.2=H c.3=comov dist c.4=lum dist c.5=t c.6=theta c.7=phi c.8=psi 
    
    
    ## Add in additional column to the data for the angle of binary orbital angular momentum, l
    ## Again here I have used a random generator between 0-20 degrees, as per constraint 
    l=np.random.rand(data[:,1].size,1)*20*np.pi/180
    
    data=np.append(data,l,axis=1)
    ## c.1=z c.2=H c.3=comov dist c.4=lum dist c.5=t c.6=theta c.7=phi c.8=psi c.9=l
    
    
    ## Add in additional column of the masses, m1 and m2, of the binary system 
    ## Then compute the total mass, M_tot, symmetric mass ratio, M_sym, and chrip mass, M_c.
    ## All masses are in units of solar mass 
    
    
    ## m_NS: Minimum is 1.17 from arXiv: 1808.02328 Max is 2.9 from arXiv: astro-ph/9608059 
    ## m_BH: maxiumum currently used value from arXiv: 1608.08008
    
    m1=np.random.randint(117,290,size=(data[:,1].size,1))/100
    ## Ratio of BHNS and BNS systems is 0.03 from LIGO, used for our simulation  
    ratio=int(data[:,1].size*0.03)       #forcing integer values so all arrays align
    
    m2_ns=np.random.randint(117,290,size=(data[:,1].size-ratio,1))/100
    m2_bh=np.random.randint(291,1000,size=(ratio,1))/100
    m2=(np.concatenate((m2_ns,m2_bh),axis=0)) 
    np.random.shuffle(m2)  ##Shuffle the masses throughout the redshift
    
    M_tot=m1+m2 ##total mass
    M_c=(1+data[:,0,None])*M_tot[:,0,None]*( (((m1[:,0,None]*m2[:,0,None])/M_tot[:,0,None]**2)) **(3/5)) * unit  ##Chirp mass converted to units of seconds
    print(M_c.shape)
    
    ##############################################################################
    ###             SECTION 3 - SIMULATING THE ERRORS
    
    """
    
    We have two errors to simulate here
    1. The error due to the distance of the GW
    2. The error from the instruments, in this case from the Einstein Telescope, ET. 
    """
    
    
    
    ## Antenna pattern functions observed by ET arXiv:1608.0800
    ## plus, cross are the orientations seen 1,2,3 are the inferometers used 
    def plus(theta,phi,psi):
        return np.sqrt(3)*( 0.5*(1+np.cos(theta)**2) * np.cos(2*phi)*np.cos(2*psi) - np.cos(theta)*np.sin(2*phi)*np.sin(2*psi)) /2 
    def cross(theta,phi,psi):
        return np.sqrt(3)*( 0.5*(1+np.cos(theta)**2) * np.cos(2*phi)*np.cos(2*psi) + np.cos(theta)*np.sin(2*phi)*np.sin(2*psi)) /2 
    F1plus=plus(theta,phi,psi)
    F1cross=cross(theta,phi,psi)
    F2plus=plus(theta,phi+2*np.pi/3,phi)
    F2cross=cross(theta,phi+2*np.pi/3,phi)
    F3plus=plus(theta,phi+4*np.pi/3,phi)
    F3cross=cross(theta,phi+4*np.pi/3,phi)
        
    
    ### Define a funtion to compute the signal to noise ratio, SNR, for any inferometer
    
    def SNR(Fplus,Fcross):
        i=0
        snr=[]
        while i<=data[:,1].size-1:  ## Go through all of the simulated redshifts
            
        
            ##Overal phase of FT based on arXiv: 0903.0338
            #gamma=0.5772156649  #Euler's constant
            
        
            
        
            ## FT of strain/relative distance moved
            Amp=np.sqrt((Fplus[i,0]*(1+np.cos(l[i,0])**2))**2 + (2*Fcross[i,0]*np.cos(l[i,0]))**2) * np.sqrt(5*np.pi/96) * np.pi**(-7/6) * M_c[i,0]**(5/6)/data[i,3]
            if np.isnan(Amp)==True:
                Amp=np.sqrt(-((Fplus[i,0]*(1+np.cos(l[i,0])**2))**2 + (2*Fcross[i,0]*np.cos(l[i,0]))**2)) * np.sqrt(5*np.pi/96) * np.pi**(-7/6) * M_c[i,0]**(5/6)/data[i,3]
                Amp=Amp*1j
                

               
            def PSD(f):   ## One sided power spectral density arXiv: 1009.0206
                f0=(200)
                x=(f)/f0
                S0=(1.449*1e-52)
                p1=(-4.05)
                p2=(-0.69)
                a1=(185.62)
                a2=(232.56)
                b1=(31.18)
                b2=(-64.72)
                b3=(52.24)
                b4=(-42.16)
                b5=(10.17)
                b6=(11.53)
                c1=(13.58)
                c2=(-36.46)
                c3=(18.56)
                c4=(27.43 ) 
                return (S0*(x**(p1) + a1*x**(p2) 
                       +a2*( (1+ b1*x + b2*x**2 + b3*x**3 + b4*x**4 + b5*x**5 +b6*x**6) / (1+ c1*x + c2*x**2 + c3*x**3 +c4*x**4) )  ))
        
            ## Here we state the bounds of frequency that ET can measure
            f_lower=1 ##current lowest detection 
            f_upper=((c**3)/(6**(3/2)*np.pi*G*(1+data[i,0])*M_tot[i,0]*SolM)) ##current higest detection from the last stable orbit freq
        
            
            ## Finally to compute the signal to noise ratio we need to compute 
            ## The intergral shown below and in arXiv: 1009.0206
            ## Simplified intergrand  - due to exp vanishing
            def Integrand(f):
                return( ((f**(-7/3)) * (np.absolute(Amp))**2 ) / PSD(f))
            
            
            I=4*integrate.quad(Integrand,f_lower,f_upper)
            snr.append((I))
            i=i+1  ## Close the while loop and start on the next simulated redshift
            
        snr=array(snr) ## Create an array of the SNR for each simulated redshift
        
        ## Combine the integral error with the integral to provide full SNR error
        snr=snr[:,0,None]+snr[:,1,None]
        return snr
    ##############################################################################
       
    ## At this point we are now ready to use our defined SNR function to produce ET error
        
    
    s2n1=s2n2=s2n3=[] ## We can then produce the SNR for each inferometer
    k=1
    while k<=3:
        if k==1:
            s2n1=(SNR(F1plus,F1cross))
        elif k==2:
            s2n2=(SNR(F2plus,F2cross))
        else:
            s2n3=(SNR(F3plus,F3cross))
        k=k+1
    
    ## We now have our SNR squared for each inferometer 
    ##sum over all and take sqrt to determine overall SNR, s2n
    s2n=np.sqrt((s2n1) + s2n2 + (s2n3))
    
    ## Add the SNR to the data, this allows us to remove any signals, with corresponding data, 
    ## that fall beneth the allowed minimum of 8
    data=np.append(data,s2n,axis=1)
    ## c.1=z c.2=H c.3=comov dist c.4=lum dist c.5=t c.6=theta c.7=phi c.8=psi c.9=l c.10=SNR
    data = data[
        np.logical_not(data[:,9] < 8)]
    return data

data=simulation(NoGW,data) ## Pull the array form the function


## The loop here makes sure that we at least obtain NoGW sources specified
while data[:,0].size<1000:
    print(data[:,0].size)
    data=np.append(data,simulation(NoGW,data),axis=0)
    

## Due to the manipulation of the data to remove SNR<8 use general data to deine terms
## c.1=z c.2=H c.3=comov dist c.4=lum dist c.5=t c.6=theta c.7=phi c.8=psi c.9=l c.10=SNR

## Define the error due to lens, lens_error, and ET, ET_error

## lum_error=0.05*z*lum_dist
lens_error=0.5*data[:,3]*0.066*((1-(1+data[:,0])**(-0.25))/0.25)**1.8
## ET_error= 2*lum_dist/SNR
ET_error=2*data[:,3]/data[:,9]

## Now we can compute the error in luminosity distance
sigma_dl=np.sqrt((ET_error)**2 + lens_error**2).reshape(-1,1)
data=np.append(data,sigma_dl,axis=1)
## c.0=z c.1=H c.2=comov dist c.3=lum dist c.4=t c.5=theta c.6=phi c.7=psi c.8=l c.9=SNR c.10=sigma_dl

## We then shift the luminosity points given by a guassian distribution
## This is done in bins

I=0
lumdist = []
while I<len(sigma_dl):
    lum_b = data[:,3][I] #Create temp bin for guassian curve to be looped over
    sigma_b = data[:,10][I]
    # Determine the guassian, and populate with new guassian shifted data
    mu = lum_b
    std = sigma_b
    lum = np.random.normal(mu,std) 
    lumdist = np.append(lumdist,lum)
    I += 1
# Replace data
data[:,3] = lumdist[:]

plt.figure()
title=data[:,0].size
plt.errorbar(data[:,0],data[:,3]*c/(3.08567758149137e22),data[:,10]*c/(3.08567758149137e22),fmt='.')
plt.xlabel('z')
plt.ylabel('$d_l [Mpc]$')
# plt.title('Number of sources used=' +str(title))

##############################################################################
end=time.time()
print(end-start)

# temp_data=np.transpose([data[:,0],data[:,3],data[:,10]])
# np.savetxt('GWSS_data_LambdaCDM_CPL_compare.dat',temp_data)



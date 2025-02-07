# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime

def freq_SL ():
    w_SL = np.zeros(n_SL_freq)
    for i in range(n_SL_freq):
        w_SL[i] = (1+ 17*i/(n_SL_freq-1))*1000*2*np.pi
    return w_SL
def J_spec_dens(w,tau):
    #w: freq S:order parameter, tau: relaxation time
    return (1-S2)* tau/(1+(w*tau)**2)
def R1_rho_CSA (w_SL, tau):
    R_1= (sigma_w_N**2/135)* (2 * J_spec_dens(w_SL- 2*w_MAS, tau)+ 4 * J_spec_dens(w_SL- w_MAS, tau)+ 4 * J_spec_dens(w_SL+ w_MAS, tau) + 2 * J_spec_dens(w_SL+ 2*w_MAS, tau))
    return R_1
def R1_rho_dip (w_SL, tau):
    R1 = (C_NH_2/60)* (2* J_spec_dens(w_SL+ w_H+ w_MAS, tau)+ 2* J_spec_dens(w_SL+ w_H- w_MAS, tau)+J_spec_dens(w_SL+ w_H+ 2*w_MAS, tau) + J_spec_dens(w_SL+ w_H- 2*w_MAS, tau) +2* J_spec_dens(w_SL- w_H+ w_MAS, tau)+ 2* J_spec_dens(w_SL+ w_H- w_MAS, tau)+J_spec_dens(w_SL- w_H+ 2*w_MAS, tau)+ J_spec_dens(w_SL- w_H- 2*w_MAS, tau))
    return R1
def print_res(w_SL, signal, addressfile):
    f = open(addressfile + "\\simulated R1_rho signal.txt","w")
    res = "\n".join("{} {}".format(x, y) for x, y in zip(w_SL, signal))
    print("#w_H=", w_H/(1000*2*np.pi), file=f)
    print(res, file=f)
    f.close()
    return
def spin_lock_dist (V_SL, tau):
    minim= V_SL*(1-3*k)
    maxim= V_SL* (1+ 0.3*k)
    V = np.zeros (n_dist)
    for i in range (n_dist):
        V[i]= minim + (maxim-minim)*i/(n_dist-1)
    P_inhom= np.exp(-(V_SL-V)**2/(2*(k*V_SL)**2))
    P_inhom = P_inhom/np.sum(P_inhom)
    R1_rho_tot = P_inhom* (R1_rho_CSA (V, tau) + R1_rho_dip (V, tau))
    return np.sum(R1_rho_tot)
def truncated_log_norm (med, sigma, limit_l, limit_u, n):
    lnmed= np.log(med)
    gauss = np.zeros(n+1)
    lnmedtemp = np.zeros(n+1)
    medtemp = np.zeros(n+1)
    for i in range(0,n+1):
        lnmedtemp[i]= (lnmed+ ((i-(n/2))*3*sigma/(n/2)))
        medtemp[i]= np.exp(lnmedtemp[i])
        gauss[i]= np.exp(-0.5*((i-(n/2))*3/(n/2))**2)/(np.sqrt(2*np.pi*(sigma**2)))
        if medtemp[i]<limit_l:
            gauss [i]=0
            medtemp[i]=0
        if medtemp[i]>limit_u:
            gauss [i]=0
            medtemp[i]=0
    return medtemp[medtemp!=0], gauss[gauss!=0]/np.sum(gauss)
def calibration_factor(tau):
    A1= -0.19157
    A2= 1.00676
    x0 = 10.98998
    p= 1.08783
    f = A2+ (A1-A2)/(1+(tau/x0)**p)
    
    return f
def real_fraction():
    tau, amp = truncated_log_norm (T_med, sigma, limit_l, limit_u, 50)
    factor= calibration_factor(tau*10**6)
    factor[factor>1]=1
    fraction = amp/factor
    fraction = fraction/np.sum(fraction)
    return tau, fraction
def final_signal ():
    sim_R1_rho= np.zeros(n_SL_freq)
    for i in range(n_SL_freq):
        R1_rho_mix = np.zeros(len(tau))
        for j in range(len(tau)):
            R1_rho_mix[j]= fraction[j]* spin_lock_dist (w_SL[i], tau[j])
        sim_R1_rho[i] = np.sum(R1_rho_mix)
    return sim_R1_rho
n_dist = 51
n_SL_freq= 500
k = 0.1
C_NH_2 = 5.2*10**9
sigma_w_N = 9.72*2*np.pi *1000 #Hz
w_H = 100 * 1000 *2*np.pi
w_MAS = 18000*2*np.pi #Hz
S2= 0.9927 # order parameter
loc = input("please give me the address of the destination folder")
#################################################
#log normal distribution parameters
T_med = 4.52707*10**-5
sigma = 2.8
limit_l = 3.99*10**-6
limit_u = 5* 10**-3
#################################################
w_SL = freq_SL ()
tau, fraction = real_fraction()
#tau = np.array([8.53309, 93.5777, 1210])*10**-6
#fraction = np.array([0.64665, 0.21029, 0.14306])

signal_R1_rho = final_signal ()
print_res(w_SL/(1000*2*np.pi), signal_R1_rho, loc)
plt.plot(w_SL,signal_R1_rho, 'k.')
plt.yscale('log')
plt.show
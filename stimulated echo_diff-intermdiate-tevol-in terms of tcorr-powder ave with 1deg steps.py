# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 13:17:48 2023

@author: Farhad
"""

import numpy as np
import os
from datetime import datetime
from scipy.optimize import curve_fit
# Function to print the powder-averaged correlation function to a relevant folder
def print_result(output_address, c_t_subtr, t_evol, tau_corr, No_stepsize_in1us):
    multi_address_file= output_address +"\subtraction correlation function"
    if not os.path.exists(multi_address_file):
            os.makedirs(multi_address_file)
    f = open(multi_address_file + "\corr_func-powder_ave-stimecho-t_corr-"+ str(tau_corr)+"-t-evol-"+ str(t_evol)+".txt","w")
    print("#data provided at", start_time, file= f)
    print("#Q_c=40kHz-time is in us", file= f)
    res = "\n".join("{} {}".format(x, y) for x, y in zip(corr_time/No_stepsize_in1us, c_t_subtr))
    print(res, file=f)
    f.close()
# Function to read the trajectory file (saved in npy format), rotate the trajectory according to angles_list,
# and save it in the folder to avoid running the rotation part again. Calls the powder averaged correlation function calculator.
def data_preparation(loc) :    
    for file in os.listdir(loc):
        if file.endswith(".npy"):
            inputaddress = os.path.join(loc, file)
            namefile = os.path.join(file)
            namefile = namefile[:-4]
            addressfile = loc[:] + "\\" + namefile
            if not os.path.exists(addressfile):
                os.makedirs(addressfile)
            # READING THE TEXT FILE AND COPY IT TO DATA
            data1 = np.load(inputaddress)
            theta_0deg = data1[:,1]
            theta_rot = np.zeros(len(theta_0deg))
            global angles_address
            angles_address = addressfile + "\\angles list"
            if not os.path.exists(angles_address):
                os.makedirs(angles_address)
                for i in range (len(angles_list)):
                    theta_rot = rotation_axis(data1, angles_list[i]*rad)
                    angle_output = angles_address + "\\" +str(angles_list[i])
                    np.save(angle_output, theta_rot)
            tau_relax_step= calc_relaxation_time (data1)
            #it makes different folders to save the correlation function in different folders according to the Quadrupolar coupling
            for sigma in sigma_list:
                address_qudro = addressfile + "\\" + str(sigma)+ "kHz"
                if not os.path.exists(address_qudro):
                    os.makedirs(address_qudro)
                #from kHz to Hz
                sigma = sigma *1000
                for tau_corr in tau_corr_list:
                    for tt in t_evol_list:
                        powder_ave_corr_func (tt, tau_corr, tau_relax_step, address_qudro, sigma)

    return
def bi_exp_decay (t, A1, tau1, A2, tau2, y0):
    return A1*np.exp(-t/tau1)+ A2*np.exp(-t/tau2)+ y0
def single_exp_decay (t, A0, tau, y0):
    return A0*np.exp(-t/tau)+ y0
# Function to rotate the symmetry axis of the trajectory. 
# The original trajectory is simulated around the Z axis and is rotated according to the angles list.
def rotation_axis(data, beta):
    theta_old = data[:,1]*rad
    fi_old = data[:,2]*rad
    x_old = np.sin(theta_old)*np.cos(fi_old)
    z_old = np.cos(theta_old)
    y_old = np.sin(theta_old)*np.sin(fi_old)
    z_new = -x_old*np.sin(beta) + z_old*np.cos(beta)
    x_new = x_old*np.cos(beta) + z_old*np.sin(beta)
    theta_new = np.arccos(z_new)/rad
    return theta_new
# Function to calculate the powder averaged correlation function based on newly rotated trajectories
def powder_ave_corr_func (t_evol, tau_corr, relax_time_stepsize, output_address, sigma):
    corr_func_diff_sym_axis = np.zeros(shape= (len(corr_time), len(angles_list)))
    No_stepsize_in1us= relax_time_stepsize/tau_corr
    i=0
    sum_sin = 0
    for file in os.listdir(angles_address):
        if file.endswith(".npy"):
            ang_address = os.path.join(angles_address, file)
            namefile = os.path.join(file)
            namefile = namefile[:-4]
            ang = int(namefile)
            theta_rotated = np.load(ang_address)
            corr_func= calc_corr_func(theta_rotated, t_evol, output_address, relax_time_stepsize, tau_corr, sigma)
            corr_func_diff_sym_axis[:,i] = corr_func *np.sin(ang*rad)
            sum_sin = sum_sin + np.sin(ang*rad)
            i+=1

    correlataion_function_multi = np.sum(corr_func_diff_sym_axis, axis=1)/sum_sin
    print_result(output_address, correlataion_function_multi, t_evol, tau_corr, No_stepsize_in1us)
# Function to calculate the relaxation time of the motion, used for defining the time axis of the final correlation function    
def calc_relaxation_time (data):
    c_t= np.ones(len(corr_time))
    skip = 200
    for j in range(len(corr_time)):
        i=0
        theta= data[:,1]*rad
        fi = data[:,2]*rad
        t = corr_time[j]
        theta_i= theta[i::skip]
        theta_it = theta[i+t::skip]
        fi_i = fi[i::skip]
        fi_it = fi[i+t::skip]
        theta_ini = theta_i[:len(theta_it)]
        fi_ini = fi_i[:len(fi_it)]
        x = np.cos(theta_ini)*np.cos(theta_it)+ np.sin(theta_ini)*np.sin(theta_it)*np.cos(fi_ini-fi_it)
        cos_2 = (np.average(x**2))
        c_t[j]= 0.5*(3*cos_2 -1)
    fittedParameters, pcov= curve_fit(single_exp_decay, corr_time[:-10], c_t[:-10])
    print (fittedParameters)
    return fittedParameters[1]
# Function to calculate the effect of the intermediate motion in each block of evolution period/time
def fi_intermediate_motion (index, t_evol_step):
    step = t_evol_step//200
    if step == 0:
        step= 1
    return np.average(3*((np.cos(theta[index:index+t_evol_step:step]))**2)-1)
# Function to calculate the correlation function of a single trajectory in a specific orientation with respect to the Z axis
def calc_corr_func (theta_rotated, t_evol_us, output_address, relax_time_stepsize, tau_corr, sigma):
    global theta
    theta= theta_rotated*rad
    i=0
    No_stepsize_in1us= relax_time_stepsize/tau_corr
    t_evol_step= int(round(t_evol_us*No_stepsize_in1us))
    
    #c_t_multi= np.ones(len(corr_time))
    c_t_subtr= np.ones(len(corr_time))
    #c_t_sum= np.ones(len(corr_time))
    for j in range(len(corr_time)):
        t = corr_time[j]
        No_ave_points = (len(theta) -t -2*t_evol_step-1)//skip_value
        #initial and end points of the calculation of correlation function are first saved in an array to avoid using any loop
        index_ini = np.arange(No_ave_points)* skip_value
        index_end = index_ini + t+ t_evol_step
        
        fi_1 = np.pi*sigma*t_evol_us* 10**-6* v_fi_intermediate_motion (index_ini, t_evol_step)
        fi_2 = np.pi*sigma*t_evol_us* 10**-6* v_fi_intermediate_motion (index_end, t_evol_step)
        #x_1 = np.cos(fi_1)*np.cos(fi_2)
        x_2 = np.cos(fi_1-fi_2)
        #x_3 = np.cos(fi_1+fi_2)
        #c_t_multi[j]= np.average(x_1)
        c_t_subtr[j]= np.average(x_2)
        #c_t_sum[j]= np.average(x_3)
        return c_t_subtr

    
loc = input("please give me the address of the folder which contains trajectory data")
start_time = datetime.now()
rad= np.pi/180
#sigma_list is the list of quadropular coupling constants. it is possible to study the effect of different quadrupolar coupling constant
sigma_list = np.array([40])
#x is the time in logarithmic scale but in the unit of steps. the trajectory is in the unit of steps. by assigning a defined time to each step it becoms in the unit of real time. for example ms or us
x = np.logspace(1,6,110)
#corr_time is the time exactly as x. since np.logspace does not contain zero, I added 0 to 10 manually
corr_time = np.append([0,1,2,3,5,7], x.astype(int))
#after considering a sets of trajectory points starting from i the next round of calculation of correlation function starts from i+skip_value
skip_value = 800
# use vectrorized version of the following functions t avoid using for loops and mske calculation much faster
v_rotation_axis = np.vectorize(rotation_axis)
v_fi_intermediate_motion = np.vectorize(fi_intermediate_motion)
#tau_corr_list is the array containing of different correlation times of the motion. 
tau_corr_list = np.array([200])
#angles_list is an array containg angles that powder averaging is applied in these angles 
angles_list = np.arange(0,91,1)
#t_evol_list is an array containing different evolution times that correlation function is calculated based on those.
t_evol_list= [11]
data_preparation(loc)

end_time = datetime.now()
print('Duration: {}'.format(end_time - start_time))
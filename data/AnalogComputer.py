# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 12:28:16 2017

@author: William Ngana
"""
from scipy.stats import linregress
from pathlib import Path
import numpy as np
import scipy.signal as sc
import scipy.optimize as opt
import scipy.interpolate as si
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from pathlib import Path
import pandas as pd
 #you can also use df['column_name']

# get the data for the summer
def findDataSum(directory,filename):
    df = pd.read_csv(filename,skiprows=[1], dtype = np.float64)
    time = df["x-axis"]
    v_1 = df["1"]
    v_2 = df["2"]
    v_0 = df["4"]
    return time, v_1, v_2, v_0

#get data for differ and inter
def findDataDiff(directory,filename):
    df = pd.read_csv(filename,skiprows=[1], dtype = np.float64)
    time = df["x-axis"]
    v_1 = df["1"]
    v_0 = df["4"]
    return time, v_1, v_0  

#get data for Damped Harmoic Motion
def findDataDamped(directory,filename):
    df = pd.read_csv(filename,skiprows=[1], dtype = np.float64)
    time = df["x-axis"]
    v_1 = df["1"]
    v_2 = df["2"]
    v_3 = df["3"]
    v_4 = df["4"]
    return time, v_1, v_2, v_3, v_4 
    
#numerical integraation with Reimann sum
def Reimann(x,y): # x and y are arrays of the same length and returns a litst
    ans =[]      
    delta = [x[0]*x[1] for x in zip(np.diff(time),v_1) ]             
    for i in range(len(delta)):
        ans.append(sum(delta[:i])*-1)
    return ans
        

    
if __name__ == '__main__':
    scope_number = 0
    filename= "scope_" + str(scope_number) +".csv"
    dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    '''
    for i in range(12):
        scope_number = i
        filename= "scope_" + str(scope_number) +".csv"
        time, v_1, v_2, v_0 = findDataSum(dir_path,filename)
        plt.plot(time,v_1, color = 'blue')
        plt.plot(time,v_2, color = 'red')
        plt.plot(time,v_0, color = 'green')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        axes.set_xlim([-.25,.25])
        axes.set_ylim([-.5,.5])
        plt.grid()
        green_patchr = mpatches.Patch(color='green', label='v0')
        blue_patchr = mpatches.Patch(color='blue', label='v_1')
        red_patchr = mpatches.Patch(color='red', label='v_2')
        plt.legend(handles=[green_patchr,blue_patchr,red_patchr])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Summing/scope_'+ str(scope_number)+'raw.png')
        plt.show()
        
        plt.plot(time,v_0, color = 'green')
        sv_1 = v_1*(103.78/104.24)
        sv_2 = v_2*(103.78/110.6)
        sum_plot = [-sum(x) for x in zip(sv_1,sv_2)]
        delta = [x[1]-x[0] for x in zip(sum_plot, v_0)]
        plt.plot(time, sum_plot , color = 'red')
        plt.plot(time, delta, color = 'blue')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        axes.set_xlim([-.25,.25])
        axes.set_ylim([-.5,.5])
        plt.grid()
        green_patch = mpatches.Patch(color='green', label='V0')
        blue_patch = mpatches.Patch(color='blue', label='delta')
        red_patch = mpatches.Patch(color='red', label='-(v_1+v_2)')
        plt.legend(handles=[green_patch,blue_patch,red_patch])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Summing/scope_'+ str(scope_number)+'.png')
        plt.show()
        '''
    #Diff
    '''
    for i in range(15,22):
        scope_number = i
        filename= "scope_" + str(scope_number) +".csv"
        time, v_1, v_0 = findDataDiff(dir_path,filename)
        plt.plot(time,v_1, color = 'blue')
        plt.plot(time,v_0, color = 'green')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        #axes.set_xlim([-.25,.25])
        #axes.set_ylim([-.5,.5])
        plt.grid()
        green_patchr = mpatches.Patch(color='green', label='v0')
        blue_patchr = mpatches.Patch(color='blue', label='v_1')
        plt.legend(handles=[green_patchr,blue_patchr])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Differentiation/scope_'+ str(scope_number)+'raw.png')
        plt.show()
        
        der = sc.medfilt(-np.diff(v_1)*(103.78*10**3*1.0254*10**(-7))/.0005,25)
        #delta = [x[1]-x[0] for x in zip(der, v_0)]
        plt.plot(time[1:],der, color = "red")
        plt.plot(time,v_0, color = 'green')
        #plt.plot(time[1:], delta, color = 'blue')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        #axes.set_xlim([-.25,.25])
        #axes.set_ylim([-.5,.5])
        plt.grid()
        green_patchr = mpatches.Patch(color='green', label='v0')
        #blue_patchr = mpatches.Patch(color='blue', label='delta')
        red_patchr = mpatches.Patch(color='red', label='Der')
        plt.legend(handles=[green_patchr,red_patchr])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Differentiation/scope_'+ str(scope_number)+'.png')
        plt.show()'''
    
    
    #integration
    '''
    for i in range(24,26):
        scope_number = i
        filename= "scope_" + str(scope_number) +".csv"
        time, v_1, v_0 = findDataDiff(dir_path,filename)
        plt.plot(time,v_1, color = 'blue')
        plt.plot(time,v_0, color = 'green')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        #axes.set_xlim([-.25,.25])
        #axes.set_ylim([-.5,.5])
        plt.grid()
        green_patchr = mpatches.Patch(color='green', label='v0')
        blue_patchr = mpatches.Patch(color='blue', label='v_1')
        plt.legend(handles=[green_patchr,blue_patchr])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Integration/scope_'+ str(scope_number)+'raw.png')

        plt.show()
        Int = Reimann(time,v_1)
        sInt = np.array(Int)*((97.05*110.6)*10**(-5))**(-1)
        plt.plot(time[:-1],sInt, color = "blue")
        plt.plot(time,v_0, color = 'green')
        plt.xlabel('Time(s)')
        plt.ylabel('VoltS(V)')
        axes =plt.gca()
        #axes.set_xlim([-.25,.25])
        #axes.set_ylim([-.5,.5])
        plt.grid()
        green_patchr = mpatches.Patch(color='green', label='v0')
        blue_patchr = mpatches.Patch(color='blue', label='Intergral')
        plt.legend(handles=[green_patchr,blue_patchr])
        plt.savefig('C:/Users/William Ngana/Desktop/360L A/Experiment 17/graphs/Integration/scope_'+ str(scope_number)+'.png')

        plt.show()
        '''
    
    # Exponetial
    
    for i in range(26,32):
        scope_number = i
        filename= "scope_" + str(scope_number) +".csv"
        cof = [0.0,2.0,4.0,6.0,8.0,9.91]
        V_0 =2.28 #V  
        R3 = 104.24*10**3
        R_0 = 9.91*10**3
        R = cof[i-26]*10**3
        Beta = (R/R_0)  
        C = 97.9*10**(-9)
        alpha = (Beta/R3*C)
        time, v_1, clock, v_0 = findDataSum(dir_path,filename)
        V = [-V_0*np.exp(-alpha*((x+0.504))) for x in time]
        V_n = [x*R3*C for x in v_0]
        derV = np.diff(V)/np.diff(time)
        plt.plot(time,v_1, color = 'red')
        #plt.plot(time,v_0, color = 'blue')
        plt.plot(time,V, color ='blue')
        axes =plt.gca()
        #axes.set_xlim([-.5,-.4])
        #axes.set_ylim([-2.8,0])
        plt.show()
    
    #Damped Harmonic Motion    
    '''
    for i in range(36,37):
        scope_number = i
        filename= "scope_" + str(scope_number) +".csv"
        V_0 =2.28 #V  
        R3 = 104.24*10**3
        R_0 = 9.91*10**3
        R = 2.00*10**3
        C = 97.9*10**(-9)
        alpha = (R/(R_0*R3*C))
        #r = 0.5*(-alpha+ np.sqrt(alpha*alpha - 4*(1.0/(R3*C)**2)))
        time, v_1, v_2, v_3, v_4 = findDataDamped(dir_path,filename)
        plt.plot(time,v_1, color = 'green')
        plt.plot(time,v_2, color = 'blue')
        plt.plot(time,v_3, color = 'red')
        plt.plot(time,v_4, color = 'black') 
        axes =plt.gca()
        axes.set_xlim([-4,-2])
        #axes.set_ylim([-.5,.5])
        plt.show()
        V = [V_0*np.exp(-alpha/2.0*((x))) for x in time]
        derV = np.diff(V)/np.diff(time)     
        a = np.diff(derV)/np.diff(time[1:])
        b = [-(x)/(R3*C) for x in derV]
        c = [x/(R*C)**2 for x in V]
        d = [-x/(R*C)**2 for x in V]
        plt.plot(time[2:],a, color = 'green')
        plt.plot(time[1:],b, color = 'blue')
        plt.plot(time,c, color = 'red')
        plt.plot(time,d, color = 'black')        
        axes =plt.gca()
        #axes.set_xlim([0,1])
        #axes.set_ylim([-4,3])
        plt.show()
    '''
    

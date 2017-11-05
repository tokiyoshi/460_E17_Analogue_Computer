# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 12:28:16 2017

@author: William Ngana
"""
from scipy.stats import linregress
from pathlib import Path
import numpy as np
import shutil
import scipy.signal as sc
import scipy.optimize as opt
import scipy.interpolate as si
from scipy.optimize import leastsq
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
from pathlib import Path
import pandas as pd

R_1 = 103.78*10**3
R_2 = 110.6*10**3
R_3 = 104.24*10**3
R_4 = 99.15*10**3
R_5 = 99.29*10**3
R_6 = 99.24*10**3
C_1 = 102.54*10**(-9)
C_2 = 97.9*10**(-9)
C_3 = 97.05*10**(-9)

V_0 = 2.28

R_0 = 9.91*10**3

#you can also use df['column_name']

# get the data for the summer
def findDataSum(filename):
    df = pd.read_csv(filename, skiprows=[1], dtype=np.float64)
    time = np.array(df["x-axis"])
    v_1 = np.array(df["1"])
    v_2 = np.array(df["2"])
    v_0 = np.array(df["4"])
    return time, v_1, v_2, v_0

#get data for differ and inter
def findDataDiff(filename):
    df = pd.read_csv(filename,skiprows=[1], dtype = np.float64)
    time = np.array(df["x-axis"])
    v_1 = np.array(df["1"])
    v_0 = np.array(df["4"])
    return time, v_1, v_0  

#get data for Damped Harmoic Motion
def findDataDamped(filename):
    df = pd.read_csv(filename, skiprows=[1], dtype = np.float64)
    time = np.array(df["x-axis"])
    v_1 = np.array(df["1"])
    v_2 = np.array(df["2"])
    v_3 = np.array(df["3"])
    v_4 = np.array(df["4"])
    return time, v_1, v_2, v_3, v_4 
    
#numerical integraation with Reimann sum
def Reimann(x,y): # x and y are arrays of the same length and returns a litst
    ans = []
    delta = [x[0]*x[1] for x in zip(np.diff(time),v_1) ]             
    for i in range(len(delta)):
        ans.append(sum(delta[:i])*-1)
    return ans

def boot_plots(x, y_s, labels, line_width = None, x_range = None, y_range = None, colours = None):
    colours_all = ['green', 'blue', 'red', 'black']
    if line_width is None:
        line_width = [0.5, 0.5, 0.5, 0.5]
    if colours is None:
        colours = colours_all[:len(y_s)]
    if x_range is None:
        x_range = [np.min(x), np.max(x)]
    legend_entries = []
    for y, label, colour, lw in zip(y_s, labels, colours, line_width):
        plt.plot(x, y, color=colour, lw=lw)
        legend_entries.append(mpatches.Patch(color=colour, label=label))
    plt.legend(handles=legend_entries)
    axes = plt.gca()
    axes.set_xlim(x_range)
    if y_range is not None:
        axes.set_ylim(y_range)
    plt.xlabel('Time(s)')
    plt.ylabel('Volts(V)')
    plt.grid()
    return axes

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

if __name__ == '__main__':
    summing = False
    differentiate = False
    intergrate = False
    exponential = False
    HO = False
    DHO = True

    scope_number = 0
    filename = "scope_" + str(scope_number) +".csv"
    dir_path = Path(os.path.dirname(os.path.realpath(__file__)))
    graphs_dir = dir_path.parent.joinpath('Graphs')
    if summing:
        save_dir = graphs_dir.joinpath('Summing')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)
        for scope_number in range(12):
            time, v_1, v_2, v_0 = findDataSum("scope_%s.csv" % scope_number)
            boot_plots(time, [v_1, v_2, v_0], ['Input 1', 'Input 2', 'Summed'], y_range=[-.5, .5])
            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            scaled_v_1 = v_1*(R_1/R_2)
            scaled_v_2 = v_2*(R_1/R_3)
            sum_plot = [-sum(x) for x in zip(scaled_v_1, scaled_v_2)]
            delta = [x[1]-x[0] for x in zip(sum_plot, v_0)]
            #print('Median of scope %s is: %s with var: %s' % (scope_number, np.median(delta), np.var(delta)))
            boot_plots(time, [v_0, sum_plot, delta], ['Summed', 'Calculated', 'Delta'], y_range=[-.5, .5])
            plt.savefig(str(save_dir.joinpath('scope_%s.png' % scope_number)))
            plt.close()

    #Diff
    if differentiate:
        save_dir = graphs_dir.joinpath('Differentiation')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)
        for scope_number in range(15, 22):
            time, v_1, v_0 = findDataDiff("scope_%s.csv" % scope_number)
            boot_plots(time, [v_1, v_0], ['Input', 'Output'])
            calc_differ = np.append((np.diff(v_1) * ((-R_1 * C_1)))/.0005,0) # adding in an empty entry to keep length consistant
            calc_differ_filtered = sc.medfilt(calc_differ, 41)

            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            boot_plots(time, [v_1], ['Input'], colours=['green'])
            plt.savefig(str(save_dir.joinpath('scope_%sraw_input.png' % scope_number)))
            plt.close()

            boot_plots(time, [v_0], ['Output'], colours=['blue'])
            plt.savefig(str(save_dir.joinpath('scope_%sraw_out.png' % scope_number)))
            plt.close()

            boot_plots(time, [calc_differ, calc_differ_filtered], ['Calculated', 'Filtered'], line_width=[.25, 1], colours=['blue', 'green'])
            plt.savefig(str(save_dir.joinpath('scope_%sdiff_filtered.png' % scope_number)))
            plt.close()

            boot_plots(time, [calc_differ_filtered, calc_differ], ['Filtered', 'Calculated'], line_width=[1, .25], colours=['green', 'blue'])
            plt.savefig(str(save_dir.joinpath('scope_%sdiff_filtered_calc_front.png' % scope_number)))
            plt.close()

            boot_plots(time, [calc_differ_filtered, v_0], ['Calculated', 'Output'], colours=['red', 'green'])
            #delta = [x[1]-x[0] for x in zip(der, v_0)]
            plt.savefig(str(save_dir.joinpath('scope_%s_calc_behind.png' % scope_number)))
            plt.close()

            boot_plots(time, [v_0, calc_differ], ['Output', 'Calculated'], colours=['green', 'red'])
            # delta = [x[1]-x[0] for x in zip(der, v_0)]
            plt.savefig(str(save_dir.joinpath('scope_%s_calc_front.png' % scope_number)))
            plt.close()

    #integration
    if intergrate:
        save_dir = graphs_dir.joinpath('Integration')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)
        for scope_number in range(22, 26):
            time, v_1, v_0 = findDataDiff("scope_%s.csv" % scope_number)
            boot_plots(time, [v_1, v_0], ['Input', 'Output'])
            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            Int = Reimann(time,v_1)
            sInt = np.array(np.append(Int, 0))*(C_3*R_2)**(-1)
            boot_plots(time, [sInt, v_0], ['Numerical', 'Output'])
            plt.savefig(str(save_dir.joinpath('scope_%s.png' % scope_number)))
            plt.close()
    
    # Exponetial
    if exponential:
        save_dir = graphs_dir.joinpath('Exponential')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)

        for scope_number in range(26, 32):
            cof = [0.0, 2.0, 4.0, 6.0, 8.0, 9.91]
            Rf = R_3
            R_0 = 9.91*10**3
            R_pot = cof[scope_number-26]*10**3
            C = C_2

            #Beta = (R_pot / (R_pot + Rf))

            Beta = R_pot/R_0
            alpha = (Beta/((Rf+R_pot)*C))

            time, v_1, clock, v_0 = findDataSum("scope_%s.csv" % scope_number)

            V = [-V_0*np.exp(-alpha*(x + abs(np.min(time)))) for x in time]
            V_n = [x*Rf*C for x in v_0]

            derV = np.append(np.diff(V)/np.diff(time), 0)

            boot_plots(time, [v_1, v_0, clock], ['v_1', 'v_0', 'Clock'])
            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            boot_plots(time, [v_1, derV], ['v_1', 'Calculated'], y_range=[-2.8, 0])
            #x_range=[-.5, -.4]
            plt.savefig(str(save_dir.joinpath('scope_%s_calc.png' % scope_number)))
            plt.close()

            #x_range=[-.5, -.4]
            boot_plots(time, [v_1, V], ['v_1', 'Calculated_Voltage'], y_range=[-2.8, 0])
            plt.savefig(str(save_dir.joinpath('scope_%s_calc.png' % scope_number)))
            plt.close()
    
    #Harmonic Motion
    if HO:
        save_dir = graphs_dir.joinpath('HO')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)

        for scope_number in range(36, 37):
            time, v_1, v_2, v_3, v_4 = findDataDamped("scope_%s.csv" % scope_number)
            Rf = R_3
            C = C_2
            R = 2.00 *10**3 # INCORRECT CHANGE

            alpha = (R/(R_0*Rf*C))
            #r = 0.5*(-alpha+ np.sqrt(alpha*alpha - 4*(1.0/(R3*C)**2)))

            boot_plots(time, [v_1, v_2, v_3, v_4], ['a', 'b', 'c', 'd'], x_range=[-4,-2])
            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            V = [V_0*np.exp(-alpha/(2.0*x)) for x in time]
            derV = np.diff(V)/np.diff(time)
            a = np.append(np.diff(derV)/np.diff(time[1:]), [0, 0])
            b = np.append([-(x)/(Rf*C) for x in derV], 0)
            c = [x/(R*C)**2 for x in V]
            d = [-x/(R*C)**2 for x in V]

            boot_plots(time, [a, b, c, d], ['calc_a', 'calc_b', 'calc_c', 'calc_d'], x_range=[-4, -2])
            plt.savefig(str(save_dir.joinpath('scope_%scalc.png' % scope_number)))
            plt.close()
    if DHO:
        save_dir = graphs_dir.joinpath('DHO')
        shutil.rmtree(save_dir, ignore_errors=True)  # Avoiding errors if folder already is gone
        save_dir.mkdir(parents=True)
        freq = [4.9, 10, 12, 14, 16, 18, 20, 25, 50]
        amplitude_list = []
        for scope_number in range(38, 47):
            time, v_1, v_2, v_3, v_4 = findDataDamped("scope_%s.csv" % scope_number)

            data = v_2[400:600]
            data_mid, data_max, data_min = np.median(data), np.max(np.sort(data)[:len(v_2)-20]), np.min(np.sort(data)[20:])
            amplitude = ( np.abs(data_mid - data_max) + np.abs(data_mid - data_min) )/2  # sketchy way to get amplitude without fitting data
            amplitude_list.append(amplitude)

            boot_plots(time, [v_1, v_2, v_3, v_4], ['a', 'b', 'c', 'd'])
            plt.savefig(str(save_dir.joinpath('scope_%sraw.png' % scope_number)))
            plt.close()

            boot_plots(time, [v_2], ['b'], colours=['blue'])
            plt.savefig(str(save_dir.joinpath('scope_%sv_2.png' % scope_number)))
            plt.close()

        plt.scatter(freq, amplitude_list)
        y = np.array(amplitude_list)
        x = np.array(freq)
        n = len(x)  # the number of data
        mean = sum(x * y) / n  # note this correction
        sigma = sum(y * (x - mean) ** 2) / n  # note this correction
        popt, pcov = curve_fit(gaus, x, y, p0=[1, mean, sigma])

        FWHM = 2.355 * np.abs(popt[2])

        plt.plot(x, y, 'b+:', label='data')
        x_fill = np.linspace(0, 50, 100)
        #plt.plot(x, gaus(x_fill, *popt), 'ro:', label='fit')
        plt.plot(x_fill, gaus(x_fill, *popt), color='red', label='fit')
        plt.legend()
        plt.title('Peak = %6.2f FWHM = %6.2f' % (popt[1], FWHM))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Amplitude (V)')
        plt.savefig(str(save_dir.joinpath('q_value.png')))



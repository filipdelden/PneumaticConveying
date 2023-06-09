#!/usr/bin/env pythons
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import csv
import tkinter as tk


class PressureCalc():
    def __init__(self):
        # System parameters
        sysParameters = {'Diam': ('Pipe Diameter (m)', 0.1),
                     'H_Length': ('Horizontal Length excl bends (m)', 8.156),
                     'V_Length': ('Vertical Length excl bends (m)', 17.541),
                     'nbH': ('Number of Horizontal Bends',4),
                     'nbV': ('Number of Vertical Bends',1),
                     'Br': ('Bend Ratio', 7.5),
                     'ksiB': ('Bend friction coefficient', 0.5),
                     'g': ('Gravitational acceleration (m/s2)', 9.81),
                     'R': ('Gas constant (kJ/(kg*K))',0.287),
                     'T': ('Temperature (K)',293),
                     'pres': ('Inlet pressure (Pa)',400000),
                     'dynVisc': ('dynamic viscosity (Pa*s)',18.1e-6),
                     'rho_p': ('Particle Density (kg/m3)',2000),
                     'dp': ('Particle Diameter (um)',1e-5),
                     'Mfs': ('Mass flow solids (kg/s)',0.3,5.3,11),
                     'Mfa': ('Mass flow air (kg/s)',0.05, 1.5, 10000)}
        
        self.pressureDrop = self.pressureDropCalc(self,sysParameters)    
        

    def pressureDropCalc(self,**sysParam):
        pressureDrop=[]
        mindeltaPx=[]
        mindeltaPy=[]
        mindeltaSLR=[]

        lenMfa = len(sysParam['Mfa'][1])
        lenMfs = len(sysParam['Mfs'][1])

        Ac = np.pi*(sysParam['Diam'][1]**2) / 4  # Cross-sectional area (m^2)
        delta_L_bend = 2*np.pi * sysParam['Br'][1] * sysParam['Diam'][1]/4  # Length of 90deg bend(m)
        deltaL = (sysParam['H_Length'][1]) + (sysParam['V_Length'][1]) + delta_L_bend*(sysParam['nbH'][1] + sysParam['nbV'][1])  # Pipe length (m)
        deltaZ = sysParam['Br'][1] * sysParam['Diam'][1] * sysParam['nbV'][1] + sysParam['V_Length'][1]  # Pipe elevation change (m)

        sysParam['deltaL'] = ('Pipe length (m)',deltaL)

        for i in range(lenMfs):
            tempList = []
            for j in range(lenMfa):
                mfs = sysParam['Mfs'][1][i]
                mfa = sysParam['Mfa'][1][j]

                rho_air = (sysParam['pres'][1]/(sysParam['R'][1]*1000*sysParam['T'][1])+1.89)/2  # kg/m3 Average density
                Vfa = mfa/(rho_air)  # Average volumetric flow rate air
                va = Vfa/Ac  # Average velocity air
                slr = mfs/mfa  # Solids loading ratio

                # solid and air friction factor - lambdaS/L
                Fri = va/(np.sqrt(sysParam['g'][1]*sysParam['Diam'][1]))  # Froude number
                lambdaS = 15.62 * slr**-0.6 * Fri**-1.93
                lambdaL = 0.316/(va*sysParam['Diam'][1]*rho_air/sysParam['dynVisc'][1])**0.25
                velRatio = 1-0.008*((1000*sysParam['dp'][1])**(0.3))*(sysParam['rho_p'][1]*0.65)**0.5

                # Pressure drops
                deltaPair = lambdaL*rho_air*va**2*deltaL/(2*sysParam['Diam'][1])
                deltaPacc = slr*va*rho_air*0.9*va
                deltaPsolids = slr*lambdaS*rho_air*va**2*deltaL/(2*sysParam['Diam'][1])
                deltaPbends = (1+slr)*(sysParam['nbH'][1]+sysParam['nbV'][1])*sysParam['ksiB'][1]*rho_air*va**2 / 2
                deltaPvertical = slr*rho_air*sysParam['g'][1]*deltaZ/velRatio

                deltaPtot = deltaPacc+deltaPair+deltaPsolids+deltaPbends+deltaPvertical

                tempList.append(deltaPtot)

            min_value = min(tempList)
            min_index = tempList.index(min_value)

            mindeltaPx.append(sysParam['Mfa'][1][min_index])
            mindeltaSLR.append(sysParam['Mfs'][1][i]/sysParam['Mfa'][1][min_index])
            mindeltaPy.append(min_value/deltaL)

            pressureDrop.append(tempList)

        return pressureDrop


    def plotOnlyPressure(self,calcData,**sysParam):
        delta_L_bend = 2*np.pi * sysParam['Br'][1] * sysParam['Diam'][1]/4  # Length of 90deg bend(m)
        deltaL = (sysParam['H_Length'][1]) + (sysParam['V_Length'][1]) + delta_L_bend*(sysParam['nbH'][1] + sysParam['nbV'][1])

        Mfa = np.linspace(sysParam['Mfa'][1], sysParam['Mfa'][2], sysParam['Mfa'][3])
        Mfs = np.linspace(sysParam['Mfs'][1], sysParam['Mfs'][2], sysParam['Mfs'][3])


        fig = Figure(figsize=(5,4), dpi=200)
        ax = fig.add_subplot(111)
        
        for i in range(len(calcData)):
            newList = [p/deltaL for p in calcData[i]]
            ax.plot(Mfa, newList, lw=1)
        #plt.scatter(mindeltaPx, mindeltaPy, c='r', s=6)
        ax.set_ylabel("Pressure drop per meter (Pa/m)")
        ax.set_xlabel("Mass flow air (kg/s)")
        #newList = [x/1000 for x in pres]
        fig.legend(Mfs, loc=5, title='Mass flow \nsolids (kg/s)')
        plt.rcParams.update({'font.size': 12})
        return fig

if __name__ == '__main__':
    sim_gui = PressureCalc()
    # plotPressureDrop(pressureDropCalc())
    # plotOnlyPressure(pressureDropCalc())
    # plotFromCSV()
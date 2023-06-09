#!/usr/bin/env python
from scipy.optimize import minimize
import copy
import tkinter as tk
import tkinter.ttk as ttk

from Pressure import PressureCalc
import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter.messagebox import showinfo
from tkinter import Menu
from tkinter.filedialog import asksaveasfilename

font1 = ("Arial", 10)


class SimulationGui():
    sysParameters = {'Diam': ('Pipe Diameter (m)', 0.1),
                     'H_Length': ('Horizontal Length excl bends (m)', 8.156),
                     'V_Length': ('Vertical Length excl bends (m)', 17.541),
                     'nbH': ('Number of Horizontal Bends', 4),
                     'nbV': ('Number of Vertical Bends', 1),
                     'Br': ('Bend Ratio', 7.5),
                     'ksiB': ('Bend friction coefficient', 0.5),
                     'g': ('Gravitational acceleration (m/s2)', 9.81),
                     'R': ('Gas constant (kJ/(kg*K))', 0.287),
                     'T': ('Temperature (K)', 293),
                     'pres': ('Inlet pressure (Pa)', 400000),
                     'dynVisc': ('dynamic viscosity (Pa*s)', 18.1e-6),
                     'rho_p': ('Particle Density (kg/m3)', 2000),
                     'dp': ('Particle Diameter (um)', 1e-5),
                     'Mfs': ('Mass flow solids (kg/s)', 0.3, 5.3, 11),
                     "Mfa": ('Mass flow air (kg/s)', 0.05, 1.5, 10000)}
    entries = {}
    tempSysParam = {}
    pressureDrop = []

    def __init__(self):
        # Create gui window called root
        self.root = tk.Tk()
        self.root.wm_title("Pneumatic Conveying GUI")

        # self.root.geometry('800x600')
        self.lbl0 = tk.Label(self.root, text="Pneumatic Conveying Simulation - Transportation Pipeline",
                             font=("Arial", 15))
        self.lbl0.grid(row=0, columnspan=2, sticky=tk.W+tk.E)

        self.createHeader()
        self.createMenu()
        self.createLeftFrame()

        # Create a matplotlib canvas
        self.fig = Figure(figsize=(10, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(row=2, column=1,padx=30, sticky=tk.W)

        # Initialization run. Start an iterative loop and event loop
        self.afterId = self.root.after(150, self.gui_update)
        self.root.mainloop()
    
    def update(self):
        if len(self.pressureDrop) != 0:
            self.plotData(self.v.get())
        else:          
            showinfo(title='Incorrect', message='Start by computing, then press update to change the figure')
            return
        
    def save_fig(self):
        # shows dialog box and return the path to save fig
        path = asksaveasfilename(title='Select Folder',filetypes=[("PNG file", ".png")], defaultextension=".png", confirmoverwrite=True) 
        self.fig.savefig(fname=path,dpi=200)
        return

    def gui_update(self):
        # loop that runs iteratively after a given time in ms
        self.afterId = self.root.after(150, self.gui_update)

    def createHeader(self):
        self.headerFrame = tk.Frame(self.root)
        self.headerFrame.grid(row=1, column=1, padx=20, pady=20)

        self.v = tk.IntVar(self.headerFrame, 1)
        self.rb1 = tk.Radiobutton(self.headerFrame, text='Pressure drop',
                                  variable=self.v, value=1).grid(row=0, column=1, padx=10)
        self.rb2 = tk.Radiobutton(self.headerFrame, text='Power consumption',
                                  variable=self.v, value=2).grid(row=0, column=2, padx=10)
        self.rb3 = tk.Radiobutton(self.headerFrame, text='SLR', variable=self.v, value=3).grid(
            row=0, column=3, padx=10)

        self.updateButton = ttk.Button(
            self.headerFrame, text='Update', command=self.update).grid(row=0, column=4, padx=10)
        
        self.updateButton = ttk.Button(
            self.headerFrame, text='Save Figure', command=self.save_fig).grid(row=0, column=5, padx=20,sticky=tk.E)


    def createMenu(self):
        # Create Menu bar
        menubar = Menu(self.root, background='#ff8000', foreground='black', activebackground='white', activeforeground='black')  
        
        # File tab
        file = Menu(menubar, tearoff=0, foreground='black') 
        fileMenu = {'New':None, 'Open': None, 'Save': None, 'Save as': None} 
        for label,command in fileMenu.items():
            file.add_command(label=label,command=command)  

        file.add_separator()  
        file.add_command(label="Exit", command=self.root.quit)  
        menubar.add_cascade(label="File", menu=file)  

        # Edit tab
        edit = Menu(menubar, tearoff=0)
        editMenu = {'Undo':None, 'Cut': None, 'Copy': None, 'Paste': None, } 
        for label,command in editMenu.items():
            edit.add_command(label=label,command=command)    

        menubar.add_cascade(label="Edit", menu=edit)  

        # Help tab
        help = Menu(menubar, tearoff=0)  
        help.add_command(label="About", command=None)  
        menubar.add_cascade(label="Help", menu=help)  
            
        self.root.config(menu=menubar)


    def createLeftFrame(self):
        # Create the left frame
        self.leftFrame = tk.Frame(self.root)
        self.leftFrame.grid(row=2, column=0, padx=20, pady=20)

        counter = 2

        for param in self.sysParameters:
            entriesList = []
            if len(self.sysParameters[param]) == 4:
                self.massFlows = tk.Frame(self.leftFrame)
                self.massFlows.grid(row=counter, column=1)
                if param == 'Mfs':

                    l = ["Low", "High", "No. of Values"]
                    for i in range(len(l)):
                        tk.Label(self.massFlows, text=l[i]).grid(
                            row=counter, pady=2, column=i+1)

                tk.Label(self.leftFrame, text=self.sysParameters[param][0]).grid(
                    row=counter, sticky=tk.W+tk.S)
                counter += 1

                for i in range(1, len(self.sysParameters[param])):
                    value = tk.StringVar()
                    entry = tk.Entry(
                        self.massFlows, textvariable=value, width=8)
                    entry.grid(row=counter, column=i,
                               padx=10, pady=2, sticky='W')
                    entry.insert(0, self.sysParameters[param][i])
                    entriesList.append(entry)

            elif len(self.sysParameters[param]) == 2:
                tk.Label(self.leftFrame, text=self.sysParameters[param][0]).grid(
                    row=counter, pady=2, sticky='W')
                value = tk.StringVar()
                entry = tk.Entry(self.leftFrame, textvariable=value)
                entry.grid(row=counter, column=1, padx=10, sticky='W')
                entry.insert(0, self.sysParameters[param][1])
                entriesList.append(entry)

            self.entries[param] = entriesList
            counter += 1

        self.label = tk.Label(self.leftFrame, text='').grid(row=counter)

        # Calculate button
        counter += 1
        self.btn = tk.Button(self.leftFrame, text="Calculate",
                             bg='lightgray', fg='red',
                             command=self.calculateButton_callback)
        self.btn.grid(row=counter, column=0, sticky=tk.W+tk.E)


        # Clear plot button
        self.btn2 = tk.Button(self.leftFrame, text="Clear plot",
                              command=self.clearplot_callback)
        self.btn2.grid(row=counter, column=1, sticky=tk.W+tk.E)


    def calculateButton_callback(self):
        
        # Create a copy of the system parameters dictionary to modify
        self.tempSysParam = copy.deepcopy(self.sysParameters)

        # Update the system parameter values based on the entries
        for key in self.entries:
            tempList = []
            tempList.append(self.tempSysParam[key][0])
            for i in range(len(self.entries[key])):
                try:
                    tempList.append(float(self.entries[key][i].get()))
                except:
                    showinfo(title='Incorrect', message='Only integers and decimals are allowed!')
                    return
            self.tempSysParam[key] = tuple(tempList)

        self.tempSysParam['delta_L_bend'] = (
            'Length of 90deg bend(m)', 2*np.pi * self.tempSysParam['Br'][1] * self.tempSysParam['Diam'][1]/4)
        self.tempSysParam['deltaL'] = ('Total length (m)', (self.tempSysParam['H_Length'][1]) + (
            self.tempSysParam['V_Length'][1]) + self.tempSysParam['delta_L_bend'][1]*(self.tempSysParam['nbH'][1] + self.tempSysParam['nbV'][1]))

        self.tempSysParam['Mfa'] = (self.tempSysParam['Mfa'][0], np.linspace(
            self.tempSysParam['Mfa'][1], self.tempSysParam['Mfa'][2], int(self.tempSysParam['Mfa'][3])))
        self.tempSysParam['Mfs'] = (self.tempSysParam['Mfs'][0], np.linspace(
            self.tempSysParam['Mfs'][1], self.tempSysParam['Mfs'][2], int(self.tempSysParam['Mfs'][3])))

        self.pressureDrop = PressureCalc.pressureDropCalc(
            self, **self.tempSysParam)
        
        self.plotData(self.v.get())

    def plotData(self, val):
        # Executed when Pressure drop is clicked
        self.clearplot_callback()

        if val == 1:  # Plot pressure drop per meter
            for i in range(len(self.pressureDrop)):
                newList = [p/self.tempSysParam['deltaL'][1]
                           for p in self.pressureDrop[i]]
                self.ax.plot(self.tempSysParam['Mfa'][1], newList, lw=1)

            self.ax.set_ylabel("Pressure drop per meter (Pa/m)")
            self.ax.set_xlabel("Mass flow air (kg/s)")
            MfsLegend = [round(x, 2) for x in self.tempSysParam['Mfs'][1]]

        elif val == 2:  # Plot power consumption
            for i in range(len(self.pressureDrop)):

                power = [2*self.tempSysParam['Mfa'][1][j]*self.tempSysParam['R'][1]*self.tempSysParam['T'][1]*np.log(self.tempSysParam['pres'][1]/(
                    self.tempSysParam['pres'][1]-self.pressureDrop[i][j]))/self.tempSysParam['deltaL'][1] for j in range(len(self.pressureDrop[i]))]
                self.ax.plot(self.tempSysParam['Mfa'][1], power, lw=1)

            self.ax.set_ylabel("Power consumption per meter (kW/m)")
            self.ax.set_xlabel("Mass flow air (kg/s)")
            MfsLegend = [round(x, 2) for x in self.tempSysParam['Mfs'][1]]


        elif val == 3:
            for i in range(len(self.pressureDrop)):
                delta_P = [p/self.tempSysParam['deltaL'][1]
                           for p in self.pressureDrop[i]]
                SLR = [self.tempSysParam['Mfs'][1][i]/mfa for mfa in self.tempSysParam['Mfa'][1]]
                self.ax.plot(SLR, delta_P, lw=1)

            self.ax.set_ylabel("Pressure drop per meter (Pa/m)")
            self.ax.set_xlabel("SLR")
            MfsLegend = [round(x, 2) for x in self.tempSysParam['Mfs'][1]]
        
        self.ax.legend(MfsLegend, bbox_to_anchor=(0.995, 0.9), title='Mass flow \nsolids (kg/s)')
        self.canvas.draw()

    def powerConsumptionButton_callback(self):
        pass


    def clearplot_callback(self):
        # Clears the plot window
        self.ax.cla()
        self.canvas.draw()

    def quit_callback(self):
        # Quits the GUI Application
        self.root.after_cancel(self.afterId)
        self.root.destroy()


if __name__ == '__main__':
    sim_gui = SimulationGui()

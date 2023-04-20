import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import os
import sys
import re


class BandStructurePlotter:
    
    def __init__(self, **kwargs):
        

        self.fig_dir = kwargs.get('fig_dir', 'data_fig')
        self.input_filename = kwargs.get('input_filename', 'BandStructure.txt')
        self.output_filename = kwargs.get('output_filename', 'BandStructure')
        self.energy_range = kwargs.get('energy_range', [-2.0,2.0])
        self.max_value = kwargs.get('max_value', 1.0)
        self.min_value = kwargs.get('min_value', 0.0)
        self.need_sum = kwargs.get('need_sum', False)

        self.PlotBandStructure()


    def read_file(self,filename):
       
        
        
        # Open the file for reading
        with open(filename, 'r') as f:
            # Read the contents of the file
            contents = f.readlines()

            # Initialize an empty list to store the integers
            Spectrum_Data = []
            temp = []

            # Iterate over each line in the contents
            for line in contents:
                # Split the line into individual values
                values = line.split()
                
                # Convert the values to integers and add them to the list
                temp.append([float(v) for v in values])

            for i in range(len(temp[0])):
                Spectrum_Data.append(np.array(temp)[:,i])
            # Print the resulting list of integers
            print(len(Spectrum_Data),len(Spectrum_Data[0]))
        
        X_Data = np.arange(len(Spectrum_Data[0]))
        Y_Data = np.linspace(self.energy_range[0],self.energy_range[1],len(Spectrum_Data))

        Segment_labels = ['\u0393', 'K', 'M']
        Segments = [0,int(len(X_Data)/2),len(X_Data)]

        return X_Data, Y_Data, Spectrum_Data, Segments, Segment_labels


    def PlotBandStructure(self):
        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)
        
        spin_name = ["Up","Down"]
        
        for s in range(len(spin_name)):

            X_Data, Y_Data, Spectrum_Data, Segments, Segment_labels = self.read_file(spin_name[s]+self.input_filename)
            Route = X_Data
            Energy = Y_Data
            Energy = np.flip(Energy)
            band_diagram = np.flip(np.array(Spectrum_Data))
            for i in range(len(band_diagram)):
                band_diagram[i] = np.flip(band_diagram[i])

            # plot the data
            fig, ax = plt.subplots(figsize=(10,8))
            data_content = ax.imshow(band_diagram, cmap='Blues', vmax = self.max_value, vmin = self.min_value,aspect = len(Route)/len(Energy)*8/10.)
            
            # set y-label
            Energy_pos = [0,int(len(Energy)/4.),int(len(Energy)/2.),int(len(Energy)*3./4.),len(Energy)-1]
            ax.set_yticks(Energy_pos)
            ax.set_yticklabels(np.array(Energy)[Energy_pos])

            # set x-label
            route_pos = Segments
            ax.set_xticks(route_pos)
            ax.set_xticklabels(Segment_labels)
            
            
            # add a label to the point 
            ax.scatter(int(len(Route)/2.), int(len(Energy)/2.), marker='*', s=100, color='red',label='(K-point, Fermi-energy)')
            ax.legend()

            ax.set_xlabel('')
            ax.set_ylabel('Energy (eV)')
            cbar = fig.colorbar(data_content, label = 'band structure (1/eV)',ax=ax)

            fig.savefig(self.fig_dir + '/BandStructure_' + spin_name[s] + self.output_filename)

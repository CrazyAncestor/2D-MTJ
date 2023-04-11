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
        self.max_value = kwargs.get('max_value', 1.0)
        self.ext = kwargs.get('ext', '.png')

        self.PlotBandStructure()


    def read_file(self,filename):
        
        def find_xyn(string,pattern):

            # Find all matches of the pattern in the string
            matches = re.findall(pattern, string)
            if matches:
                if type(matches[0]) == tuple:
                    num = int(matches[0][0])
                else:
                    num = int(matches[0])
                return num
            else:
                return 0
        
        def appending_data(lines,line_id,tot_n):
            data = []
            for i in range(tot_n):
                line = lines[line_id+2+i]
                temp = [float(x.strip()) for x in line.split(',')]
                if len(temp)==1:
                    data.append(temp[0])
                else:
                    data.append(temp)
            return data

        #   Read the file
        with open(filename, 'r') as f:
            file_lines = f.readlines()
        
        #   Clean out nothing lines
        data_lines = []
        
        for i in range (len(file_lines)):
            line = file_lines[i]
            if line!='\n':
                data_lines.append(line)
        
        #   Extact spectrum data

        def extract_data(head,pattern):
            line_id = []
            row_num = []
            for i in range (len(data_lines)):
                line = data_lines[i]
                if line.startswith(head):
                    line_id.append(i)
                    row_num.append(find_xyn(str(line),pattern))
            data =[]
            if len(line_id)!=3:
                print('Data Error! Must contain all three spectra!')
                return
            for i in range(len(line_id)):
                data.append(appending_data(data_lines,line_id[i],row_num[i]))
            
            return data
        
        def extract_segment_x(x_list,val):

            # Find the index of the closest element in x_list to val
            idx = np.abs(np.array(x_list) - val).argmin()

            return idx
        
        X_Data = extract_data('xData [', r"\[(\d+)\]") 
        Y_Data = extract_data('yData [', r"\[(\d+)\]") 
        Spectrum_Data = extract_data('zData',r"\[(\d+) by (\d+)\]")

        G = 0.
        X = 0.3333
        Y = 1.0

        Segment_labels = []
        Segments = []
        for i in range(len(X_Data)):
            Segment_labels.append(['G', 'K', 'M'])
            Segments.append([extract_segment_x(X_Data[i],G), extract_segment_x(X_Data[i],X), extract_segment_x(X_Data[i],Y)])

        return X_Data, Y_Data, Spectrum_Data, Segments, Segment_labels


    def PlotBandStructure(self):
        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)
        
        X_Data, Y_Data, Spectrum_Data, Segments, Segment_labels = self.read_file(self.input_filename)
        spin_name = ["Sum","Up","Down"]
        
        for s in range(len(spin_name)):

            Route = X_Data[s]
            Energy = Y_Data[s]
            Energy = np.flip(Energy)
            band_diagram = np.flip(np.array(Spectrum_Data[s]))
            for i in range(len(band_diagram)):
                band_diagram[i] = np.flip(band_diagram[i])

            # plot the data
            fig, ax = plt.subplots(figsize=(6, 6*int(len(Energy)/len(Route))))
            data_content = ax.imshow(band_diagram, cmap='viridis', vmax = self.max_value)
            
            # set x-label
            route_pos = Segments[s]
            ax.set_xticks(route_pos)
            ax.set_xticklabels(Segment_labels[s])
            
            # set y-label
            Energy_pos = [0,int(len(Energy)/4.),int(len(Energy)/2.),int(len(Energy)*3./4.),len(Energy)-1]
            ax.set_yticks(Energy_pos)
            ax.set_yticklabels(np.array(Energy)[Energy_pos])

            # add a label to the point (3,6)
            ax.scatter(int(len(Route)/2.), int(len(Energy)/2.), marker='*', s=100, color='red',label='(K-point, Fermi-energy)')
            ax.legend()

            ax.set_xlabel('')
            ax.set_ylabel('Energy (eV)')
            cbar = fig.colorbar(data_content, label = 'band structure (1/eV)',ax=ax)
            fig.savefig(self.fig_dir + '/' + spin_name[s] + ' Spin Band Diagram' + self.ext)
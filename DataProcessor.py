import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import os
import sys


class DataProcessor:
    
    def __init__(self, **kwargs):
        #   Define physical constant
        self.g = 2.   #   Lande g-factor
        self.muB = 9.273e-24 #   Bohr magneton
        self.hbar = 1.054571817e-34  #   Reduced Planck Constant

        self.fig_dir = kwargs.get('fig_dir', 'data_fig')
        self.input_filenames = kwargs.get('input_filenames', [])
        self.original = kwargs.get('original', False)
        self.ext = kwargs.get('ext', '.png')

        self.plot_func = kwargs.get('plot_func', 'plot_MR')
        self.LockingRatio = kwargs.get('LockingRatio', 1.0e6)

        self.labels = kwargs.get('labels', ['forward','backward'])
        self.styles = kwargs.get('styles', ['dot','dot'])

        self.bool_out_of_plane = kwargs.get('bool_out_of_plane', True)
        self.BI_conver_ratio = kwargs.get('BI_conver_ratio', 200.)
        self.Hanle_Signal_Range = kwargs.get('Hanle_Signal_Range', 1000.)
        
        self.V_col = kwargs.get('V_col',0)
        self.I_col = kwargs.get('I_col',1)
        self.T_col = kwargs.get('T_col',2)
        self.R_col = kwargs.get('R_col',5)
        self.rs_col = kwargs.get('rs_col',0)
        self.int_col = kwargs.get('int_col',1)

        self.ProcessData()


    def read_data(self, filenames):

        def read_file(filename):
            with open(filename, encoding='utf-8')as file:
                # Skip the first line
                next(file)
                data = []

                # Read each line
                for line in file:

                    # Split the line into a list of strings
                    string_numbers = line.strip().split()
                    numbers=[]

                    # Convert each string to a number and store in a list
                    for i in range(len(string_numbers)):
                        num = 0.
                        try:                    
                            num = float(string_numbers[i])
                        except:                    
                            num = 0.
                        numbers.append(num)
                    data.append(numbers)
                    
                return np.array(data)

        result = []
        for i in range(len(filenames)):
            data = np.array([])
            for j in range(len(filenames[i])):
                temp = read_file(filenames[i][j])
                if j ==0:
                    data = temp
                else:
                    data = np.concatenate((data, temp), axis=0)
            result.append(data)
        return result
        

    #   Result plotting function
    def plot_figure(self, xs, ys, labels, x_label, y_label, title, figname, style):
        # Create a figure and axis object
        fig, ax = plt.subplots(figsize=(10, 8))

        # Set colors
        color_candidate = ['blue','red','black','green','brown','purple']
        color_num = len(color_candidate)

        # Plot the data on the axis object
        for i in range(len(xs)):
            if style[i] == 'dot':
                ax.scatter(xs[i], ys[i],label = labels[i],color=color_candidate[i%color_num],s=10)
            elif style[i] == 'line':
                ax.plot(xs[i], ys[i],label = labels[i],color=color_candidate[i%color_num])

        # Set the x and y axis labels
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        # Set the title of the plot
        ax.set_title(title)

        # Set legend
        ax.legend()

        # Show the plot
        if self.original:
            plt.savefig(self.fig_dir+'/original '+figname+self.ext)
        else:
            plt.savefig(self.fig_dir+'/'+figname+self.ext)

    #   Filter out outliers in data
    def eliminate_outliers(self, x, y, threshold=3.):
        def find_outliers(data,threshold=3.):
            z_scores = np.abs((data - np.mean(data)) / np.std(data))
            outliers = np.where(z_scores > threshold)[0]
            return outliers
        outlier_indices = find_outliers(y,threshold)
        x_new = np.delete(x, outlier_indices)
        y_new = np.delete(y, outlier_indices)
        return np.array(x_new),np.array(y_new)
    
    def sort_x_y(self, x, y):
        sort_indices = np.argsort(x)
        x_sorted = x[sort_indices]
        y_sorted = y[sort_indices]
        return np.array(x_sorted),np.array(y_sorted)

    def give_some_range_of_data(self, x, y, value_range, reverse=False):
        data_range = np.where(np.logical_and(x > value_range[0], x < value_range[1]))[0]
        if reverse:
            data_range = np.where(np.logical_or(x < value_range[0], x > value_range[1]))[0]
        x_new = []
        y_new = []
        
        for i in range(len(data_range)):
            x_new.append(x[data_range[i]])
            y_new.append(y[data_range[i]])

        x_new = np.array(x_new)
        y_new = np.array(y_new)
        return x_new, y_new
    
    #   Define the function to fit
    def inverse(self,x):
        return 1/x

    def log(self,x):
        return np.log(x)

    def identity(self,x):
        return x

    def linear(self,x,*para):
        a, b = para
        return a*x + b
        
    def parabolic_func(self,x, *para):
        a, b, c = para
        return a*x**2 + b*x + c
    
    def deviation_func(self,x, *para):
        a = para
        return a

    def Hanle_effect(self,Bz, *para):
        taus,r0 = para
        omega_L = (Bz /10000.)*self.g*self.muB/self.hbar # Turning Bz from Gauss to Tesla
        y2 = r0 /(1+(taus*omega_L)**2)
        return y2


    
    #   Data processing and figure making
    def plot_ACR_Monitor(self):
        unit_convert = [1.0,self.LockingRatio]
        self.x_col = self.T_col
        self.y_col = self.R_col

        def text_func(x):
            return f'Band Gap = {(x/11606*1e3):.2e} meV'
        
        self.plot_data(unit_convert,self.inverse,self.log,self.linear,[1.,1.],['Temp(K)','1/T(1/K)'],['R(Ohm)','log R'],['R-Temp curve','Boltzman Model fitting of R-Temp'],['R-Temp data','R-T fitting'],text_func=text_func)
    
    def plot_IV(self):
        unit_convert = [1.0,1.0]
        self.x_col = self.V_col
        self.y_col = self.I_col

        def text_func(x):
            return f'Resistance = {(1/x):.2e} Ohm'
        
        self.plot_data(unit_convert,self.identity,self.identity,self.linear,[1.,1.],['V(V)','V(V)'],['I(A)','I(A)'],['I-V data','I-V curve'],['I-V data','I-V fitting'],text_func=text_func)

    def plot_MR(self):
        self.x_col = self.I_col
        self.y_col = self.R_col
        unit_convert = [200.0,self.LockingRatio]

        def subtraction(x,y):
            z = []
            for i in range(len(x)):
                z.append(x[i] - y[i])
            return z

        if self.original:
            self.plot_data(unit_convert,self.identity,self.identity,self.parabolic_func,[1.,1.,1.],['B(B)','B(G)'],['R(Ohm)','R(Ohm)'],['B-MR curve','B-MR curve'],['B-MR data','B-MR fitting'])
        
        elif not self.bool_out_of_plane:
            self.plot_data(unit_convert,self.identity,self.identity,self.parabolic_func,[1.,1.,1.],['B(B)','B(G)'],['R(Ohm)','R(Ohm)'],['B-MR curve','B-MR curve'],['B-MR data','B-MR fitting'],data_point_style='line')

        elif self.bool_out_of_plane and not self.original:
            B, MR_raw, MR_para = self.plot_data(unit_convert,self.identity,self.identity,self.parabolic_func,[1.,1.,1.],['B(B)','B(G)'],['R(Ohm)','R(Ohm)'],['B-MR curve','B-MR curve'],['B-MR data','B-MR fitting'])
        
            #   Fit Hanle signal
            Hanle_data = subtraction(MR_raw , MR_para)
            unit_convert = [1.,1.]

            def text_func(x):
                return f'relaxation time = {x:.2e}sec'
            
            self.plot_data(unit_convert,self.identity,self.identity,self.Hanle_effect,[1e-10,1.],['B(B)','B(G)'],['R(Ohm)','R(Ohm)'],['Hanle Signal','Hanle Signal & Fitting'],['Hanle Signal','Hanle Signal & Fitting'],give_other_data=True,datax=B,datay=Hanle_data,text_func=text_func)

    def plot_Raman(self):
        self.x_col = self.rs_col
        self.y_col = self.int_col
        unit_convert = [1.0,1.0]
        
        self.original = True
        self.plot_data(unit_convert,self.identity,self.identity,self.identity,[1.],['raman shift(1/cm)'],['Intensity'],['Raman Spectrum'],['RamanSpectrum'],data_point_style='line')

    #   Data processing and figure making
    def plot_data(self,unit_convert,x_func,y_func,z_func,p0,x_labels,y_labels,titles,fignames,give_other_data=False,datax=[],datay=[],text_func=None,data_point_style='dot'):
        x_data_all = []
        y_data_all = []
        N = len(self.data)
        for i in range(N):
            x_data = self.data[i][:,self.x_col]*unit_convert[0]
            y_data = self.data[i][:,self.y_col]*unit_convert[1]
            x_data_all.append(x_data)
            y_data_all.append(y_data)
        
        if give_other_data:
            x_data_all = datax
            y_data_all = datay

        # Plot original data
        if self.original:
            self.plot_figure(x_data_all,y_data_all,self.labels,x_labels[0],y_labels[0],title=titles[0],figname=fignames[0],style=[data_point_style]*N)
            return

        x_new = []
        y_new = []
        y_fit = []
        fit_labels = []
        
        
        for i in range(N):
            # Eliminate outliers
            x_data_all[i] , y_data_all[i] = self.sort_x_y(x_data_all[i],y_data_all[i])
            x_data_all[i] , y_data_all[i] = self.eliminate_outliers(x_data_all[i],y_data_all[i])

            # Fit the ACR function of T-R curve
            x_new0 = x_func(x_data_all[i])
            y_new0 = y_func(y_data_all[i])

            if self.plot_func == 'plot_MR' and z_func == self.parabolic_func and np.max(np.abs(x_new0))>=2000:
                xf, yf = self.give_some_range_of_data(x_new0,y_new0,value_range=[-7500,7500])
                xf, yf = self.give_some_range_of_data(xf,yf,value_range=[-1000,1000],reverse=True)
                
                popt, pcov = curve_fit(lambda x, *para: z_func(x, *para),xf,yf,p0=p0)
                y_fit0 = z_func(x_new0,*popt)

            elif self.plot_func == 'plot_MR' and z_func == self.parabolic_func and np.max(np.abs(x_new0))<2000:
                y_fit0 = np.min(y_new0)*np.ones(len(y_new0))
                

            elif self.plot_func == 'plot_MR' and z_func == self.Hanle_effect:
                x_new0,y_new0 = self.give_some_range_of_data(x_new0,y_new0,value_range=[-1000,1000])
                popt, pcov = curve_fit(lambda x, *para: z_func(x, *para),x_new0,y_new0,p0=p0)
                y_fit0 = z_func(x_new0,*popt)
            
            elif self.plot_func == 'plot_ACR_Monitor':
                x_new0,y_new0 = self.give_some_range_of_data(x_new0,y_new0,value_range=[1/150.,1/20.])
                popt, pcov = curve_fit(lambda x, *para: z_func(x, *para),x_new0,y_new0,p0=p0)
                y_fit0 = z_func(x_new0,*popt)

            else:
                popt, pcov = curve_fit(lambda x, *para: z_func(x, *para),x_new0,y_new0,p0=p0)
                y_fit0 = z_func(x_new0,*popt)

            y_fit.append(y_fit0)
            x_new.append(x_new0)
            y_new.append(y_new0)
            
            # Text message telling band gap energy
                        
            try:
                fit_labels.append(self.labels[i]+' fitting, '+text_func(np.abs(popt[0])))
            except:
                fit_labels.append('')
                

        # Plot figures
        self.plot_figure(x_new+x_new,y_new+y_fit,self.labels+fit_labels,x_labels[1],y_labels[1],title=titles[1],figname=fignames[1],style=[data_point_style]*N+['line']*N)
        self.plot_figure(x_data_all,y_data_all,self.labels,x_labels[0],y_labels[0],title=titles[0],figname=fignames[0],style=[data_point_style]*N)

        return x_new,y_new,y_fit

    def ProcessData(self):
        #   Preliminary setup
        if self.plot_func == 'plot_MR':
            self.R_col = 5
        elif self.plot_func == 'plot_ACR_Monitor':
            self.R_col = 6
        
        if self.original:
            self.fig_dir = 'Original Data' + self.fig_dir

        if not os.path.exists(self.fig_dir):
            os.makedirs(self.fig_dir)
        else:
            print("Directory already exists")   

        #   Read data
        self.data = self.read_data(self.input_filenames)
        
        #   Plot figure
        if self.plot_func == 'plot_MR':
            self.plot_MR()
        elif self.plot_func == 'plot_ACR_Monitor':
            self.plot_ACR_Monitor()
        elif self.plot_func == 'plot_IV':
            self.plot_IV()
        elif self.plot_func == 'plot_Raman':
            self.plot_Raman()
        else:
            print('Plotting function error!')
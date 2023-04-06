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

        self.bool_out_of_plane = kwargs.get('bool_out_of_plane', True)
        self.BI_conver_ratio = kwargs.get('BI_conver_ratio', 200.)
        self.Hanle_Signal_Range = kwargs.get('Hanle_Signal_Range', 1000.)
        
        self.V_col = kwargs.get('V_col',0)
        self.I_col = kwargs.get('I_col',1)
        self.T_col = kwargs.get('T_col',2)
        self.R_col = kwargs.get('R_col',5)

        self.ProcessData()


    def read_data(self, filenames):
        def read_file(filename):
            with open(filename, encoding='utf-8')as file:

                try:
                    # Create a CSV reader object
                    reader = csv.reader(file)
                    for i in range(1):
                        next(reader)
                    data = [[float(num) for num in row] for row in reader]
                    data = np.array(data)
                    return data

                except:
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

        result = np.array([])

        for i in range(len(filenames)):
            data = read_file(filenames[i])
            if i ==0:
                result = data
            else:
                result = np.concatenate((result, data), axis=0)
        return result

    #   Result plotting function
    def plot_figure(self, xs, ys, labels, x_label, y_label, title, plot_filename, style,text_message=''):
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

        # Print text messages
        ax.text(0.95,0.05, text_message, fontsize=12, ha='right', va='bottom', transform=ax.transAxes, bbox={'facecolor': 'white', 'pad': 0.5, 'edgecolor': 'black', 'alpha': 0.5})


        # Show the plot
        if self.original:
            plt.savefig(self.fig_dir+'/original '+plot_filename+self.ext)
        else:
            plt.savefig(self.fig_dir+'/'+plot_filename+self.ext)

    #   Filter out outliers in data
    def eliminate_outliers(self, x, y, threshold=3.):
        def find_outliers(data,threshold=3.):
            z_scores = np.abs((data - np.mean(data)) / np.std(data))
            outliers = np.where(z_scores > threshold)[0]
            return outliers
        outlier_indices = find_outliers(y,threshold)
        x_new = np.delete(x, outlier_indices)
        y_new = np.delete(y, outlier_indices)
        return x_new,y_new
    
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
    
    def fit_hanle_signal(self, B_field, R_data, left_or_right='both'):

        #   Define the function to fit
        def parabolic_func(x, a, b, c):
            return a*x**2 + b*x + c

        def Hanle_effect(Bz,r0,taus):
            omega_L = (Bz /10000.)*self.g*self.muB/self.hbar # Turning Bz from Gauss to Tesla
            y2 = r0 /(1+(taus*omega_L)**2)
            return y2

        if np.max(B_field)<2000.:
            Bz = B_field
            R_parabolic_background = np.ones(len(R_data))*np.min(R_data)
            R_Hanle_signal = R_data -R_parabolic_background
            
            # Fit the value of relax time in Hanle signals with direct model
            if left_or_right=='both':
                Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[-self.Hanle_Signal_Range,self.Hanle_Signal_Range])
            elif left_or_right=='left':
                Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[-self.Hanle_Signal_Range,0.])
            elif left_or_right=='right':
                Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[0.,self.Hanle_Signal_Range])

            r0i = np.max(R_Hanle)
            tausi = 1/(self.g*self.muB/self.hbar)/np.std(Bz/10000.)
            p0 = [r0i,tausi]
            popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle,p0=p0)
            r0, taus = popt2
            R_fit_directModel = Hanle_effect(Bz,r0,taus)
            taus_message = f'relaxation time = {taus:.2e}sec'

            return Bz, R_parabolic_background, R_Hanle, R_fit_directModel, taus_message

        # Filter out the parabolic background
        SR1 = 7500.
        SR2 = 1000.

        B_para, R_para = self.give_some_range_of_data(B_field,R_data,[-SR1,SR1])
        B_para, R_para = self.give_some_range_of_data(B_para, R_para,[-SR2,SR2],reverse=True)
        popt, pcov = curve_fit(parabolic_func, B_para, R_para)

        a, b, c = popt
        R_parabolic_background = parabolic_func(B_field,a,b,c)
        R_Hanle_signal = R_data - R_parabolic_background

        # Fit the value of relax time in Hanle signals with direct model
        if left_or_right=='both':
            Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[-self.Hanle_Signal_Range,self.Hanle_Signal_Range])
        elif left_or_right=='left':
            Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[-self.Hanle_Signal_Range,0.])
        elif left_or_right=='right':
            Bz, R_Hanle = self.give_some_range_of_data(B_field,R_Hanle_signal,[0.,self.Hanle_Signal_Range])

        r0i = np.max(R_Hanle)
        tausi = 1/(self.g*self.muB/self.hbar)/np.std(Bz/10000.)
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle,p0=p0)
        r0, taus = popt2
        R_fit_directModel = Hanle_effect(Bz,r0,taus)
        taus_message = f'relaxation time = {taus:.2e}sec'

        return Bz, R_parabolic_background, R_Hanle, R_fit_directModel, taus_message
    
    #   Data processing and figure making
    def plot_ACR_Monitor(self):
        T_data = self.data[:,self.T_col]
        R_data = self.data[:,self.R_col]*self.LockingRatio
        T_data , R_data = self.give_some_range_of_data(T_data,R_data,value_range)
        # Plot original data
        if self.original:
            self.plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title='R-Temp curve',plot_filename='R-Tempcurve',style=['dot'])
            return
        
        # Eliminate outliers
        T_data , R_data = self.sort_x_y(T_data,R_data)
        T_data , R_data = self.eliminate_outliers(T_data,R_data)

        def ACR_func(T_conv,a,b):
            return a*T_conv + b
        
        # Fit the ACR function of T-R curve
        popt, pcov = curve_fit(ACR_func,1/T_data,np.log(R_data))
        a, b = popt
        LOGRVR_fit = ACR_func(1/T_data,a,b)
        
        # Text message telling band gap energy
        text_message = f'Band Gap = {(a/11606*1e3):.2e} meV'

        # Plot figures
        self.plot_figure([1/T_data,1/T_data],[np.log(R_data),LOGRVR_fit],['R-T data','R-T fitting'],'1/T(1/K)','lnR(Ohm)',title='R-Temp curve',plot_filename="R_T_fitting",style=['dot','line'],text_message=text_message)
        self.plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title='R-Temp curve',plot_filename='R-Tempcurve',style=['dot'],text_message=text_message)

    def plot_IV(self):
        I_data = self.data[:,self.I_col]
        V_data = self.data[:,self.V_col]

        # Plot original data
        if self.original:
            self.plot_figure([V_data],[I_data],['I-V data'],'V(V)','I(A)',title='I-V curve',plot_filename='I-Vcurve',style=['dot'])
            return

        # Eliminate outliers
        V_data , I_data = self.eliminate_outliers(V_data,I_data)
        V_data , I_data = self.sort_x_y(V_data,I_data)

        def linear_func(x,a,b):
            return a*x + b
        
        # Fit the linear function of I-V curve
        popt, pcov = curve_fit(linear_func,V_data,I_data)
        a, b = popt
        I_fit = linear_func(V_data,a,b)
        
        # Text message telling resistance
        text_message = f'Resistance = {(1/a):.2e} Ohm'

        # Plot figures
        self.plot_figure([V_data,V_data],[I_data,I_fit],['I-V data','I-V fitting'],'V(V)','I(A)',title='I-V curve',plot_filename='I-Vcurve',style=['dot','line'],text_message=text_message)

    def plot_MR(self):
        I_data = self.data[:,self.I_col]
        R_data = self.data[:,self.R_col]

        # Plot original data
        if self.original:
            self.plot_figure([I_data],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title='B-MR curve',plot_filename='B-MRcurve',style=['dot'])
            return
        
        # Eliminate outliers
        I_data, R_data = self.eliminate_outliers(I_data,R_data,threshold = 3.)
        I_data, R_data = self.eliminate_outliers(I_data,R_data,threshold = 3.)
        I_data, R_data = self.sort_x_y(I_data,R_data)

        # Unit conversion (B_field: Gauss, R_data: Ohm)
        B_field = I_data * self.BI_conver_ratio
        R_data = R_data * self.LockingRatio

        if self.bool_out_of_plane:
            #   Fit the data
            choose_plot_range = ['both','left','right']
            for i in range(len(choose_plot_range)):
                Bz, R_parabolic_background, R_Hanle, R_fit, taus_message = self.fit_hanle_signal(B_field,R_data,left_or_right=choose_plot_range[i])
                text_message = taus_message
                if i ==0:
                    self.plot_figure([B_field,B_field],[R_data,R_parabolic_background],['Original Data','Parabolic Background'],'Bz(Gauss)','Resistance(Omega)',title='B-MR curve',plot_filename='B-MRcurve',style=['dot','line'],text_message=text_message)
                self.plot_figure([Bz,Bz],[R_Hanle,R_fit],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)',title='Hanle signal vs. Hanle fitting',plot_filename='HanleSignalFitting_'+choose_plot_range[i],style=['dot','line'],text_message=text_message)

        else:
            self.plot_figure([B_field],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title='B-MR curve',plot_filename='B-MRcurve',style=['dot'])

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
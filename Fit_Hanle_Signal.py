import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit
from scipy.integrate import quad
import os

#   Define physical constant
g = 2.   #   Lande g-factor
muB = 9.273e-24 #   Bohr magneton
hbar = 1.054571817e-34  #   Reduced Planck Constant

#   Data preprocessing function
def sort_x_y(x, y):
    sort_indices = np.argsort(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    return np.array(x_sorted),np.array(y_sorted)

def eliminate_outliers(x, y, threshold=3.):
    def find_outliers(data,threshold=3.):
        z_scores = np.abs((data - np.mean(data)) / np.std(data))
        outliers = np.where(z_scores > threshold)[0]
        return outliers
    outlier_indices = find_outliers(y,threshold)
    x_new = np.delete(x, outlier_indices)
    y_new = np.delete(y, outlier_indices)
    return np.array(x_new),np.array(y_new)

#   Data reading function
def read_data(filename,data_start_row):
    # Open the file in read mode
    with open(filename, 'r',encoding="utf-8-sig") as file:

        # Create a CSV reader object
        reader = csv.reader(file)
        for i in range(data_start_row):
            next(reader)
        data = [[float(num) for num in row] for row in reader]
        data = np.array(data)
        return data

#   Data fitting function
def fit_hanle_signal(data,fig_dir,ext,Bz_colnum=2,R_colnum=3):

    #   Define the function to fit
    def parabolic_func(x, a, b, c):
        return a*x**2 + b*x + c

    def Hanle_effect(Bz,r0,taus):
        omega_L = (Bz/10000.)*g*muB/hbar # Turning Bz from Gauss to Tesla
        y2 = r0 /(1+(taus*omega_L)**2)
        return y2

    def give_some_range_of_data(x, y, value_range, reverse=False):
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
    
    # Give the magnetic field and signal
    Bz = data[:,Bz_colnum]
    R_OriData = data[:,R_colnum]
    Signal_range = 2000

    # Sort X,Y
    Bz, R_OriData = sort_x_y(Bz, R_OriData)

    # Filter outliers
    Bz, R_OriData = eliminate_outliers(Bz, R_OriData)

    # Filter out the parabolic background
    if np.max(Bz)>2000:
        x1,y1 = give_some_range_of_data(Bz,R_OriData,[-1000,1000],reverse=True)
        x1,y1 = give_some_range_of_data(x1,y1,[-7500,7500])
        popt, pcov = curve_fit(parabolic_func, x1, y1)

        a, b, c = popt
        R_parabolic_background = parabolic_func(Bz,a,b,c)
        R_Hanle_signal = R_OriData - R_parabolic_background

        # Fit the value of relax time in Hanle signals with Lorentzian model
        Bz_hanle,R_Hanle_signal = give_some_range_of_data(Bz,R_Hanle_signal,[-Signal_range,Signal_range])
        r0i = 1.
        tausi = 1e-10
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz_hanle, R_Hanle_signal,p0=p0)
        r0, taus = popt2
        taus = np.abs(taus)
        Bz_fit = np.linspace(Bz_hanle[0],Bz_hanle[-1],100)
        R_fit_LorentzianModel = Hanle_effect(Bz_fit,r0,taus)
        
        # Plot the results
        plot_figure([Bz,Bz],[R_OriData*1e-1,R_parabolic_background*1e-1],['Raw Data','Parabolic Background Signal'],'B(G)',r'$V(\mu V)$',fig_dir+' Out-of-plane MR',fig_dir+'/'+fig_dir+'OriData_ParabolicBackground'+ext,-1)
        plot_figure([Bz_hanle,Bz_fit],[R_Hanle_signal*1e-1,R_fit_LorentzianModel*1e-1],['Hanle Signal','Hanle Fitting'],'B(G)',r'$\Delta V(\mu V)$','Hanle signal & Fitting',fig_dir+'/'+fig_dir+'HanleSignalFitting'+ext,taus)

    else:
        R_background = np.min(R_OriData)*np.ones(len(R_OriData))
        R_Hanle_signal = R_OriData - R_background

        # Fit the value of relax time in Hanle signals with Lorentzian model
        Bz_hanle,R_Hanle_signal = give_some_range_of_data(Bz,R_Hanle_signal,[-Signal_range,Signal_range])
        r0i = 1.
        tausi = 1e-10
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz_hanle, R_Hanle_signal,p0=p0)
        r0, taus = popt2
        taus = np.abs(taus)
        Bz_fit = np.linspace(Bz_hanle[0],Bz_hanle[-1],100)
        R_fit_LorentzianModel = Hanle_effect(Bz_fit,r0,taus)

        # Plot the results
        plot_figure([Bz,Bz],[R_OriData*1e-1,R_background*1e-1],['Raw Data','Parabolic Background Signal'],'B(G)',r'$V(\mu V)$',fig_dir+' Out-of-plane MR',fig_dir+'/'+fig_dir+'OriData_Background'+ext,-1)
        plot_figure([Bz_hanle,Bz_fit],[R_Hanle_signal*1e-1,R_fit_LorentzianModel*1e-1],['Hanle Signal','Hanle Fitting'],'B(G)',r'$\Delta V(\mu V)$','Hanle signal & Fitting',fig_dir+'/'+fig_dir+'HanleSignalFitting'+ext,taus)



#   Result plotting function
def plot_figure(xs,ys,labels,x_label,y_label,title,plot_filename,taus,inplane=False):
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the data on the axis object
    for i in range(len(ys)):
        if i ==0:
            ax.scatter(xs[i], ys[i],label = labels[i],color='blue')
        else:
            ax.plot(xs[i], ys[i],label = labels[i],color='red')

    # Set the x and y axis labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


    # Set the title of the plot
    ax.set_title(title)

    # Set legend
    ax.legend()

    # Tell the relaxation time
    if taus!=-1:
        taus_message = f'relaxation time = {taus:.2e}sec'
        ax.text(0.02,0.98, taus_message, fontsize=12, ha='left', va='top', transform=ax.transAxes, bbox={'facecolor': 'white', 'pad': 0.5, 'edgecolor': 'black', 'alpha': 0.5})

    # Show the plot
    plt.savefig(plot_filename)


def DataProcessing(datafile,fig_dir,xcol,ycol,ext='.png',inplane=False):
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    #   Give the data filename
    data_start_row = 1
    Bz_colnum = xcol
    R_colnum = ycol

    #   Read the data
    data = read_data(filename=datafile,data_start_row=data_start_row)

    if inplane:
        # Give the magnetic field and signal
        Bz = data[:,Bz_colnum]
        R_OriData = data[:,R_colnum]

        # Sort X,Y
        Bz, R_OriData = sort_x_y(Bz, R_OriData)

        # Filter outliers
        Bz, R_OriData = eliminate_outliers(Bz, R_OriData,threshold=2.)
        Bz, R_OriData = eliminate_outliers(Bz, R_OriData,threshold=2.)

        R_parallel = R_OriData[3]
        R_ratio=(R_OriData - R_parallel*np.ones(len(R_OriData)) )/R_parallel 

        def Hanle_effect(Bz,r0,taus):
            omega_L = (Bz/10000.)*g*muB/hbar # Turning Bz from Gauss to Tesla
            y2 = r0 /(1+(taus*omega_L)**2)
            return y2

        r0i = 1.
        tausi = 1e-10
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_ratio,p0=p0)
        r0, taus = popt2
        taus = np.abs(taus)
        Bz_fit = np.linspace(Bz[0],Bz[-1],100)
        R_fit_LorentzianModel = Hanle_effect(Bz_fit,r0,taus)

        # Plot the results
        plot_figure([Bz,Bz_fit],[R_ratio*100,R_fit_LorentzianModel*100],['Raw Data',''],'B(G)','MR Ratio (%)','In-plane MR',fig_dir+'/'+fig_dir+'InplaneMR'+ext,-1,inplane=inplane)
        
    else:
        #   Fit the data
        fit_hanle_signal(data,fig_dir,ext=ext,Bz_colnum=Bz_colnum,R_colnum=R_colnum)

#   Main function

if __name__ == '__main__':
    extension = '.svg'
    DataProcessing('TC_Gr_MR_data.csv','SLGr',0,1,ext=extension)
    DataProcessing('TC_Gr_MR_data.csv','MLGr',2,3,ext=extension)
    DataProcessing('TC_Gr_MR_data.csv','SLFET',4,5,ext=extension)
    DataProcessing('TC_Gr_MR_data.csv','WHGr',6,7,ext=extension)

    DataProcessing('TCDATA_Inplane_MLGR.csv','Inplane_MLGr',0,1,ext=extension,inplane=True)

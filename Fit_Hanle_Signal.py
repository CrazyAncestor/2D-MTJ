import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit
from scipy.integrate import quad

#   Define physical constant
g = 2.   #   Lande g-factor
muB = 9.273e-24 #   Bohr magneton
hbar = 1.054571817e-34  #   Reduced Planck Constant

#   Data reading function
def read_data(filename):
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
def fit_hanle_signal(data,Bz_colnum=2,R_colnum=3):

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
    # Filter out the parabolic background
    if np.max(Bz)>2000:
        x1,y1 = give_some_range_of_data(Bz,R_OriData,[-2500,2500],reverse=True)
        popt, pcov = curve_fit(parabolic_func, x1, y1)

        a, b, c = popt
        R_parabolic_background = parabolic_func(Bz,a,b,c)
        R_Hanle_signal = R_OriData - R_parabolic_background

        # Fit the value of relax time in Hanle signals with Lorentzian model
        Bz_hanle,R_Hanle_signal = give_some_range_of_data(Bz,R_Hanle_signal,[-Signal_range,Signal_range])
        r0i = 1.
        tausi = 1/(g*muB/hbar)
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz_hanle, R_Hanle_signal,p0=p0)
        r0, taus = popt2
        taus = np.abs(taus)
        R_fit_LorentzianModel = Hanle_effect(Bz_hanle,r0,taus)
        
        # Plot the results
        plot_figure(Bz,[R_OriData,R_parabolic_background],['Original Data','Parabolic Background'],'Bz(Gauss)','Resistance(Omega)','Ori_Data vs. Para Background','OriData_ParabolicBackground.png',-1)
        plot_figure(Bz_hanle,[R_Hanle_signal,R_fit_LorentzianModel],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)','Hanle signal vs. Hanle fitting','HanleSignalFitting.png',taus)

    else:
        R_background = np.min(R_OriData)*np.ones(len(R_OriData))
        R_Hanle_signal = R_OriData - R_background

        # Fit the value of relax time in Hanle signals with Lorentzian model
        Bz_hanle,R_Hanle_signal = give_some_range_of_data(Bz,R_Hanle_signal,[-Signal_range,Signal_range])
        r0i = 1.
        tausi = 1/(g*muB/hbar)
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz_hanle, R_Hanle_signal,p0=p0)
        r0, taus = popt2
        taus = np.abs(taus)
        R_fit_LorentzianModel = Hanle_effect(Bz_hanle,r0,taus)
        
        # Plot the results
        plot_figure(Bz,[R_OriData,R_background],['Original Data','Background'],'Bz(Gauss)','Resistance(Omega)','Ori_Data vs. Background','OriData_Background.png',-1)
        plot_figure(Bz_hanle,[R_Hanle_signal,R_fit_LorentzianModel],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)','Hanle signal vs. Hanle fitting','HanleSignalFitting.png',taus)



#   Result plotting function
def plot_figure(x,ys,labels,x_label,y_label,title,plot_filename,taus):
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the data on the axis object
    for i in range(len(ys)):
        if i ==0:
            ax.scatter(x, ys[i],label = labels[i])
        else:
            ax.plot(x, ys[i],label = labels[i])

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


#   Main function

if __name__ == '__main__':

    #   Give the data filename
    filename = 'data.csv'
    data_start_row = 1
    Bz_colnum = 2
    R_colnum = 3

    #   Read the data
    data = read_data(filename=filename)

    #   Fit the data
    fit_hanle_signal(data,Bz_colnum=Bz_colnum,R_colnum=R_colnum)

    

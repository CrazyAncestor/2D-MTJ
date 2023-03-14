import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit
from scipy.integrate import quad

#   Define physical constant

g = 2.   #   Lande g-factor
muB = 9.273e-24 #   Bohr magneton
hbar = 1.054571817e-34  #   Reduced Planck Constant

#   Define device parameter
D = 1.5
delta_x = 50e-6

#   Give the data filename
filename = 'data.csv'

#   Data reading function
def read_data(filename):
    # Open the file in read mode
    with open(filename, 'r') as file:

        # Create a CSV reader object
        reader = csv.reader(file)
        next(reader)
        data = [[float(num) for num in row] for row in reader]
        data = np.array(data)
        return data

#   Data fitting function
def fit_hanle_signal(data,physical_model='direct'):

    #   Define the function to fit
    def parabolic_func(x, a, b, c):
        return a*x**2 + b*x + c

    def Hanle_effect(Bz,r0,taus):
        omega_L = (Bz/10000.)*g*muB/hbar # Turning Bz from Gauss to Tesla
        y2 = r0 /(1+(taus*omega_L)**2)
        return y2

    def Hanle_diffusion_model(Bz,r0,taus):
        result = []
        for i in range(len(Bz)):
            #   Give the diffusion spin-injection model
            def f_diffusion(t,dx,Bz,D):
                x = 1/(4*np.pi*D*t)**0.5 *np.exp(-(dx**2)/(4*np.pi*D*t)) * np.cos(g*muB*Bz*t/hbar) *np.exp(-t/taus)
                return x
            
            # define the integration limits
            ti = 0.
            tf = 10.*taus
            params = (delta_x,Bz[i]/1e4,D)
            answer, error = quad(f_diffusion, ti, tf, args=params)
            result.append(answer)
        
        result = result/np.max(result)*r0
        return result

    # Give the magnetic field and signal
    Bz = data[:,2]
    R_OriData = data[:,3]

    # Filter out the parabolic background
    popt, pcov = curve_fit(parabolic_func, Bz, R_OriData)

    a, b, c = popt
    R_parabolic_background = parabolic_func(Bz,a,b,c)
    R_Hanle_signal = R_OriData - R_parabolic_background

    # Fit the value of relax time in Hanle signals with direct model
    r0i = 1.
    tausi = 1/(g*muB/hbar)
    p0 = [r0i,tausi]
    popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle_signal,p0=p0)
    r0, taus = popt2
    R_fit_directModel = Hanle_effect(Bz,r0,taus)
    
    if physical_model=='direct':
        return Bz, R_OriData, R_parabolic_background, R_Hanle_signal, R_fit_directModel, taus

    # Fit the value of relax time in Hanle signals with diffusion model
    r0i = r0
    tausi = taus*1.0
    p0 = [r0i,tausi]
    popt3, pcov3 = curve_fit(Hanle_diffusion_model, Bz, R_Hanle_signal,p0=p0, bounds=([r0*0.7, 0.05*taus], [r0*1.5, 2.0*taus]))
    r0, taus = popt3
    R_fit_diffusionModel = Hanle_diffusion_model(Bz,r0,taus)
    
    if physical_model=='diffusion model':
        return Bz, R_OriData, R_parabolic_background, R_Hanle_signal, R_fit_diffusionModel, taus

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
    taus_message = f'relaxation time = {taus:.2e}sec'
    ax.text(0.02,0.98, taus_message, fontsize=12, ha='left', va='top', transform=ax.transAxes, bbox={'facecolor': 'white', 'pad': 0.5, 'edgecolor': 'black', 'alpha': 0.5})

    # Show the plot
    plt.savefig(plot_filename)


#   Main function

if __name__ == '__main__':
    #   Read the data
    data = read_data(filename=filename)

    #   Fit the data
    Bz, R_OriData, R_parabolic_background, R_Hanle_signal, R_fit, taus = fit_hanle_signal(data,physical_model='diffusion model')

    # Plot the results
    plot_figure(Bz,[R_OriData,R_parabolic_background],['Original Data','Parabolic Background'],'Bz(Gauss)','Resistance(Omega)','Ori_Data vs. Para Background','OriData_ParabolicBackground.png',taus)
    plot_figure(Bz,[R_Hanle_signal,R_fit],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)','Hanle signal vs. Hanle fitting','HanleSignalFitting.png',taus)



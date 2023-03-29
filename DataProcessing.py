import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad
import os


#   Define physical constant
g = 2.   #   Lande g-factor
muB = 9.273e-24 #   Bohr magneton
hbar = 1.054571817e-34  #   Reduced Planck Constant

def read_data(filenames):
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
def plot_figure(xs,ys,labels,x_label,y_label,title,plot_filename,style,text_message=' ',dot_size=10):
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(10, 8))

    # Set colors
    color_candidate = ['blue','red','black','green','brown','purple']
    color_num = len(color_candidate)

    # Plot the data on the axis object
    for i in range(len(xs)):
        if style[i] == 'dot':
            ax.scatter(xs[i], ys[i],label = labels[i],color=color_candidate[i%color_num],s=dot_size)
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
    plt.savefig(plot_filename+'.png')

#   Filter out outliers in data
def eliminate_outliers(x,y,threshold=3.):
    def find_outliers(data,threshold=3.):
        z_scores = np.abs((data - np.mean(data)) / np.std(data))
        outliers = np.where(z_scores > threshold)[0]
        return outliers
    outlier_indices = find_outliers(y,threshold)
    x_new = np.delete(x, outlier_indices)
    y_new = np.delete(y, outlier_indices)
    return x_new,y_new
def sort_x_y(x,y):
    sort_indices = np.argsort(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    return np.array(x_sorted),np.array(y_sorted)

def give_some_range_of_data(x,y,value_range,reverse=False):
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

def fit_hanle_signal(B_field, R_data, Hanle_Signal_Range = 1000., left_or_right='both'):

    #   Define the function to fit
    def parabolic_func(x, a, b, c):
        return a*x**2 + b*x + c

    def Hanle_effect(Bz,r0,taus):
        omega_L = (Bz /10000.)*g*muB/hbar # Turning Bz from Gauss to Tesla
        y2 = r0 /(1+(taus*omega_L)**2)
        return y2

    if np.max(B_field)<2000.:
        Bz = B_field
        R_parabolic_background = np.ones(len(R_data))*np.min(R_data)
        R_Hanle_signal = R_data -R_parabolic_background
        
        # Fit the value of relax time in Hanle signals with direct model
        if left_or_right=='both':
            Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Hanle_Signal_Range,Hanle_Signal_Range])
        elif left_or_right=='left':
            Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Hanle_Signal_Range,0.])
        elif left_or_right=='right':
            Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[0.,Hanle_Signal_Range])

        r0i = np.max(R_Hanle)
        tausi = 1/(g*muB/hbar)/np.std(Bz/10000.)
        p0 = [r0i,tausi]
        popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle,p0=p0)
        r0, taus = popt2
        R_fit_directModel = Hanle_effect(Bz,r0,taus)
        taus_message = f'relaxation time = {taus:.2e}sec'

        return Bz, R_parabolic_background, R_Hanle, R_fit_directModel, taus_message

    # Filter out the parabolic background
    SR1 = 7500.
    SR2 = 1000.

    B_para, R_para = give_some_range_of_data(B_field,R_data,[-SR1,SR1])
    B_para, R_para = give_some_range_of_data(B_para, R_para,[-SR2,SR2],reverse=True)
    popt, pcov = curve_fit(parabolic_func, B_para, R_para)

    a, b, c = popt
    R_parabolic_background = parabolic_func(B_field,a,b,c)
    R_Hanle_signal = R_data - R_parabolic_background

    # Fit the value of relax time in Hanle signals with direct model
    if left_or_right=='both':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Hanle_Signal_Range,Hanle_Signal_Range])
    elif left_or_right=='left':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Hanle_Signal_Range,0.])
    elif left_or_right=='right':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[0.,Hanle_Signal_Range])

    r0i = np.max(R_Hanle)
    tausi = 1/(g*muB/hbar)/np.std(Bz/10000.)
    p0 = [r0i,tausi]
    popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle,p0=p0)
    r0, taus = popt2
    R_fit_directModel = Hanle_effect(Bz,r0,taus)
    taus_message = f'relaxation time = {taus:.2e}sec'

    return Bz, R_parabolic_background, R_Hanle, R_fit_directModel, taus_message

#   Data processing and figure making
def plot_ACR_Monitor(data,fig_dir,original=False,T_col=2,R_col=6):
    T_data = data[:,T_col]
    R_data = data[:,R_col]

    # Plot original data
    if original:
        plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title='R-Temp curve',plot_filename=fig_dir+'/original '+'R-Tempcurve',style=['dot'],text_message='',dot_size=10)
        return
    
    # Eliminate outliers
    T_data , R_data = eliminate_outliers(T_data,R_data)
    T_data , R_data = sort_x_y(T_data,R_data)

    def ACR_func(T_conv,a,b):
        return a*T_conv + b
    
    # Fit the ACR function of T-R curve
    popt, pcov = curve_fit(ACR_func,1/T_data,np.log(R_data))
    a, b = popt
    LOGRVR_fit = ACR_func(1/T_data,a,b)
    
    # Text message telling band gap energy
    text_message = f'Band Gap = {(a/11606*1e3):.2e} meV'

    # Plot figures
    plot_figure([1/T_data,1/T_data],[np.log(R_data),LOGRVR_fit],['R-T data','R-T fitting'],'1/T(1/K)','lnR(Ohm)',title='R-Temp curve',plot_filename=fig_dir+'/'+"R_T_fitting",style=['dot','line'],text_message=text_message,dot_size=10)
    plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title='R-Temp curve',plot_filename=fig_dir+'/'+'R-Tempcurve',style=['dot'],text_message=text_message,dot_size=10)

def plot_IV(data,fig_dir,original=False,I_col=1,V_col=0):
    I_data = data[:,I_col]
    V_data = data[:,V_col]

    # Plot original data
    if original:
        plot_figure([V_data],[I_data],['I-V data'],'V(V)','I(A)',title='I-V curve',plot_filename=fig_dir+'/original '+'I-Vcurve',style=['dot'],text_message='',dot_size=10)
        return

    # Eliminate outliers
    V_data , I_data = eliminate_outliers(V_data,I_data)
    V_data , I_data = sort_x_y(V_data,I_data)

    def linear_func(x,a,b):
        return a*x + b
    
    # Fit the linear function of I-V curve
    popt, pcov = curve_fit(linear_func,V_data,I_data)
    a, b = popt
    I_fit = linear_func(V_data,a,b)
    
    # Text message telling resistance
    text_message = f'Resistance = {(1/a):.2e} Ohm'

    # Plot figures
    plot_figure([V_data,V_data],[I_data,I_fit],['I-V data','I-V fitting'],'V(V)','I(A)',title='I-V curve',plot_filename=fig_dir+'/'+'I-Vcurve',style=['dot','line'],text_message=text_message,dot_size=10)

def plot_MR(data,fig_dir,original=False,bool_out_of_plane=True,BI_conver_ratio = 200.,LockingRatio=1.0e6,Hanle_Signal_Range=1000.,I_col=1,R_col=5):
    I_data = data[:,I_col]
    R_data = data[:,R_col]

    # Plot original data
    if original:
        plot_figure([I_data],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title='B-MR curve',plot_filename=fig_dir+'/original'+'B-MRcurve',style=['dot'],text_message='',dot_size=10)
        return
    
    # Eliminate outliers
    I_data, R_data = eliminate_outliers(I_data,R_data,threshold = 3.)
    I_data, R_data = eliminate_outliers(I_data,R_data,threshold = 3.)
    I_data, R_data = sort_x_y(I_data,R_data)

    # Unit conversion (B_field: Gauss, R_data: Ohm)
    B_field = I_data * BI_conver_ratio
    R_data = R_data * LockingRatio

    if bool_out_of_plane:
        #   Fit the data
        choose_plot_range = ['both','left','right']
        for i in range(len(choose_plot_range)):
            Bz, R_parabolic_background, R_Hanle, R_fit, taus_message = fit_hanle_signal(B_field,R_data,Hanle_Signal_Range=Hanle_Signal_Range,left_or_right=choose_plot_range[i])
            text_message = taus_message
            if i ==0:
                plot_figure([B_field,B_field],[R_data,R_parabolic_background],['Original Data','Parabolic Background'],'Bz(Gauss)','Resistance(Omega)',title='B-MR curve',plot_filename=fig_dir+'/'+'B-MRcurve',style=['dot','line'],text_message=text_message,dot_size=10)
            plot_figure([Bz,Bz],[R_Hanle,R_fit],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)',title='Hanle signal vs. Hanle fitting',plot_filename=fig_dir+'/'+'HanleSignalFitting_'+choose_plot_range[i],style=['dot','line'],text_message=text_message,dot_size=10)

    else:
        plot_figure([B_field],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title='B-MR curve',plot_filename=fig_dir+'/'+'B-MRcurve',style=['dot'],text_message='',dot_size=10)

def data_processing(**kwargs):
    
    fig_dir = kwargs.get('fig_dir', 'data_fig')
    input_filenames = kwargs.get('input_filenames', [])
    plot_func = kwargs.get('plot_func', plot_MR)
    original = kwargs.get('original', False)

    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    else:
        print("Directory already exists")   

    #   Read data
    data = read_data(input_filenames)
    
    #   Plot figure
    if plot_func is plot_MR:
        bool_out_of_plane = kwargs.get('bool_out_of_plane', True)
        BI_conver_ratio = kwargs.get('BI_conver_ratio', 200.)
        LockingRatio = kwargs.get('LockingRatio', 1.0e6)
        Hanle_Signal_Range = kwargs.get('Hanle_Signal_Range', 1000.)
        I_col = kwargs.get('I_col',1)
        R_col = kwargs.get('R_col',5)
        plot_func(data,fig_dir,original,BI_conver_ratio=BI_conver_ratio,LockingRatio=LockingRatio,Hanle_Signal_Range=Hanle_Signal_Range,I_col=I_col,R_col=R_col)
    else:
        plot_func(data,fig_dir,original=original)

   

#   Main function

if __name__ == '__main__':

    #   Plot Ting-Chun's data
    TC_file = './ori_data/TC_Gr_MR_data.csv'
    data_processing(fig_dir='TC Data Single layer',input_filenames=[TC_file],plot_func=plot_MR,BI_conver_ratio=1.0,LockingRatio=1.0,Hanle_Signal_Range=2500.,I_col=0,R_col=1)
    data_processing(fig_dir='TC Data Multilayer',input_filenames=[TC_file],plot_func=plot_MR,BI_conver_ratio=1.0,LockingRatio=1.0,Hanle_Signal_Range=2500.,I_col=2,R_col=3)
    data_processing(fig_dir='TC Data 1L FET Gr',input_filenames=[TC_file],plot_func=plot_MR,BI_conver_ratio=1.0,LockingRatio=1.0,Hanle_Signal_Range=2500.,I_col=4,R_col=5)
    
    #   Plot R-T
    RTfilename_total = './ori_data/ACR_monitor_total.txt'
    data_processing(fig_dir='ACR Monitor Total',input_filenames=[RTfilename_total],plot_func=plot_ACR_Monitor)
    data_processing(fig_dir='Original Data ACR Monitor Total',input_filenames=[RTfilename_total],plot_func=plot_ACR_Monitor,original=True)

    RTfilename_part = './ori_data/ACR_monitor_part.txt'
    data_processing(fig_dir='ACR Monitor part',input_filenames=[RTfilename_part],plot_func=plot_ACR_Monitor)
    data_processing(fig_dir='Original Data ACR Monitor part',input_filenames=[RTfilename_part],plot_func=plot_ACR_Monitor,original=True)

    #   1T out-of-plane MR
    Outplane_MR_1T_filename1 = './ori_data/AC_Scan_KA_1T_0T_T10K_1MOhm_wt3s_tc1s.txt'
    Outplane_MR_1T_filename2 = './ori_data/AC_Scan_KA_m1T_0T_T10K_1MOhm_wt3s_tc1s.txt'
    data_processing(fig_dir='1T Out-of-plane MR',input_filenames=[Outplane_MR_1T_filename1,Outplane_MR_1T_filename2],plot_func=plot_MR)
    data_processing(fig_dir='Original Data 1T Out-of-plane MR',input_filenames=[Outplane_MR_1T_filename1,Outplane_MR_1T_filename2],plot_func=plot_MR,original=True)
    
    #   1000G out-of-plane MR
    Outplane_MR_1000G_filename1 = './ori_data/AC_Scan_KA_1000G_0G_T10K_1MOhm_wt3s_tc1s.txt'
    Outplane_MR_1000G_filename2 = './ori_data/AC_Scan_KA_m1000G_0G_T10K_1MOhm_wt3s_tc1s.txt'
    data_processing(fig_dir='1000G Out-of-plane MR',input_filenames=[Outplane_MR_1000G_filename1,Outplane_MR_1000G_filename2],plot_func=plot_MR)
    data_processing(fig_dir='Original Data 1000G Out-of-plane MR',input_filenames=[Outplane_MR_1000G_filename1,Outplane_MR_1000G_filename2],plot_func=plot_MR,original=True)
    
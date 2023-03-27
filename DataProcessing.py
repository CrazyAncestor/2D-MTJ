import numpy as np
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
    plt.savefig(plot_filename)

#   Filter out outliers in data
def find_outliers(data,threshold=3.):
    z_scores = np.abs((data - np.mean(data)) / np.std(data))
    outliers = np.where(z_scores > threshold)[0]
    return outliers

def eliminate_outliers(x,y,threshold=3.):
    outlier_indices = find_outliers(y,threshold)
    x_new = np.delete(x, outlier_indices)
    y_new = np.delete(y, outlier_indices)
    return x_new,y_new

def plot_outliers(data,x_col,y_col,title='Test plotting outliers',plot_filename='TEST_PLOT_OUTLIERS.png',dot_size=10):
    x = data[:,x_col]
    y = data[:,y_col]
    outlier_indices = find_outliers(y)
    x_out = []
    y_out = []
    for i in range (len(outlier_indices)):
        x_out.append(x[outlier_indices[i]])
        y_out.append(y[outlier_indices[i]])
    
    # filter out outliers
    x_new,y_new = eliminate_outliers(x,y)

    

    plot_figure([x_new],[y_new],['filtered_data'],'AU','AU',title='plot filtered data',plot_filename='filtered_data.png',style=['dot'],text_message=text_message,dot_size=dot_size)
    plot_figure([x,x_out],[y,y_out],['data','outliers'],'AU','AU',title=title,plot_filename=plot_filename,style=['dot','dot'],text_message=text_message,dot_size=dot_size)

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

def fit_hanle_signal(B_field, R_data, left_or_right='both'):

    #   Define the function to fit
    def parabolic_func(x, a, b, c):
        return a*x**2 + b*x + c

    def Hanle_effect(Bz,r0,taus):
        omega_L = (Bz /10000.)*g*muB/hbar # Turning Bz from Gauss to Tesla
        y2 = r0 /(1+(taus*omega_L)**2)
        return y2

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
    Signal_Range = 1000.
    if left_or_right=='both':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Signal_Range,Signal_Range])
    elif left_or_right=='left':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[-Signal_Range,0.])
    elif left_or_right=='right':
        Bz, R_Hanle = give_some_range_of_data(B_field,R_Hanle_signal,[0.,Signal_Range])

    r0i = np.max(R_Hanle)
    tausi = 1/(g*muB/hbar)/np.std(Bz/10000.)
    p0 = [r0i,tausi]
    popt2, pcov2 = curve_fit(Hanle_effect, Bz, R_Hanle,p0=p0)
    r0, taus = popt2
    R_fit_directModel = Hanle_effect(Bz,r0,taus)
    taus_message = f'relaxation time = {taus:.2e}sec'

    return Bz, R_parabolic_background, R_Hanle, R_fit_directModel, taus_message

#   Data processing and figure making
def plot_ARC_Monitor(data,T_col=2,R_col=6,title='R-Temp curve',plot_filename='R-Tempcurve.png',text_message=' ',dot_size=10):
    T_data = data[:,T_col]
    R_data = data[:,R_col]
    T_data , R_data = eliminate_outliers(T_data,R_data)

    def ARC_func(T_conv,a,b):
        return a*T_conv + b
    
    # Fit the ARC function of T-R curve
    popt, pcov = curve_fit(ARC_func,1/T_data,np.log(R_data))
    a, b = popt
    LOGRVR_fit = ARC_func(1/T_data,a,b)
    
    text_message = f'Band Gap = {(a/11606*1e3):.2e} meV'
    plot_figure([1/T_data,1/T_data],[np.log(R_data),LOGRVR_fit],['R-T data','R-T fitting'],'1/T(1/K)','lnR(Ohm)',title=title,plot_filename="R_T_fitting.png",style=['dot','line'],text_message=text_message,dot_size=dot_size)

    plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)

def plot_IV(data,I_col=1,V_col=0,title='I-V curve',plot_filename='I-Vcurve.png',dot_size=10):
    I_data = data[:,I_col]
    V_data = data[:,V_col]
    V_data , I_data = eliminate_outliers(V_data,I_data)

    def linear_func(x,a,b):
        return a*x + b
    
    # Fit the linear function of I-V curve
    popt, pcov = curve_fit(linear_func,V_data,I_data)
    a, b = popt
    I_fit = linear_func(V_data,a,b)
    
    text_message = f'Resistance = {(1/a):.2e} Ohm'
    plot_figure([V_data,V_data],[I_data,I_fit],['I-V data','I-V fitting'],'V(V)','I(A)',title=title,plot_filename=plot_filename,style=['dot','line'],text_message=text_message,dot_size=dot_size)

def plot_MR(data,bool_out_of_plane=False,LockingRatio=1.0e6,I_col=1,R_col=5,title='B-MR curve',plot_filename='B-MRcurve.png',style=['dot','line'],text_message=' ',dot_size=10):
    I_data = data[:,I_col]
    R_data = data[:,R_col]

    I_data, R_data = eliminate_outliers(I_data,R_data,threshold = 3.)
    I_data, R_data = eliminate_outliers(I_data,R_data,threshold = 3.)

    # Unit conversion (B_field: Gauss, R_data: Ohm)
    B_field = I_data *200.
    R_data = R_data *LockingRatio

    if bool_out_of_plane:
        #   Fit the data
        choose_plot_range = ['both','left','right']
        for i in range(len(choose_plot_range)):
            Bz, R_parabolic_background, R_Hanle, R_fit, taus_message = fit_hanle_signal(B_field,R_data,left_or_right=choose_plot_range[i])
            text_message = taus_message
            if i ==0:
                plot_figure([B_field,B_field],[R_data,R_parabolic_background],['Original Data','Parabolic Background'],'Bz(Gauss)','Resistance(Omega)',title=title,plot_filename=plot_filename,style=style,text_message=text_message,dot_size=dot_size)
            plot_figure([Bz,Bz],[R_Hanle,R_fit],['Hanle Signal','Hanle Fitting'],'Bz(Gauss)','Resistance(Omega)',title='Hanle signal vs. Hanle fitting',plot_filename='HanleSignalFitting_'+choose_plot_range[i]+'.png',style=style,text_message=text_message,dot_size=dot_size)

    else:
        plot_figure([B_field],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)


#   Main function

if __name__ == '__main__':

    figure_dir = 'fitting_data_figure'
    
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    else:
        print("Directory already exists")

    RTemp_filename = './ori_data/ACR_monitor_KA.txt'
    IV_filename = './ori_data/IV_KA_1Vtom1V_6.6K_wt100ms.txt'
    
    Inplane_MR_400G_filename1 = '.\ori_data\AC_Scan_KA_m400GTO0G_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_400G_filename2 = '.\ori_data\AC_Scan_KA_0GTO400G_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_400G_filename3 = '.\ori_data\AC_Scan_KA_400GTO0G_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_400G_filename4 = '.\ori_data\AC_Scan_KA_0GTOm400G_T10K_10MOhm_wt3s_tc1s.txt'

    Inplane_MR_1T_filename1 = '.\ori_data\AC_Scan_KA_m1TO0T_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_1T_filename2 = '.\ori_data\AC_Scan_KA_0TO1T_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_1T_filename3 = '.\ori_data\AC_Scan_KA_1TO0T_T10K_10MOhm_wt3s_tc1s.txt'
    Inplane_MR_1T_filename4 = '.\ori_data\AC_Scan_KA_0TOm1T_T10K_10MOhm_wt3s_tc1s.txt'

    Outplane_MR_filename = './ori_data/OutOfPlane_AC_Scan_KA_m1TTO1T_T10K_10MOhm_wt3s_tc1s.txt'

    
    
    #   R vs. Temp
    RTemp_data =read_data([RTemp_filename])
    plot_ARC_Monitor(RTemp_data,plot_filename= figure_dir +'/'+ 'R-Temp Curve.png')
    
    #   I-V curve
    IV_data = read_data([IV_filename])
    plot_IV(IV_data,plot_filename= figure_dir +'/'+ 'I-V Curve.png')

    #   In-plane MR 1T
    Inplane_MR_1T_data1 = read_data([Inplane_MR_1T_filename1,Inplane_MR_1T_filename2])
    plot_MR(Inplane_MR_1T_data1,title='In-plane MR 1T forward',plot_filename= figure_dir +'/'+ 'In-plane MR data 1T forward.png')

    Inplane_MR_1T_data2 = read_data([Inplane_MR_1T_filename3,Inplane_MR_1T_filename4])
    plot_MR(Inplane_MR_1T_data2,title='In-plane MR 1T backward',plot_filename= figure_dir +'/'+ 'In-plane MR data 1T  backward.png')

    #   In-plane MR 400G
    Inplane_MR_400G_data1 = read_data([Inplane_MR_400G_filename1,Inplane_MR_400G_filename2])
    plot_MR(Inplane_MR_400G_data1,title='In-plane MR 400G forward',plot_filename= figure_dir +'/'+ 'In-plane MR data 400G forward.png')

    Inplane_MR_400G_data2 = read_data([Inplane_MR_400G_filename3,Inplane_MR_400G_filename4])
    plot_MR(Inplane_MR_400G_data2,title='In-plane MR 400G backward',plot_filename= figure_dir +'/'+ 'In-plane MR data 400G  backward.png')

    #   Out-of-plane MR
    Outplane_MR_data = read_data([Outplane_MR_filename])
    plot_MR(Outplane_MR_data,bool_out_of_plane=True,title='Out-of-plane MR',plot_filename= figure_dir +'/'+ 'Out-of-plane MR data.png')
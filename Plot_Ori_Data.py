import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

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

#   Plot original data
def plot_ARC_Monitor(data,T_col=2,R_col=6,title='R-Temp curve',plot_filename='R-Tempcurve.png',text_message=' ',dot_size=10):
    T_data = data[:,T_col]
    R_data = data[:,R_col]

    plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)

def plot_IV(data,I_col=1,V_col=0,title='I-V curve',plot_filename='I-Vcurve.png',dot_size=10):
    I_data = data[:,I_col]
    V_data = data[:,V_col]
    
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

    B_field = I_data *200.
    R_data = R_data *LockingRatio

    plot_figure([B_field],[R_data],['B-MR data'],'B(Gauss)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)


#   Main function

if __name__ == '__main__':

    figure_dir = 'ori_data_figure'
    
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)
    else:
        print("Directory already exists")

    RTemp_filename = './ori_data/ACR_monitor_KA.txt'
    IV_filename = './ori_data/IV_(2023_03_21_21_24)NO.2_KA.txt'
    
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
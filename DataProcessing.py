import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def read_data(filename):
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

#   Result plotting function
def plot_figure(xs,ys,labels,x_label,y_label,title,plot_filename,style,text_message='HelloWorld',dot_size=10):
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

    text_message = ''

    plot_figure([x_new],[y_new],['filtered_data'],'AU','AU',title='plot filtered data',plot_filename='filtered_data.png',style=['dot'],text_message=text_message,dot_size=dot_size)
    plot_figure([x,x_out],[y,y_out],['data','outliers'],'AU','AU',title=title,plot_filename=plot_filename,style=['dot','dot'],text_message=text_message,dot_size=dot_size)

def plot_ARC_Monitor(data,T_col=2,R_col=6,title='R-Temp curve',plot_filename='R-Tempcurve.png',dot_size=10):
    T_data = data[:,T_col]
    R_data = data[:,R_col]
    T_data , R_data = eliminate_outliers(T_data,R_data)

    text_message = ''
    plot_figure([T_data],[R_data],['R-Temp data'],'Temp(K)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)

#   Data processing and figure making
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

def plot_MR(data,LockingRatio=1.0e6,I_col=1,R_col=5,title='B-MR curve',plot_filename='B-MRcurve.png',dot_size=10):
    I_data = data[:,I_col]
    R_data = data[:,R_col]

    I_data, R_data = eliminate_outliers(I_data,R_data)

    # Unit conversion (B_field: Gauss, R_true: Ohm)
    B_field = I_data *200.
    R_true = R_data *LockingRatio
    
    text_message = ''
    plot_figure([B_field],[R_true],['B-MR data'],'B(Gauss)','R(Ohm)',title=title,plot_filename=plot_filename,style=['dot'],text_message=text_message,dot_size=dot_size)


#   Main function

if __name__ == '__main__':

    RTemp_filename = 'ACR_monitor_KC_2.txt'
    IV_filename = 'IV_KC_6.6K.txt'
    Inplane_MR_filename = 'Inplane_AC_Scan_KC_10K_1T_1sec.txt'
    
    RTemp_data =read_data(filename=RTemp_filename)
    plot_ARC_Monitor(RTemp_data)
    
    IV_data = read_data(filename=IV_filename)
    plot_IV(IV_data)

    Inplane_MR_data = read_data(Inplane_MR_filename)
    plot_MR(Inplane_MR_data)


   
import numpy as np
import csv
import matplotlib.pyplot as plt
import sys
import os


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
        
        print(data)
        return np.array(data)

#   Result plotting function
def plot_figure(xs,ys,labels,x_label,y_label,title,plot_filename):
    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(10, 8))

    # Set colors
    color_candidate = ['blue','blue','red','red','green','green']
    color_num = len(color_candidate)

    # Plot the data on the axis object
    for i in range(len(xs)):
        ax.scatter(xs[i], ys[i],label = labels[i],color=color_candidate[(2*i)%color_num],s=10)
        ax.plot(xs[i], ys[i],color=color_candidate[(2*i)%color_num])

    # Set the x and y axis labels and scales
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_yscale('log')

    # Set the title of the plot
    ax.set_title(title)

    # Set legend
    ax.legend()

    # Show the plot
    plt.savefig(plot_filename)


def sort_x_y(x,y):
    sort_indices = np.argsort(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    return np.array(x_sorted),np.array(y_sorted)


#   Data processing and figure making

def plot_conductivity(data,labels,x_col=0,y_col=2,plot_filename='conductivity'):
    xs = []
    ys = []

    for i in range(len(data)):
        x_data = data[i][:,x_col]
        y_data = data[i][:,y_col]
        xs.append(x_data)
        ys.append(y_data)


    # Plot figures
    plot_figure(xs,ys,labels,'layer number','conductivity (e^2/h)',title='Conductivity',plot_filename=plot_filename)

#   Main function

if __name__ == '__main__':

    #   Data Reading
    data_minor = read_data('minor.txt')
    data_major = read_data('major.txt')
    data_anti = read_data('anti.txt')
    data = [data_minor,data_major,data_anti]
    labels=['minor channel in parallel config.','major channel in parallel config.','anti-parallel config.']

    #   Plot figure
    plot_filename = sys.argv[1]
    plot_conductivity(data,labels,x_col=0,y_col=2,plot_filename=plot_filename)

    
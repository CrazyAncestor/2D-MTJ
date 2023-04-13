import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_trans(filename,data_dir,plot_dir,ext='.png',value_max=1.0):
    
    spin_name = ["Up","Down"]
    
    for s in range(len(spin_name)):
        data = []
        with open(data_dir+spin_name[s]+filename+'.txt', "r") as file:
            # Read the file into a list of strings
            lines = file.readlines()
            # Exclude the first three lines
            lines = lines[3:]
            
            # Read the remaining lines line by line
            for line in lines:
                # Split the line into a list of numbers
                numbers = line.split()
                # Convert the numbers from strings to floats (or ints, if appropriate)
                numbers = [float(num) for num in numbers]
                # Do something with the numbers
                data.append(numbers)

        data = np.array(data)
        ka = data[:,0]
        kb = data[:,1]
        transmission = data[:,2]

        # plot the data
        fig, ax = plt.subplots()
        data_content = ax.scatter(ka, kb, c=transmission, cmap='Reds', vmin = 0., vmax = value_max)
        ax.set_xlabel('ka')
        ax.set_ylabel('kb')
        cbar = fig.colorbar(data_content, label = 'transmission',ax=ax)
        fig.savefig(plot_dir+spin_name[s]+filename+ext)

if __name__ == '__main__':
    data_dir = './'
    plot_dir = 'Transmission_Spectra/'
    structure = sys.argv[1]

    #   Plot the transmission results
    trans_max = sys.argv[2]
    ext = sys.argv[3]
    
    filename1 = "transmission_Para"+structure
    filename2 = "transmission_Anti"+structure

    plot_trans(filename1,data_dir,plot_dir,ext,trans_max)
    plot_trans(filename2,data_dir,plot_dir,ext,trans_max)
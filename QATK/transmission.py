import numpy as np
import matplotlib.pyplot as plt
def write_transmissionfile(filename):
    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    transmission_spectrum = nlread(filename+'.hdf5', TransmissionSpectrum)[0]
    # Print out all k-dependent transmission coefficients
    spin_name = ["Up","Down"]
    spins = [Spin.Up,Spin.Down] 

    #print data
    for s in range(len(spin_name)):
        with open(spin_name[s]+'transmission_'+filename+'.txt', "w") as file:

            data = transmission_spectrum.transmission(spin = spins[s])
            energies = transmission_spectrum.energies()
            kpoints = transmission_spectrum.kpoints()

            file.write("Transmission coefficient of spin component " + spin_name[s] +"\n")
            for i in range(data.shape[0]):
                file.write('Transmission at energy = %12.6f eV \n' % (energies[i].inUnitsOf(eV)))
                file.write('      kx         ky      transmission \n')
                for j in range(data.shape[1]):
                    file.write('%10.4f %10.4f %16.6e\n' % (kpoints[j][0],kpoints[j][1],data[i][j]))

def plot_trans(filename):
    
    spin_name = ["Up","Down"]
    
    for s in range(len(spin_name)):
        data = []
        with open(spin_name[s]+filename+'.txt', "r") as file:
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
        print(data.shape)
        ka = data[:,0]
        kb = data[:,1]
        transmission = data[:,2]

        fig,ax = plt.subplots()
        data_content = ax.scatter(ka, kb, c=transmission, cmap='Reds')
        ax.set_xlabel('ka')
        ax.set_ylabel('kb')
        fig.colorbar(data_content, label = 'transmission')
        fig.savefig(spin_name[s]+filename+'.png')

parafile = 'ParaCobaltGrCobalt'
#antifile = 'AntiCobaltGrCobalt'

write_transmissionfile(parafile)
#write_transmissionfile(antifile)

filename1 = "transmission_ParaCobaltGrCobalt"
#filename2 = "transmission_AntiCobaltGrCobalt"

plot_trans(filename1)
#plot_trans(filename2)

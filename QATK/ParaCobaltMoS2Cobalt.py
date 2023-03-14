# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Two-probe Configuration
# -------------------------------------------------------------

# -------------------------------------------------------------
# Left Electrode
# -------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np

# Set up lattice
vector_a = [5.48, 0.0, 0.0]*Angstrom
vector_b = [0.0, 3.165, 0.0]*Angstrom
vector_c = [0.0, 0.0, 6.12]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
left_electrode_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt,  Cobalt]

# Define coordinates
left_electrode_coordinates = [[0.45675, 2.37375, 1.02], [3.19725, 0.79125, 1.02], [2.284, 2.37375, 3.06], [5.0245, 0.79125, 3.06], [1.37025, 0.79125, 5.1], [4.11125, 2.37375, 5.1]]*Angstrom

# Set up configuration
left_electrode = BulkConfiguration(
    bravais_lattice=left_electrode_lattice,
    elements=left_electrode_elements,
    cartesian_coordinates=left_electrode_coordinates
    )

# Add tags
left_electrode.addTags('Selection 0')

# -------------------------------------------------------------
# Right Electrode
# -------------------------------------------------------------

# Set up lattice
vector_a = [5.48, 0.0, 0.0]*Angstrom
vector_b = [0.0, 3.165, 0.0]*Angstrom
vector_c = [0.0, 0.0, 6.12]*Angstrom
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
right_electrode_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt,  Cobalt]

# Define coordinates
right_electrode_coordinates = [[1.37025, 0.79125, 1.02], [4.11125, 2.37375, 1.02], [2.284, 2.37375, 3.06], [5.0245, 0.79125, 3.06], [0.45675, 2.37375, 5.1], [3.19725, 0.79125, 5.1]]*Angstrom

# Set up configuration
right_electrode = BulkConfiguration(
    bravais_lattice=right_electrode_lattice,
    elements=right_electrode_elements,
    cartesian_coordinates=right_electrode_coordinates
    )

# Add tags
right_electrode.addTags('Selection 0')

# -------------------------------------------------------------
# Central Region
# -------------------------------------------------------------

# Set up lattice
vector_a = [5.48, 0.0, 0.0]*Angstrom
vector_b = [0.0, 3.165, 0.0]*Angstrom
vector_c = [0.0, 0.0, 17.439999999999998]*Angstrom
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
central_region_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, 
                           Sulfur, Sulfur, Molybdenum,  Sulfur, Sulfur, Molybdenum,
                           Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, Cobalt]

# Define coordinates
central_region_coordinates = [[0.45675, 2.37375, 1.02], [3.19725, 0.79125, 1.02], [2.284, 2.37375, 3.06], [5.0245, 0.79125, 3.06], [1.37025, 0.79125, 5.1], [4.11125, 2.37375, 5.1], [0.45675, 2.37375, 7.139999999999999], [0.45675, 2.37375, 10.299999999999999], [1.37025, 0.79125, 8.719999999999999], [3.19725, 0.79125, 7.139999999999999], [3.19725, 0.79125, 10.299999999999999], [4.11125, 2.37375, 8.719999999999999], [1.37025, 0.79125, 12.34], [4.11125, 2.37375, 12.34], [2.284, 2.37375, 14.38], [5.0245, 0.79125, 14.38], [0.45675, 2.37375, 16.42], [3.19725, 0.79125, 16.42]]*Angstrom

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

# Add tags
central_region.addTags('Selection 0', [0, 1, 16, 17])

device_configuration = DeviceConfiguration(
    central_region,
    [left_electrode, right_electrode]
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------
#----------------------------------------
# Basis Set
#----------------------------------------
basis_set = [
    GGABasis.Molybdenum_DoubleZetaPolarized,GGABasis.Sulfur_DoubleZetaPolarized,
    GGABasis.Cobalt_DoubleZetaPolarized,
    ]

#----------------------------------------
# Exchange-Correlation
#----------------------------------------
exchange_correlation = SGGA.PBE

#----------------------------------------
# Numerical Accuracy Settings
#----------------------------------------
left_electrode_k_point_sampling = MonkhorstPackGrid(
    na=7,
    nb=7,
    nc=100,
    )
left_electrode_numerical_accuracy_parameters = NumericalAccuracyParameters(
    electron_temperature=1200.0*Kelvin,
    k_point_sampling=left_electrode_k_point_sampling,
    )

right_electrode_k_point_sampling = MonkhorstPackGrid(
    na=7,
    nb=7,
    nc=100,
    )
right_electrode_numerical_accuracy_parameters = NumericalAccuracyParameters(
    electron_temperature=1200.0*Kelvin,
    k_point_sampling=right_electrode_k_point_sampling,
    )

device_k_point_sampling = MonkhorstPackGrid(
    na=7,
    nb=7,
    nc=100,
    )
device_numerical_accuracy_parameters = NumericalAccuracyParameters(
    electron_temperature=1200.0*Kelvin,
    k_point_sampling=device_k_point_sampling,
    )
#----------------------------------------
# Poisson Solver Settings
#----------------------------------------
left_electrode_poisson_solver = FastFourier2DSolver(
    boundary_conditions=[[PeriodicBoundaryCondition(),PeriodicBoundaryCondition()],
                         [PeriodicBoundaryCondition(),PeriodicBoundaryCondition()],
                         [PeriodicBoundaryCondition(),PeriodicBoundaryCondition()]]
    )

right_electrode_poisson_solver = FastFourier2DSolver(
    boundary_conditions=[[PeriodicBoundaryCondition(),PeriodicBoundaryCondition()],
                         [PeriodicBoundaryCondition(),PeriodicBoundaryCondition()],
                         [PeriodicBoundaryCondition(),PeriodicBoundaryCondition()]]
    )

#----------------------------------------
# Electrode Calculators
#----------------------------------------
left_electrode_calculator = LCAOCalculator(
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
    numerical_accuracy_parameters=left_electrode_numerical_accuracy_parameters,
    poisson_solver=left_electrode_poisson_solver,
    )

right_electrode_calculator = LCAOCalculator(
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
    numerical_accuracy_parameters=right_electrode_numerical_accuracy_parameters,
    poisson_solver=right_electrode_poisson_solver,
    )

#----------------------------------------
# Device Calculator
#----------------------------------------
calculator = DeviceLCAOCalculator(
    basis_set=basis_set,
    exchange_correlation=exchange_correlation,
    numerical_accuracy_parameters=device_numerical_accuracy_parameters,
    electrode_calculators=
        [left_electrode_calculator, right_electrode_calculator],
    )

device_configuration.setCalculator(calculator)


# -------------------------------------------------------------
# Initial State
# -------------------------------------------------------------
initial_spin = InitialSpin(scaled_spins=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   
                                         1.0, 1.0, 1.0, 1.0, 1.0,  1.0 ])
device_configuration.setCalculator(
    calculator,
    initial_spin=initial_spin,
)
device_configuration.update()
nlsave('ParaCobaltMoS2Cobalt.nc', device_configuration)
nlprint(device_configuration)

# -------------------------------------------------------------
# Forces
# -------------------------------------------------------------
forces = Forces(device_configuration)
nlsave('ParaCobaltMoS2Cobalt.nc', forces)
nlprint(forces)

# -------------------------------------------------------------
# Mulliken Population
# -------------------------------------------------------------
mulliken_population = MullikenPopulation(device_configuration)
nlsave('ParaCobaltMoS2Cobalt.nc', mulliken_population)
nlprint(mulliken_population)

# -------------------------------------------------------------
# Transmission Spectrum
# -------------------------------------------------------------

def analyze_transmission(filename):
    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    device_configuration = nlread(filename+'.hdf5', DeviceConfiguration)[0]
    # -------------------------------------------------------------
    # Transmission Spectrum
    # -------------------------------------------------------------
    transmission_spectrum = TransmissionSpectrum(
        configuration=device_configuration,
        energies=numpy.linspace(0,0,1)*eV,
        kpoints=MonkhorstPackGrid(151,151),
        energy_zero_parameter=AverageFermiLevel,
        infinitesimal=1e-06*eV,
        self_energy_calculator=RecursionSelfEnergy(),
        )
    nlsave(filename+'.nc', transmission_spectrum)
    nlprint(transmission_spectrum)

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

analyze_transmission('ParaCobaltMoS2Cobalt')
write_transmissionfile('ParaCobaltMoS2Cobalt')
plot_trans('transmission_ParaCobaltMoS2Cobalt')

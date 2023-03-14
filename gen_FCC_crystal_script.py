import numpy as np
import random

# Change this
ele1 = "Cobalt"
ele2 = "Cobalt"

ele1_comma = str(ele1)+", "
ele2_comma = str(ele2)+", "

# Change this
barrier = "Gr"
ParaorAnti = "Anti"

cent_ele = ""

# Change this
Left_FMLayerNum = 3
Right_FMLayerNum = 3

# Change this
gap_width = 2.04
two_D_spacing = 3.4
Barrier_LayerNum = 3

TwoDMat_UnitCell_AtomNum = 0 
a = 0.
b = 0.
c = 0.

spincomma1 = "1.0, "
spin2 = "1.0"
spincomma2 = "1.0, "



# Set a,b,c
if barrier=="MoS2":
    a = 5.48
    b = 3.165
    cent_ele = "Sulfur, Sulfur, Molybdenum,  Sulfur, Sulfur, Molybdenum,"
    TwoDMat_UnitCell_AtomNum = 6
    BasisSet = "GGABasis.Molybdenum_DoubleZetaPolarized,GGABasis.Sulfur_DoubleZetaPolarized,"

elif barrier=="Gr":
    a = 4.26
    b = 2.46
    cent_ele = "Carbon, Carbon, Carbon, Carbon,"
    TwoDMat_UnitCell_AtomNum = 4
    BasisSet = "GGABasis.Carbon_DoubleZetaPolarized,"

if ele1 == "Cobalt":
    c = 2.04

if ParaorAnti == "Para":
    spin2 = "1.0 "
    spincomma2 = "1.0, "
elif ParaorAnti == "Anti":
    spin2 = "-1.0 "
    spincomma2 = "-1.0, "

Total_Central_AtomNum = Left_FMLayerNum*2+Right_FMLayerNum*2+TwoDMat_UnitCell_AtomNum*Barrier_LayerNum
Filename = ParaorAnti + ele1 + barrier + ele2

# Physical coordinates of atoms
Gr_base = np.array([[  0.71,     1.845,   0.],[  1.42,     0.615,   0.],[  2.84,     0.615,   0.],[  3.55,     1.845,  0.]])
MoS2_base = np.array([[0.9135,2.37375,-1.58],[0.9135,2.37375, 1.58],[  1.827,0.79125, 0.],[  3.654,0.79125,-1.58],[  3.654,0.79125, 1.58],[  4.568,2.37375, 0.]])

if barrier=="MoS2":
    barrier_atom_base = MoS2_base
elif barrier=="Gr":
    barrier_atom_base = Gr_base

# Set coordinates of left and right electrodes of FCC structure
if barrier=="Gr":
    A1 = np.copy(barrier_atom_base[0])
    A2 = np.copy(barrier_atom_base[2])
    C1 = np.copy(barrier_atom_base[1])
    C2 = np.copy(barrier_atom_base[3])
    B1 = (A1 + C2) / 2.
    B2 = A2 + (B1 - A1)

elif barrier=="MoS2":
    A1 = (np.copy(barrier_atom_base[0]) + np.copy(barrier_atom_base[1]))/2.
    A2 = (np.copy(barrier_atom_base[3]) + np.copy(barrier_atom_base[4]))/2.
    C1 = np.copy(barrier_atom_base[2])
    C2 = np.copy(barrier_atom_base[5])
    B1 = (A1 + C2) / 2.
    B2 = A2 + (B1 - A1)

z1 = np.array([0.,0.,c/2.*1.])
z2 = np.array([0.,0.,c/2.*3.])
z3 = np.array([0.,0.,c/2.*5.])

def convert_nparray_to_list(M,dim_space):
    LIST = []
    for i in range(len(M)):
        temp = []
        for j in range(dim_space):
            temp.append(M[i][j])
        LIST.append(temp)
    return LIST

left_electrode  = np.array([A1 + z1 , A2 + z1 , B1 + z2 , B2 + z2 , C1 + z3 , C2 + z3])
right_electrode = np.array([C1 + z1 , C2 + z1 , B1 + z2 , B2 + z2 , A1 + z3 , A2 + z3])

left_electrode_atoms = convert_nparray_to_list(left_electrode,3)
right_electrode_atoms =  convert_nparray_to_list(right_electrode,3)

# Set central atom coordinates
spacing_num = Barrier_LayerNum - 1
tunnel_dist = spacing_num * two_D_spacing + gap_width * 2
tunnel_vec  = np.array([0.,0.,tunnel_dist])

# Atoms on the left and right of the barrier atoms
h = np.array([0.,0.,c*3.])

central_atom_on_the_left = []
central_atom_on_the_right= []

for i in range (int(Left_FMLayerNum/3)):
    temp = np.copy(left_electrode)
    for j in range (len(left_electrode)):
        temp[j] += h * i
        central_atom_on_the_left.append(temp[j])

for i in range (int(Right_FMLayerNum/3)):
    temp = np.copy(right_electrode)
    for j in range (len(right_electrode)):
        temp[j] += h * i
        central_atom_on_the_right.append(temp[j])

for i in range (len(central_atom_on_the_right)):
    central_atom_on_the_right[i] += tunnel_vec + np.array([0.,0.,(Left_FMLayerNum-1)*c])

# Calculate lattice unit z-coord
lattice_z = (Left_FMLayerNum + Right_FMLayerNum -1) * c + tunnel_dist

# 2D-material barrier atoms
Two_D_atoms = []
for i in range(Barrier_LayerNum):
    layer_z  = (Left_FMLayerNum - 0.5) * c + gap_width + i * two_D_spacing
    for j in range(TwoDMat_UnitCell_AtomNum):
        layer_atom = barrier_atom_base[j] + np.array([0.,0.,layer_z])
        Two_D_atoms.append(layer_atom)

# Put all the coordinates into central_all_atoms
central_all_atoms_left = convert_nparray_to_list(central_atom_on_the_left,3)
central_all_atoms_right = convert_nparray_to_list(central_atom_on_the_right,3)
central_all_atoms_mid = convert_nparray_to_list(Two_D_atoms,3)
central_all_atoms = central_all_atoms_left + central_all_atoms_mid + central_all_atoms_right

# Displacement in x-direction
x_dis = barrier_atom_base[0][0]/2.

def displacement(M,displace,dir):
    for i in range(len(M)):
        M[i][dir] += displace

displacement(left_electrode_atoms,-x_dis,0)
displacement(right_electrode_atoms,-x_dis,0)
displacement(central_all_atoms,-x_dis,0)


code = f"""
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
vector_a = [{a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {c*3}]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
left_electrode_elements = [{ele1_comma*5} {ele1}]

# Define coordinates
left_electrode_coordinates = {left_electrode_atoms}*Angstrom

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
vector_a = [{a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {c*3}]*Angstrom
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
right_electrode_elements = [{ele2_comma*5} {ele2}]

# Define coordinates
right_electrode_coordinates = {right_electrode_atoms}*Angstrom

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
vector_a = [{a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {lattice_z}]*Angstrom
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
central_region_elements = [{ele1_comma *Left_FMLayerNum*2}
                           {cent_ele *Barrier_LayerNum}
                           {ele2_comma *(Right_FMLayerNum*2-1)}{ele2}]

# Define coordinates
central_region_coordinates = {central_all_atoms}*Angstrom

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

# Add tags
central_region.addTags('Selection 0', [0, 1, {Total_Central_AtomNum-2}, {Total_Central_AtomNum-1}])

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
    {BasisSet}
    GGABasis.{ele1}_DoubleZetaPolarized,
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
initial_spin = InitialSpin(scaled_spins=[{spincomma1 *Left_FMLayerNum*2} 
                                         {'0.0, ' *Barrier_LayerNum*TwoDMat_UnitCell_AtomNum}  
                                         {spincomma2 *(Right_FMLayerNum*2-1)} {spin2}])
device_configuration.setCalculator(
    calculator,
    initial_spin=initial_spin,
)
device_configuration.update()
nlsave('{Filename}.nc', device_configuration)
nlprint(device_configuration)

# -------------------------------------------------------------
# Forces
# -------------------------------------------------------------
forces = Forces(device_configuration)
nlsave('{Filename}.nc', forces)
nlprint(forces)

# -------------------------------------------------------------
# Mulliken Population
# -------------------------------------------------------------
mulliken_population = MullikenPopulation(device_configuration)
nlsave('{Filename}.nc', mulliken_population)
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

            file.write("Transmission coefficient of spin component " + spin_name[s] +"\\n")
            for i in range(data.shape[0]):
                file.write('Transmission at energy = %12.6f eV \\n' % (energies[i].inUnitsOf(eV)))
                file.write('      kx         ky      transmission \\n')
                for j in range(data.shape[1]):
                    file.write('%10.4f %10.4f %16.6e\\n' % (kpoints[j][0],kpoints[j][1],data[i][j]))

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

analyze_transmission('{Filename}')
write_transmissionfile('{Filename}')
plot_trans('transmission_{Filename}')

"""
pythonfile = Filename+".py"
# Open a file for writing
with open(pythonfile, "w") as f:
    # Write some sample code to the file
    f.write(code)

# Close the file
f.close()
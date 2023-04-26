import numpy as np
import random
import os

class QATK_CodeGenerator:
    
    def __init__(self, **kwargs):
        self.ele1 = kwargs.get('ele1', 'Cobalt')
        self.ele2 = kwargs.get('ele2', 'Cobalt')
        self.barrier = kwargs.get('barrier', 'Gr')
        self.ParaorAnti = kwargs.get('ParaorAnti', 'Para')

        self.gap_width = kwargs.get('gap_width',2.04)
        self.two_D_spacing = kwargs.get('two_D_spacing',3.33)
        self.CrystalPlaneStacking = kwargs.get('CrystalPlaneStacking',[['A','B','C'],['A','B','C']])
        self.BarrierPlaneStacking = kwargs.get('BarrierPlaneStacking',['A'])

        self.FMOxide = kwargs.get('FMOxide',0)
        self.script_dir = kwargs.get('script_dir','./')

        #   Set element specific parameters
        self.ele1_comma = str(self.ele1)+", "
        self.ele2_comma = str(self.ele2)+", "

        if self.barrier=="MoS2":
            self.lattice_const_a = 3.16
            self.lattice_const_c = 3.16
            self.cent_ele = "Sulfur, Sulfur, Molybdenum,  Sulfur, Sulfur, Molybdenum,"
            self.TwoDMat_UnitCell_AtomNum = 6
            self.BasisSet = "GGABasis.Molybdenum_DoubleZetaPolarized,GGABasis.Sulfur_DoubleZetaPolarized,"

        elif self.barrier=="Gr":
            self.lattice_const_a = 2.46
            self.lattice_const_c = 0.
            self.cent_ele = "Carbon, Carbon, Carbon, Carbon,"
            self.TwoDMat_UnitCell_AtomNum = 4
            self.BasisSet = "GGABasis.Carbon_DoubleZetaPolarized,"

        elif self.barrier=="BN":
            self.lattice_const_a = 2.50
            self.lattice_const_c = 0.
            self.cent_ele = "Boron, Nitrogen, Boron, Nitrogen,"
            self.TwoDMat_UnitCell_AtomNum = 4
            self.BasisSet = "GGABasis.Boron_DoubleZetaPolarized,GGABasis.Nitrogen_DoubleZetaPolarized,"

        # Set a,b,c
        self.a = self.lattice_const_a * 3**0.5
        self.b = self.lattice_const_a
        if self.ele1 == "Cobalt":
            self.c = 2.04

        #   Set spins
        self.spincomma1 = "1.0, "

        if self.ParaorAnti == "Para":
            self.spin2 = "1.0 "
            self.spincomma2 = "1.0, "
        elif self.ParaorAnti == "Anti":
            self.spin2 = "-1.0 "
            self.spincomma2 = "-1.0, "

        #   Set oxide related parameters
        if self.FMOxide!=0:
            self.BasisSet += "GGABasis.Oxygen_DoubleZetaPolarized,"
            
        if self.FMOxide==0:
            self.ele_surface1 = self.ele1_comma
            self.ele_surface2 = self.ele2_comma
            self.spin_surface1 = self.spincomma1
            self.spin_surface2 = self.spincomma2

        elif self.FMOxide==1:
            self.ele_surface1 = "Oxygen,"
            self.ele_surface2 = self.ele2_comma
            self.spin_surface1 = "0.0,"
            self.spin_surface2 = self.spincomma2
        
        elif self.FMOxide==2:
            self.ele_surface1 = "Oxygen,"
            self.ele_surface2 = "Oxygen,"
            self.spin_surface1 = "0.0,"
            self.spin_surface2 = "0.0,"
        
        #   Set crystal layer numbers
        self.Left_FMLayerNum = len(self.CrystalPlaneStacking[0])
        self.Right_FMLayerNum = len(self.CrystalPlaneStacking[1])
        self.Barrier_LayerNum = len(self.BarrierPlaneStacking)

        #   Set atom number in central arae
        self.Total_Central_AtomNum = self.Left_FMLayerNum*2+self.Right_FMLayerNum*2+self.TwoDMat_UnitCell_AtomNum*self.Barrier_LayerNum

        #   Set output filename
        self.Filename = self.ParaorAnti + self.ele1 + self.barrier + self.ele2
        if self.FMOxide == 1:
            self.Filename = "Oxide" + self.Filename
        elif self.FMOxide == 2:
            self.Filename = "BothOxide" + self.Filename

        #   Set atom coordinates
        self.left_electrode_atoms, self.right_electrode_atoms, self.central_atoms, self.lattice_z = self.set_coordinates()
        
        #   Generate python file
        self.generate_python_code_file(self.Filename,self.generate_code_text())

    def set_coordinates(self):
        # Physical coordinates of atoms
        if self.barrier=="Gr":
            Barrier_Type = 'GrType'
        elif self.barrier=="MoS2":
            Barrier_Type = 'MoS2Type'
        elif self.barrier=="BN":
            Barrier_Type = 'GrType'

        lc = self.lattice_const_a
        lc2 = self.lattice_const_c

        if Barrier_Type=='GrType':
            barrier_atom_base = np.array([[  lc/(3**0.5)/2,     lc/4.*3,   0.],
                                        [  lc/(3**0.5),       lc/4.,     0.],
                                        [  lc/(3**0.5)*2,     lc/4.,     0.],
                                        [  lc/(3**0.5)*2.5,   lc/4.*3,   0.]])
        elif Barrier_Type=='MoS2Type':
            barrier_atom_base = np.array([[  lc/(3**0.5)/2,     lc/4.*3,   -lc2/2.],
                                        [  lc/(3**0.5)/2,     lc/4.*3,    lc2/2.],
                                        [  lc/(3**0.5),       lc/4.,     0.],
                                        [  lc/(3**0.5)*2,     lc/4.,     -lc2/2.],
                                        [  lc/(3**0.5)*2,     lc/4.,      lc2/2.],
                                        [  lc/(3**0.5)*2.5,   lc/4.*3,   0.]])

        # Set basis coordinates

        if Barrier_Type=='GrType':
            A1 = np.copy(barrier_atom_base[0])
            A2 = np.copy(barrier_atom_base[2])
            C1 = np.copy(barrier_atom_base[1])
            C2 = np.copy(barrier_atom_base[3])
            B1 = (A1 + C2) / 2.
            B2 = A2 + (B1 - A1)

        elif Barrier_Type=='MoS2Type':
            A1 = (np.copy(barrier_atom_base[0]) + np.copy(barrier_atom_base[1]))/2.
            A2 = (np.copy(barrier_atom_base[3]) + np.copy(barrier_atom_base[4]))/2.
            C1 = np.copy(barrier_atom_base[2])
            C2 = np.copy(barrier_atom_base[5])
            B1 = (A1 + C2) / 2.
            B2 = A2 + (B1 - A1)

        central_left_FM_plane_type  = self.CrystalPlaneStacking[0]
        central_right_FM_plane_type = self.CrystalPlaneStacking[1]

        zb = np.array([0.,0.,self.c])
        z_barrier_dis = np.array([0.,0.,self.gap_width + self.c * len(central_left_FM_plane_type)])
        z_right_dis = np.array([0.,0.,(self.Barrier_LayerNum - 1) * self.two_D_spacing + self.gap_width * 2 + self.c * len(central_left_FM_plane_type)])

        def displace(x,dx):
            y = np.copy(x)
            for i in range(len(x)):
                y[i] += dx
            return y

        def set_FM_crystal_atom(plane_type):
            coord_list = []
            for i in range(len(plane_type)):
                if plane_type[i]=='A':
                    coord_list.append(A1 + zb*i)
                    coord_list.append(A2 + zb*i)
                elif plane_type[i]=='B':
                    coord_list.append(B1 + zb*i)
                    coord_list.append(B2 + zb*i)
                elif plane_type[i]=='C':
                    coord_list.append(C1 + zb*i)
                    coord_list.append(C2 + zb*i)
            return np.array(coord_list)

        def set_barrier_atom(plane_type,barrier_atom_base):
            for i in range(len(plane_type)):
                if i == 0:
                    if plane_type[i]=='A':
                        bab = np.copy(barrier_atom_base)
                        
                    elif plane_type[i]=='B':
                        bab = np.copy(-barrier_atom_base)
                        bab = displace(bab,np.array([self.a,self.b,0.]))

                    bab = displace(bab, np.array([0.,0.,i*self.two_D_spacing]))
                    coord_list = np.copy(bab)

                else:
                    if plane_type[i]=='A':
                        bab = np.copy(barrier_atom_base)
                        
                    elif plane_type[i]=='B':
                        bab = np.copy(-barrier_atom_base)
                        bab = displace(bab,np.array([self.a,self.b,0.]))
                        
                    bab = displace(bab, np.array([0.,0.,i*self.two_D_spacing]))
                    coord_list = np.concatenate((coord_list,np.copy(bab)))
                    
            return coord_list

        left_electrode_atoms = set_FM_crystal_atom([central_left_FM_plane_type[0],central_left_FM_plane_type[1],central_left_FM_plane_type[2]])
        central_left_atoms = set_FM_crystal_atom(central_left_FM_plane_type)
        central_right_atoms = set_FM_crystal_atom(central_right_FM_plane_type)
        central_barrier_atoms = set_barrier_atom(self.BarrierPlaneStacking,barrier_atom_base)
        right_electrode_atoms = set_FM_crystal_atom([central_right_FM_plane_type[-3],central_right_FM_plane_type[-2],central_right_FM_plane_type[-1]])

        left_electrode_atoms = displace(left_electrode_atoms , zb/2.)
        central_left_atoms = displace(central_left_atoms , zb/2.)
        central_right_atoms = displace(central_right_atoms , zb/2. + z_right_dis)
        central_barrier_atoms = displace(central_barrier_atoms , z_barrier_dis)
        right_electrode_atoms = displace(right_electrode_atoms , zb/2.)

        def convert_nparray_to_list(M,dim_space):
                LIST = []
                for i in range(len(M)):
                    temp = []
                    for j in range(dim_space):
                        temp.append(M[i][j])
                    LIST.append(temp)
                return LIST

        left_electrode_atoms = convert_nparray_to_list(left_electrode_atoms , 3)
        central_left_atoms = convert_nparray_to_list(central_left_atoms , 3)
        central_right_atoms = convert_nparray_to_list(central_right_atoms , 3)
        central_barrier_atoms = convert_nparray_to_list(central_barrier_atoms , 3)
        right_electrode_atoms = convert_nparray_to_list(right_electrode_atoms , 3)

        central_atoms = central_left_atoms + central_barrier_atoms + central_right_atoms
        lattice_z = self.c * (self.Left_FMLayerNum + self.Right_FMLayerNum) + (self.Barrier_LayerNum - 1) * self.two_D_spacing + self.gap_width * 2
        
        return left_electrode_atoms, right_electrode_atoms, central_atoms, lattice_z


    def generate_python_code_file(self,filename,code_text):
        pythonfile = self.script_dir + '/' + filename+".py"
        # Open a file for writing
        with open(pythonfile, "w") as f:
            # Write some sample code to the file
            f.write(code_text)

        # Close the file
        f.close()

    def generate_code_text(self):
        code_text = f"""
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
vector_a = [{self.a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {self.b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {self.c*3}]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
left_electrode_elements = [{self.ele1_comma*5} {self.ele1}]

# Define coordinates
left_electrode_coordinates = {self.left_electrode_atoms}*Angstrom

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
vector_a = [{self.a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {self.b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {self.c*3}]*Angstrom
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
right_electrode_elements = [{self.ele2_comma*5} {self.ele2}]

# Define coordinates
right_electrode_coordinates = {self.right_electrode_atoms}*Angstrom

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
vector_a = [{self.a}, 0.0, 0.0]*Angstrom
vector_b = [0.0, {self.b}, 0.0]*Angstrom
vector_c = [0.0, 0.0, {self.lattice_z}]*Angstrom
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
central_region_elements = [{self.ele1_comma *(self.Left_FMLayerNum*2-1)}{self.ele_surface1}
                        {self.cent_ele *self.Barrier_LayerNum}
                        {self.ele_surface2}{self.ele2_comma *(self.Right_FMLayerNum*2-2)}{self.ele2}]

# Define coordinates
central_region_coordinates = {self.central_atoms}*Angstrom

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

# Add tags
central_region.addTags('Selection 0', [0, 1, {self.Total_Central_AtomNum-2}, {self.Total_Central_AtomNum-1}])

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
    {self.BasisSet}
    GGABasis.{self.ele1}_DoubleZetaPolarized,
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
initial_spin = InitialSpin(scaled_spins=[{self.spincomma1 *(self.Left_FMLayerNum*2-1)} {self.spin_surface1}
                                        {'0.0, ' *self.Barrier_LayerNum*self.TwoDMat_UnitCell_AtomNum}  
                                        {self.spin_surface2}{self.spincomma2 *(self.Right_FMLayerNum*2-2)} {self.spin2}])
device_configuration.setCalculator(
    calculator,
    initial_spin=initial_spin,
)
device_configuration.update()
nlsave('{self.Filename}.nc', device_configuration)
nlprint(device_configuration)

# -------------------------------------------------------------
# Forces
# -------------------------------------------------------------
forces = Forces(device_configuration)
nlsave('{self.Filename}.nc', forces)
nlprint(forces)

# -------------------------------------------------------------
# Mulliken Population
# -------------------------------------------------------------
mulliken_population = MullikenPopulation(device_configuration)
nlsave('{self.Filename}.nc', mulliken_population)
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

analyze_transmission('{self.Filename}')
write_transmissionfile('{self.Filename}')
plot_trans('transmission_{self.Filename}')

"""
        return code_text

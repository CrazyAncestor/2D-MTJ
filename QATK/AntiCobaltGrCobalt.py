# -*- coding: utf-8 -*-
# -------------------------------------------------------------
# Two-probe Configuration
# -------------------------------------------------------------

# -------------------------------------------------------------
# Left Electrode
# -------------------------------------------------------------

# Set up lattice
vector_a = [4.26, 0.0, 0.0]*Angstrom
vector_b = [0.0, 2.46, 0.0]*Angstrom
vector_c = [0.0, 0.0, 6.12]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
left_electrode_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt,  Cobalt]

# Define coordinates
left_electrode_coordinates = [[0.355, 1.845, 1.02], [2.485, 0.615, 1.02], [1.775, 1.845, 3.06], [3.905, 0.615, 3.06], [1.065, 0.615, 5.1], [3.195, 1.845, 5.1]]*Angstrom

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
vector_a = [4.26, 0.0, 0.0]*Angstrom
vector_b = [0.0, 2.46, 0.0]*Angstrom
vector_c = [0.0, 0.0, 6.12]*Angstrom
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
right_electrode_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt,  Cobalt]

# Define coordinates
right_electrode_coordinates = [[1.065, 0.615, 1.02], [3.195, 1.845, 1.02], [1.775, 1.845, 3.06], [3.905, 0.615, 3.06], [0.355, 1.845, 5.1], [2.485, 0.615, 5.1]]*Angstrom

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
vector_a = [4.26, 0.0, 0.0]*Angstrom
vector_b = [0.0, 2.46, 0.0]*Angstrom
vector_c = [0.0, 0.0, 14.28]*Angstrom
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
central_region_elements = [Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, 
                           Carbon, Carbon, Carbon, Carbon,
                           Cobalt, Cobalt, Cobalt, Cobalt, Cobalt, Cobalt]

# Define coordinates
central_region_coordinates = [[0.355, 1.845, 1.02], [2.485, 0.615, 1.02], [1.775, 1.845, 3.06], [3.905, 0.615, 3.06], [1.065, 0.615, 5.1], [3.195, 1.845, 5.1], [0.355, 1.845, 7.14], [1.065, 0.615, 7.14], [2.485, 0.615, 7.14], [3.195, 1.845, 7.14], [1.065, 0.615, 9.18], [3.195, 1.845, 9.18], [1.775, 1.845, 11.22], [3.905, 0.615, 11.22], [0.355, 1.845, 13.26], [2.485, 0.615, 13.26]]*Angstrom

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

# Add tags
central_region.addTags('Selection 0', [0, 1, 14, 15])

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
    GGABasis.Carbon_DoubleZetaPolarized,
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
                                         0.0, 0.0, 0.0, 0.0,   
                                         -1.0, -1.0, -1.0, -1.0, -1.0,  -1.0 ])
device_configuration.setCalculator(
    calculator,
    initial_spin=initial_spin,
)
device_configuration.update()
nlsave('AntiCobaltGrCobalt.nc', device_configuration)
nlprint(device_configuration)

# -------------------------------------------------------------
# Forces
# -------------------------------------------------------------
forces = Forces(device_configuration)
nlsave('AntiCobaltGrCobalt.nc', forces)
nlprint(forces)

# -------------------------------------------------------------
# Mulliken Population
# -------------------------------------------------------------
mulliken_population = MullikenPopulation(device_configuration)
nlsave('AntiCobaltGrCobalt.nc', mulliken_population)
nlprint(mulliken_population)

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
nlsave('AntiCobaltGrCobalt.nc', transmission_spectrum)
nlprint(transmission_spectrum)

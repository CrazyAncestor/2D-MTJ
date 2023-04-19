
def give_parameters(AntiPara,FM_ele1,FM_ele2,TwoD_barrier):
    filename = AntiPara + FM_ele1 + TwoD_barrier + FM_ele2

    if TwoD_barrier == 'Gr':
        lattice_a = 2.46
        lattice_b = 4.26
    
    # K-point
    x0 = 2*np.pi/3/lattice_a
    y0 = 2*np.pi/3/lattice_a / 3**0.5
    z0 = 4*np.pi/3/lattice_a / 3**0.5

    x0 = 2*np.pi/lattice_b - x0 # folded by lattice vector
    x0 /= np.pi/lattice_b*2.
    y0 /= np.pi/lattice_a*2.
    z0 /= np.pi/lattice_a*2.
    K_points = [[x0,y0,0.],[-x0,y0,0.],[x0,-y0,0.],[-x0,-y0,0.],[0,z0,0.],[0,-z0,0.]]

    return filename , K_points

def band_analysis(filename,K_points,kp_num,energy_value,atoms=None, elements=None, angular_momenta=None):

    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    device_configuration = nlread(filename+'.hdf5', DeviceConfiguration)[0]

    #   Set projection list
    if elements==None:
        projection_list = ProjectionList(atoms=atoms,angular_momenta=angular_momenta)
    elif atoms==None:
        projection_list = ProjectionList(elements=elements,angular_momenta=angular_momenta)
    # -------------------------------------------------------------
    # Surface Bandstructure
    # -------------------------------------------------------------
    surface_bandstructure = SurfaceBandstructure(
        configuration=device_configuration,
        kpoints=[[0.,0.,0.],K_points[kp_num[0]],K_points[kp_num[1]]]
        points_per_segment=20,
        energies=numpy.linspace(-energy_value, energy_value, 101)*eV,
        projection_list=projection_list,
        contributions=All,
        self_energy_calculator=RecursionSelfEnergy(storage_strategy=NoStorage()),
        energy_zero_parameter=AverageFermiLevel,
        infinitesimal=1e-06*eV,
        )
    nlsave('BandStructure_'+filename+'.hdf5', surface_bandstructure)
    nlprint(surface_bandstructure)

if __name__=='main':
    FM_ele1 = 'Cobalt'
    FM_ele2 = 'Cobalt'
    TwoD_barrier = 'Gr'
    AntiPara = 'Para'
    kp_num = [0,5]
    energy_value = 2.0
    atoms = All
    angular_momenta = [1]

    filename , K_points = give_parameters(AntiPara,FM_ele1,FM_ele2,TwoD_barrier)
    band_analysis(filename,K_points,kp_num,energy_value,atoms=atoms, elements=None, angular_momenta=angular_momenta)

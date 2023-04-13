structure_para = ''
structure_anti = ''
Energy_value = 2.0

structures = [structure_para,structure_anti]

for structure in structures:
    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    device_configuration = nlread(structure+'.hdf5', DeviceConfiguration)[0]


    # -------------------------------------------------------------
    # Surface Bandstructure
    # -------------------------------------------------------------
    surface_bandstructure = SurfaceBandstructure(
        configuration=device_configuration,
        route=['G', 'X', 'Y'],
        points_per_segment=20,
        energies=numpy.linspace(-Energy_value, Energy_value, 101)*eV,
        projection_list=ProjectionList(atoms=All),
        contributions=All,
        self_energy_calculator=RecursionSelfEnergy(storage_strategy=NoStorage()),
        energy_zero_parameter=AverageFermiLevel,
        infinitesimal=1e-06*eV,
        )
    nlsave(structure+'Band.hdf5', surface_bandstructure)
    nlprint(surface_bandstructure)
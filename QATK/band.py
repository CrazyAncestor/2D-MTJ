import numpy as np
from BandStructurePlotter import BandStructurePlotter

def give_parameters(AntiPara,FM_ele1,FM_ele2,TwoD_barrier,kp_num,elements,angular_momenta):
    structure = AntiPara + FM_ele1 + TwoD_barrier + FM_ele2

    if TwoD_barrier == 'Gr':
        lattice_a = 2.46
        lattice_b = 4.26
    elif TwoD_barrier == 'BN':
        lattice_a = 2.5
        lattice_b = 4.3
    elif TwoD_barrier == 'MoS2':
        lattice_a = 3.16
        lattice_b = 5.47
    
    # K-point
    x0 = 2*np.pi/3/lattice_a
    y0 = 2*np.pi/3/lattice_a / 3**0.5
    z0 = 4*np.pi/3/lattice_a / 3**0.5

    x0 = 2*np.pi/lattice_b - x0 # folded by lattice vector
    x0 /= np.pi/lattice_b*2.
    y0 /= np.pi/lattice_a*2.
    z0 /= np.pi/lattice_a*2.
    K_nodes = [[x0,y0,0.],[-x0,y0,0.],[x0,-y0,0.],[-x0,-y0,0.],[0,z0,0.],[0,-z0,0.]]
    G = [0.,0.,0.]
    
    def make_segment(p1,p2,my_list,num):
        xs = np.linspace(p1[0],p2[0],num)
        ys = np.linspace(p1[1],p2[1],num)
        zs = np.linspace(p1[2],p2[2],num)
        for i in range(num):
            my_list.append([xs[i],ys[i],zs[i]])

    def mid_point(p1,p2):
        p = []
        for i in range(3):
            p.append((p1[i]+p2[i])/2.)
        return p

    K_points = []
    make_segment(G,K_nodes[0],K_points,kp_num)
    make_segment(K_nodes[0],mid_point(K_nodes[0],K_nodes[2]),K_points,kp_num)

    if elements[0]==Carbon:
        elements_name = 'Carbon'
    elif elements[0]==Cobalt:
        elements_name = 'Cobalt'
    elif elements[0]==Boron and elements[1]==Nitrogen:
        elements_name = 'BN'
    elif elements[0]==Molybdenum and elements[1]==Sulfur:
        elements_name = 'MoS2'

    if angular_momenta[0]==1:
        angular_momenta_name = 'p_orbital'
    elif angular_momenta[0]==2:
        angular_momenta_name = 'd_orbital'

    return structure , K_points, elements_name, angular_momenta_name

def band_analysis(structure,outputfile,K_points,energy_value,atoms=None, elements=None, angular_momenta=None):

    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    device_configuration = nlread(structure+'.hdf5', DeviceConfiguration)[0]

    projection_list = ProjectionList(atoms=atoms,elements=elements,angular_momenta=angular_momenta)

    # -------------------------------------------------------------
    # Surface Bandstructure
    # -------------------------------------------------------------
    surface_bandstructure = SurfaceBandstructure(
        configuration=device_configuration,
        kpoints=K_points,
        energies=numpy.linspace(-energy_value, energy_value, 101)*eV,
        projection_list=projection_list,
        contributions=All,
        self_energy_calculator=RecursionSelfEnergy(storage_strategy=NoStorage()),
        energy_zero_parameter=AverageFermiLevel,
        infinitesimal=1e-06*eV,
        )
    nlsave(outputfile, surface_bandstructure)
    nlprint(surface_bandstructure)

def write_SurfaceBandstructurefile(inputfile,outputfile):
    # -------------------------------------------------------------
    # Load device configuration
    # -------------------------------------------------------------
    SurfaceBandstructure_spectrum = nlread(inputfile, SurfaceBandstructure)[0]
    # Print out all k-dependent SurfaceBandstructure coefficients
    spin_name = ["Up","Down"]
    spins = [Spin.Up,Spin.Down] 

    print(SurfaceBandstructure_spectrum)
    #print data
    for s in range(len(spin_name)):
        with open(spin_name[s]+outputfile, "w") as file:

            data = 10#SurfaceBandstructure_spectrum.zData(spin = spins[s])
            energies = SurfaceBandstructure_spectrum.energies()
            kpoints = SurfaceBandstructure_spectrum.kpoints()
            ss = SurfaceBandstructure_spectrum.evaluate(spins[s])
            print(ss)
            for i in range(len(ss)):
                for j in range(len(ss[i])):
                    file.write('%12.6f' % (ss[i][j]))
                file.write('\n')

#      Main Function
#   Set parameters
FM_ele1 = 'Cobalt'
FM_ele2 = 'Cobalt'
TwoD_barrier = 'Gr'
AntiPara = 'Para'
kp_num = 50
energy_value = 2.0
elements = [Carbon]
angular_momenta = [1]
plotting_max = 0.2

#   Band analysis
structure, K_points, elements_name, angular_momenta_name = give_parameters(AntiPara,FM_ele1,FM_ele2,TwoD_barrier,kp_num,elements,angular_momenta)
HDF5_Output = 'BandStructure_'+structure+elements_name+angular_momenta_name+'.hdf5'
band_analysis(structure,HDF5_Output,K_points,energy_value,atoms=None, elements=elements, angular_momenta=angular_momenta)

#   Export band data
TXT_Output = 'BandStructure_'+structure+elements_name+angular_momenta_name+'.txt'
write_SurfaceBandstructurefile(HDF5_Output,K_points,TXT_Output)


#   Plot band structure
fig_name = structure+elements_name+angular_momenta_name+'.png'
a = BandStructurePlotter(fig_dir='band',input_filename='Up'+TXT_Output,output_filename='Up'+fig_name,max_value=plotting_max)
b = BandStructurePlotter(fig_dir='band',input_filename='Down'+TXT_Output,output_filename='Down'+fig_name,max_value=plotting_max)


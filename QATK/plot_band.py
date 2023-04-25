from BandStructurePlotter import BandStructurePlotter
import numpy as np
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
    
    if len(elements)==1:
        if elements[0]=='Carbon':
            elements_name = 'Carbon'
        elif elements[0]=='Cobalt':
            elements_name = 'Cobalt'
        elif elements[0]=='Molybdenum':
            elements_name = 'Molybdenum'
        elif elements[0]=='Sulfur':
            elements_name = 'Sulfur'
    elif len(elements)==2:
        if elements[0]=='Boron' and elements[1]=='Nitrogen':
            elements_name = 'BN'
        elif elements[0]=='Molybdenum' and elements[1]=='Sulfur':
            elements_name = 'MoS2'
    

    if angular_momenta[0]==1:
        angular_momenta_name = 'p_orbital'
    elif angular_momenta[0]==2:
        angular_momenta_name = 'd_orbital'

    return structure , K_points, elements_name, angular_momenta_name

def plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'Gr', AntiPara = 'Para', kp_num = 50, energy_value = 2.0,\
    elements = ['Carbon'], angular_momenta = [1], plotting_max = 0.2, dir_name = 'Gr'):
    #   Band analysis
    structure, K_points, elements_name, angular_momenta_name = give_parameters(AntiPara,FM_ele1,FM_ele2,TwoD_barrier,kp_num,elements,angular_momenta)
    TXT_Output = 'BandStructure_'+structure+elements_name+angular_momenta_name+'.txt'
    fig_name = structure+elements_name+angular_momenta_name+'.png'

    a = BandStructurePlotter(fig_dir=dir_name,data_dir=dir_name,input_filename=TXT_Output,output_filename=fig_name,energy_range=[-energy_value,energy_value],max_value=plotting_max)

if __name__=='__main__':
#      Gr
    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'Gr', AntiPara = 'Para', kp_num = 50, energy_value = 2.0,\
    elements = ['Carbon'], angular_momenta = [1], plotting_max = 0.2, dir_name = 'Gr')

    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'Gr', AntiPara = 'Para', kp_num = 50, energy_value = 2.0,\
    elements = ['Cobalt'], angular_momenta = [2], plotting_max = 20., dir_name = 'Gr')

#      BN-small
    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'BN', AntiPara = 'Para', kp_num = 50, energy_value = 2.0,\
    elements = ['Boron','Nitrogen'], angular_momenta = [1], plotting_max = 0.2, dir_name = 'BN_small')

    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'BN', AntiPara = 'Para', kp_num = 50, energy_value = 2.0,\
    elements = ['Cobalt'], angular_momenta = [2], plotting_max = 20., dir_name = 'BN_small')

#      BN-large
    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'BN', AntiPara = 'Para', kp_num = 50, energy_value = 10.0,\
    elements = ['Boron','Nitrogen'], angular_momenta = [1], plotting_max = 0.2, dir_name = 'BN_large')

    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'BN', AntiPara = 'Para', kp_num = 50, energy_value = 10.0,\
    elements = ['Cobalt'], angular_momenta = [2], plotting_max = 20., dir_name = 'BN_large')
    
#      MoS2
    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'MoS2', AntiPara = 'Para', kp_num = 50, energy_value = 4.0,\
    elements = ['Molybdenum'], angular_momenta = [2], plotting_max = 20., dir_name = 'MoS2')

    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'MoS2', AntiPara = 'Para', kp_num = 50, energy_value = 4.0,\
    elements = ['Cobalt'], angular_momenta = [2], plotting_max = 20., dir_name = 'MoS2')

    plot_band(FM_ele1 = 'Cobalt', FM_ele2 = 'Cobalt', TwoD_barrier = 'MoS2', AntiPara = 'Para', kp_num = 50, energy_value = 4.0,\
    elements = ['Sulfur'], angular_momenta = [1], plotting_max = 1., dir_name = 'MoS2')

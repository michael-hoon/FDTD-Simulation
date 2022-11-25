import matplotlib.pyplot as plt
import numpy as np
import meep as mp

#Setting up initial conds
# wvl = 0.6
wvl = [0.35, 0.4, 0.45, 0.5, 0.6]
# frq = 1/wvl 
fnumber = 4.0
t_width = 4.5
# s_width = 2.0 * wvl * fnumber / np.pi

#set up box
dpml = 1
dair = 20
s = dair + 2*dpml
cell_size = mp.Vector3(s, s)
resolution = 50
pml_layers = [mp.PML(thickness=dpml)]
symmetries = [mp.Mirror(mp.Y)] #???
#Setting up simulation

for w in wvl:

    # is_integrated=True necessary for any planewave source extending into PML
    sources = [
        mp.GaussianBeamSource(
            mp.GaussianSource(wavelength = w, width = t_width,cutoff = t_width * 2,is_integrated = True),
            center = mp.Vector3(-dair/2, 0), 
            size = mp.Vector3(0, s), #line of y = -dair/2 (vert line through entire box)
            beam_x0 = mp.Vector3(dair/2, 0), #Focus at center of box
            beam_kdir = mp.Vector3(1, 0), #pos x dir
            beam_w0 = 2.0 * w * fnumber / np.pi, #spatial width
            beam_E0 = mp.Vector3(0, 0, 1)
        )]

    #create silicon sphere at center
    r = 1.0
    lorentz_suscep = [mp.LorentzianSusceptibility(frequency = 1/0.354241, gamma = 1/5.42347, sigma = 2.37226548)]
    material = mp.Medium(epsilon = 1.0001, E_susceptibilities = lorentz_suscep)
    geometry = [
        mp.Sphere(material=material, center=mp.Vector3(), radius=r)
    ]

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=mp.Vector3(),
        symmetries=symmetries,
        geometry=geometry,
    )

    # def output(sim, todo):
    #     if todo == 'step':
    #         time = sim.meep_time()
    #         slice = sim.get_array(component = mp.Ex, center = (0,0,0), size = cell_size)
            # plt.imshow(slice)
            # plt.savefig(f'{time}.png')
            # plt.colorbar()
            # plt.clf()
        # if todo == 'finish':
        #     pass

    #Run Simulation and obtain data
    sim.use_output_directory(f"Test_{w}") #output all to default directory mie_scattering-out/
    plt.figure()
    sim.run(
        mp.at_beginning(mp.output_epsilon),
        # mp.at_every(1, output),
        # mp.to_appended("efield_z", mp.at_every(1, mp.output_efield_z)),
        mp.at_every(1, mp.output_png(mp.Ez, "-Zc /home/maikuhl/.conda/envs/mp/share/h5utils/colormaps/dkbluered -C $EPS")),
        until=60)


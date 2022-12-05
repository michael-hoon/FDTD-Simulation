import matplotlib.pyplot as plt
import numpy as np
import meep as mp

#Setting up initial conds
# wvl = 0.6
wvl = [0.35, 0.4, 0.45, 0.5, 0.6]
# frq = 1/wvl 
fnumber = 4.0
t_width = 4.5/np.sqrt(4*np.log(2))
t_cutoff = np.sqrt(4*np.log(2))*2.0
# s_width = 2.0 * wvl * fnumber / np.pi

#set up box
dpml = 1
dair = 20
s = dair + 2*dpml
cell_size = mp.Vector3(s, s)
resolution = 50
pml_layers = [mp.PML(thickness=dpml)]
#Setting up simulation

for w in wvl:

    # is_integrated=True necessary for any planewave source extending into PML
    sources = [
        mp.GaussianBeamSource(
            mp.GaussianSource(wavelength = w, width = t_width,cutoff = t_cutoff,is_integrated = True),
            center = mp.Vector3(-dair/2, 0), 
            size = mp.Vector3(0, s), #line of y = -dair/2 (vert line through entire box)
            beam_x0 = mp.Vector3(dair/2, 0), #Focus at center of box
            beam_kdir = mp.Vector3(1, 0), #pos x dir
            beam_w0 = 2.0 * w * fnumber / np.pi, #spatial width
            beam_E0 = mp.Vector3(0, 0, 1)
        )]

    #create silicon sphere at center
    r = 1.0
    lorentz_suscep = [
        mp.LorentzianSusceptibility(frequency = 1/0.354241, gamma = 1/5.42347, sigma = 2.37226548)
        # mp.LorentzianSusceptibility(frequency = 1/0.354241, gamma = 1/5.42347, sigma = 1.61254),
        # mp.LorentzianSusceptibility(frequency = 1/0.29173, gamma = 1/1.180, sigma = 10.1351)
        ]
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
        geometry=geometry,
    )

    #Run Simulation and obtain data
    sim.use_output_directory(f"Test3_{w}") #output all to default directory mie_scattering-out/
    sim.run(
        mp.at_beginning(mp.output_epsilon),
        # mp.to_appended(f"efield_z_{w}", mp.at_every(1, mp.output_efield_z)),
        # mp.at_every(1, mp.output_png(mp.Ez, "-Zc /home/maikuhl/.conda/envs/mp/share/h5utils/colormaps/dkbluered -C $EPS")), #Using h5topng --help to find path due to problem with h5utils that is compatible with meep
        mp.at_every(1, mp.output_png(mp.Ez, "-Zc /home/draco/miniconda3/envs/mp/share/h5utils/colormaps/dkblueredblack -C $EPS")),
        until=40)


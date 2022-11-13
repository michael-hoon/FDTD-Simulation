import matplotlib.pyplot as plt
import numpy as np
import os

import meep as mp

#Setting up initial conds
wvl = 0.6 #Test purposes
# wvl = [0.35, 0.4, 0.45, 0.5, 0.6]
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
symmetries = [mp.Mirror(mp.Y), mp.Mirror(mp.Z, phase=-1)] #Phase = -1 is odd mirror? Do we need this?

#Setting up simulation

#For wvl in wvl:

# is_integrated=True necessary for any planewave source extending into PML
sources = [
    mp.GaussianBeamSource(
        mp.GaussianSource(wavelength = wvl, width = t_width,cutoff = t_width * 2,is_integrated = True),
        center = mp.Vector3(-dair/2, 0), 
        size = mp.Vector3(0, s/2), #line of y = -dair/2 (vert line through entire box)
        beam_x0 = mp.Vector3(dair/2, 0), #Focus at center of box
        beam_kdir = mp.Vector3(1, 0), #pos x dir
        beam_w0 = 2.0 * wvl * fnumber / np.pi, #spatial width
    )]

#create silicon sphere at center
r = 1.0
peak = 0.354241
gamma = 5.42347
sigma = 460.149
material = mp.Medium(index = 2) #to model lorentzian resonance
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

#Run Simulation and obtain data
sim.use_output_directory() #output all to default directory mie_scattering-out/
sim.run(
    mp.to_appended("efield_x", mp.at_every(1, mp.output_efield_x)),
    # mp.to_appended("b_field", mp.at_every(1, mp.output_bfield)),
    until=50)

# os.system("convert mie_scattering-out/ez.t*.png ez.gif")

import matplotlib.pyplot as plt
import numpy as np
import meep as mp

wvls = [0.6,0.5,0.45,0.4,0.35];

# basically, 4.5 is the full-width-half-maximum intensity. See the gaussian
# beam on wiki. To convert, divide by the given factor and an extra sqrt(2),
# hence the factor of np.sqrt(4*np.log(2)) instead of np.sqrt(2*np.log(2)) on wiki.
# The reason is meep using a gaussian like
# exp( - x^2 / (2*w^2) )
# instead of that used by wikipedia (and most laser physicists)
# exp( - x^2 / w^2)
t_width = 4.5/np.sqrt(4*np.log(2))

# t_cutoff is multiplied by the width by the mp.Source object,
# need to cancel out crap to give just a factor of 2.0, and thus have
# a cut off at 9 micron
t_cutoff = np.sqrt(4*np.log(2))*2.0;

#set up box
dpml = 1
dair = 20
s = dair + 2*dpml
cell_size = mp.Vector3(s, s);
resolution = 50
pml_layers = [mp.PML(thickness=dpml)]
symmetries = [] #no symmetries, just a 2D simulation.

#loop over wavelengths
for w in wvls:
    # is_integrated=True necessary for any planewave source extending into PML
    sources = [
        mp.Source(
            mp.GaussianSource(wavelength = w, width = t_width,cutoff = t_cutoff,is_integrated = True),
            center = mp.Vector3(-dair/2, 0), 
            size = mp.Vector3(0, s),
            component = mp.Ey,
        )]
    #create single lorentzian material sphere at center
    r = 1.0
    lorentz_suscep = [mp.LorentzianSusceptibility(
        frequency = 1/0.354241, gamma = 1/5.42347, sigma = 2.37226548)]
    material = mp.Medium(epsilon = 1.001, E_susceptibilities = lorentz_suscep)
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
    sim.use_output_directory(f"mainmeep_{w:.3f}")
    sim.run(
        mp.to_appended("ey",
                       mp.at_beginning(mp.output_efield_y),
                       mp.at_every(1, mp.output_efield_y)),
        #sorry guys, if you also want to output a video, be sure to add your needed hack :)
        until=30)

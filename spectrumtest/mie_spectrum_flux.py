import matplotlib.pyplot as plt
import numpy as np
import meep as mp

# wvls = [0.6,0.5,0.45,0.4,0.35];
wvls = [0.35, 0.40, 0.45]

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

energy = {}
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
    #Set up boundary region to measure energy spectrum
    wvlmax = 0.7
    wvlmin = 0.3
    frqmax = 1 / wvlmax
    frqmin = 1 / wvlmin
    frq_cen = (frqmax + frqmin) / 2
    dfrq = frqmax - frqmin
    nfrq = 100

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=mp.Vector3(),
    )

    # flux_box = sim.add_flux(frq_cen, dfrq, nfrq,
    #     mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, 2 * r)),
    #     mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, 2 * r)),
    #     mp.FluxRegion(center = mp.Vector3(y = -r), size = mp.Vector3(2 * r, 0)),
    #     mp.FluxRegion(center = mp.Vector3(y = +r), size = mp.Vector3(2 * r, 0)))

    incident_line = sim.add_flux(frq_cen, dfrq, nfrq,
        mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, s)))

    sim.run(
        # mp.at_every(1, mp.output_png(mp.Ey, "-Zc /home/draco/miniconda3/envs/mp/share/h5utils/colormaps/dkbluered")),
        until=30)

    # flux0 = mp.get_fluxes(flux_box)
    # flux_data = sim.get_flux_data(flux_box)
    # freq = mp.get_flux_freqs(flux_box)

    incident = mp.get_fluxes(incident_line)
    incident_data = sim.get_flux_data(incident_line)

    sim.reset_meep()

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=mp.Vector3(),
        geometry=geometry,
    )

    # flux_box = sim.add_flux(frq_cen, dfrq, nfrq,
    #     mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, 2 * r), weight = -1),
    #     mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, 2 * r)),
    #     mp.FluxRegion(center = mp.Vector3(y = -r), size = mp.Vector3(2 * r, 0), weight = -1),
    #     mp.FluxRegion(center = mp.Vector3(y = +r), size = mp.Vector3(2 * r, 0)))

    # flux_box1 = sim.add_flux(frq_cen, dfrq, nfrq,
    #     mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, 2 * r), weight = -1),
    #     mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, 2 * r)),
    #     mp.FluxRegion(center = mp.Vector3(y = -r), size = mp.Vector3(2 * r, 0), weight = -1),
    #     mp.FluxRegion(center = mp.Vector3(y = +r), size = mp.Vector3(2 * r, 0)))

    tran = sim.add_flux(frq_cen, dfrq, nfrq,
        mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, s)))
    refl = sim.add_flux(frq_cen, dfrq, nfrq,
        mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, s), weight = -1.0))

    sim.load_minus_flux_data(refl, incident_data)

    # sim.use_output_directory("spectrumtest")
    sim.run(
        # mp.at_beginning(mp.output_epsilon),
        # mp.at_every(1, mp.output_png(mp.Ey, "-Zc /home/draco/miniconda3/envs/mp/share/h5utils/colormaps/dkbluered -C $EPS")),
        until=30)

    tran1 = mp.get_fluxes(tran)
    refl1 = mp.get_fluxes(refl)
    freq = mp.get_flux_freqs(tran)
    
    incident = np.asarray(incident)
    tran1 = np.asarray(tran1)
    refl1 = np.asarray(refl1)
    freq = np.asarray(freq)
    #funny maths

    plt.figure()
    plt.plot(freq, tran1, label = f"tran")
    plt.plot(freq, refl1, label = f"refl")
    plt.plot(freq, incident, label = f"incident")
    plt.grid(True, which="both", ls="-")
    plt.xlabel("Frequency")
    plt.ylabel("Flux")
    plt.legend(loc="upper right")
    plt.title(f"Fluxtest2_{w}.png")
    plt.savefig(f"Fluxtest2_{w}.png")
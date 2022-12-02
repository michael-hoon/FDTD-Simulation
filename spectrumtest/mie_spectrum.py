import matplotlib.pyplot as plt
import numpy as np
import meep as mp


#set up flux boxes
r      = 1.0 #disk radius
wvlmax = 0.7
wvlmin = 0.3
frqmin = 1 / wvlmax
frqmax = 1 / wvlmin
frq_cen = (frqmax + frqmin) / 2
dfrq = frqmax - frqmin
nfrq = 100

dpml = 1
dair = 4
s = 2*(dair + dpml + r);
cell_size = mp.Vector3(s, s);
resolution = 50
pml_layers = [mp.PML(thickness=dpml)]

#laser
las_wvl     = [0.40, 0.45, 0.50, 0.55, 0.60]
las_cutoff  = 5.0
las_width   = 0.7/2.0/np.sqrt(np.log(2))

stored = {}
for wvl in las_wvl:
    # is_integrated=True necessary for any planewave source extending into PML
    sources = [
        mp.Source(
            mp.GaussianSource(
                wavelength = wvl, 
                width = las_width,
                cutoff =las_cutoff, 
                is_integrated = True),
            center = mp.Vector3(-dair/2, 0), 
            size = mp.Vector3(0, s),
            component = mp.Ey,
        )]
    #create single lorentzian material sphere at center
    lorentz_suscep = [
        mp.LorentzianSusceptibility(frequency = 1/0.3646595, gamma = 1/5.42347, sigma = 2.37226548)
        # mp.LorentzianSusceptibility(frequency = 1/0.4, gamma = 1/5.42347, sigma = 2.37226548)
        ]
    material = mp.Medium(epsilon = 1.001, E_susceptibilities = lorentz_suscep)
    geometry = [
        mp.Sphere(material=material, center=mp.Vector3(), radius=r)
    ]
    #Set up boundary region to measure energy spectrum

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        k_point=mp.Vector3(),
        symmetries=[],
    )

    line_x1 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, 2 * r))
    )

    line_x2 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, 2 * r))
    )

    line_y1 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(y = -r), size = mp.Vector3(2 * r, 0))
    )

    line_y2 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(y = +r), size = mp.Vector3(2 * r, 0))
    )

    sim.run(until_after_sources=10)

    freqs  = np.asarray(mp.get_flux_freqs(line_x1));
    dats   = [sim.get_flux_data(line)
            for line in [line_x1,line_x2,line_y1,line_y2]]
    inputfield = mp.get_fluxes(line_x1);

    sim.reset_meep()

    sim = mp.Simulation(
        resolution=resolution,
        cell_size=cell_size,
        boundary_layers=pml_layers,
        sources=sources,
        geometry=geometry,
        k_point=mp.Vector3(),
        symmetries=[],
    )

    line_x1 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(x = -r), size = mp.Vector3(0, 2 * r, 0))
    )

    line_x2 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(x = +r), size = mp.Vector3(0, 2 * r, 0))
    )

    line_y1 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(y = -r), size = mp.Vector3(2 * r, 0, 0))
    )

    line_y2 = sim.add_flux(
        frq_cen,
        dfrq,
        nfrq,
        mp.FluxRegion(center = mp.Vector3(y = +r), size = mp.Vector3(2 * r, 0, 0))
    )
    lines = [line_x1,line_x2,line_y1,line_y2];
    for dat,line in zip(dats,lines):
        sim.load_minus_flux_data(line,dat);

    sim.run(until_after_sources=100)
    cs = [ 1.0, -1.0, 1.0, -1.0]
    end_fluxes = [
        np.asarray(mp.get_fluxes(line))*c
        for line,c in zip(lines,cs) ];
    sumflux = - sum(end_fluxes) / inputfield;
    stored[wvl] = sumflux
    # np.savez("meep_spectrum.npz", freqs=freqs, spec=sumflux);
    plt.figure()
    plt.plot(freqs, sumflux,ls='-',marker='o',label=f'{wvl * 1000}nm');
    plt.grid(True, which="both", ls='-')
    plt.xlabel("Frequency (micron$^{-1}$)")
    plt.ylabel("A.U.")
    plt.legend()
    plt.title('2D "Mie Scattering" of a Lorentzian Material Disk')
    # plt.show();
    plt.savefig(f"Spectrumtest2_{wvl}.png")

    sim.reset_meep()

plt.figure()
for wvl in stored:
    plt.plot(freqs, stored[wvl],ls='-',marker='o',label=f'{wvl * 1000}nm');
plt.grid(True, which="both", ls='-')
plt.xlabel("Frequency (micron$^{-1}$)")
plt.ylabel("A.U.")
plt.legend()
plt.title('2D "Mie Scattering" of a Lorentzian Material Disk')
plt.savefig(f"Spectrumtest2all.png")
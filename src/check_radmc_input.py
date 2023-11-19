"""
  check_radmc_input.py

  Purpose: Read radmc *.inp files.
           Make plots of input molecular abundances and gas temperatures

  Author: Stella Offner

"""

from radmc3dPy import *
import matplotlib.pylab as plb
import numpy as np
from radmc3dPy.image import radmc3dImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# Default input values (see inputs_gizmo_carver.py)
dust_to_gas = 0.01  # dust to gas ratio
hydrogen_ratio = 1.4  # nH/nHe = 10
helium_mass_fraction = 0.284
molecular_abundance = 1e-4  # 2*10**-8 # abundance of NH3 relative to H2
mh = 1.67e-24
pc = 3.09e18
box = 10.0 * pc
fac = 1  # Test / Prefactor for NH3 abundance
mname = "co"

# Plots to make
plot_mol = True
plot_temp = True
plot_dust = False
save = True  # Save plots

# Read input data used for RADMC calculation
data = analyze.readData(
    gdens=True, ddens=plot_dust, ispec=mname, gtemp=plot_temp, gvel=save, gmag=save
)
sz = len(data.ndens_mol)
dx = box / sz
print("Box (pc), sz, dx (pc)", box / pc, sz, dx / pc)
# print(data.__dict__) # Names of data structures

if save == True:
    np.save("GasVel.npy", data.gasvel)
    np.save("GasMagField.npy", data.gasmag)

# Calculate N(H2) and N(NH3)
if plot_mol:
    nh2 = data.ndens_mol / molecular_abundance  # Convert to nh2
    Nh2x = nh2.sum(0) * dx
    Nh2y = nh2.sum(1) * dx
    Nh2z = nh2.sum(2) * dx

    print(np.mean(nh2), np.median(nh2))
    print(np.mean(data.ndens_mol), np.median(data.ndens_mol))
    # print(np.shape(data.ndens_mol))
    print("Max density =", np.max(data.ndens_mol))

    # Plot Column Density distribution
    fig, ax = plt.subplots()
    ax.hist(
        np.ravel(Nh2x * molecular_abundance),
        bins=np.logspace(np.log10(1e10), np.log10(2e15), 25),
    )
    ax.set_xscale("log")
    ax.set(xlabel="N$_{" + mname + "}$ [cm$^{-2}$] ", ylabel=" N")
    if save == True:
        plt.savefig("Col" + mname + "_dist_check.png")
        np.save("Ndens_" + mname + ".npy", data.ndens_mol)
    else:
        plt.show()

    # Plot effective NH3 abundance (as compared to dust, where we assume the dust traces the molecular gas)
    if plot_dust == True:
        rho = data.rhodust / dust_to_gas
        nh2eff = rho / (2.8 * 1.67e-24)  # Effective H2 assuming all molecular
        fig, ax = plt.subplots()
        ax.hist(
            np.ravel(nh2 * molecular_abundance * fac / nh2eff),
            bins=np.logspace(np.log10(1e-10), np.log10(2e-7), 25),
        )
        ax.set_xscale("log")
        ax.set(xlabel="X$_{" + mname + "}$ [cm$^2$/g] ", ylabel=" N")
        if save == True:
            plt.savefig("X_" + mname + "_dist_check.png")
        else:
            plt.show()

    # Plot Column Density maps
    fig = plt.figure()
    gs = fig.add_gridspec(1, 3, hspace=0, wspace=0)
    ax = gs.subplots(sharex=True, sharey=True)

    ax[0].imshow(np.log10(Nh2x[:, :, 0]), extent=[0, box / pc, 0, box / pc])
    ax[1].imshow(np.log10(Nh2y[:, :, 0]), extent=[0, box / pc, 0, box / pc])
    ax[0].set(xlabel="L [pc] ", ylabel="L [pc]")
    ax[2].imshow(np.log10(Nh2z[:, :, 0]), extent=[0, box / pc, 0, box / pc])
    ax[0].set(xlabel="L [pc] ", ylabel="L [pc]")
    ax[1].set(title="H$_2$ Column Density")

    if save == True:
        plt.savefig("ColH2_dist_check.png")
    else:
        plt.show()
    print("Max Nh2", np.max(Nh2x), np.max(Nh2y), np.max(Nh2z))


if plot_temp:
    # Plot gas (and dust) temperature distribution

    fig, ax = plt.subplots()
    ax.hist(np.ravel(data.gastemp), bins=np.logspace(np.log10(1e0), np.log10(1e5), 20))
    ax.set_xscale("log")
    ax.set(xlabel="T [K]", ylabel=" N")
    if save == True:
        plt.savefig("Gastemp_dist_check.png")
        np.save("Gastemp.npy", data.gastemp)
    else:
        plt.show()

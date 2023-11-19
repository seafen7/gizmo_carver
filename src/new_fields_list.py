import yt
import inputs_gizmo_carver as inputs
from yt.units import *
import numpy as np


# Definition of the dust density. Uses dust to gas ratio from inputs file
def _DustDensity(field, data):
    return inputs.dust_to_gas * data[("PartType0", "Density")]


# Definition of the gas temperature. Uses info from inputs
def _gas_temp(field, data):
    try:
        return yt.YTArray(data[("PartType0", "Temperature")], "K")
    except:
        y_helium = inputs.helium_mass_fraction / (
            4.0 * (1 - inputs.helium_mass_fraction)
        )
        mu = (1 + 4.0 * y_helium) / (1 + y_helium)

        # In molecular gas electron abundance < 10^-5
        # mu = (1+4.0*y_helium)/(1+y_helium+data[(‘PartType0’,‘ElectronAbundance’)]
        # const = ((1.2*mh.in_mks()*(gamma-1))/kboltz.in_mks())

        const = (mu * mh.in_mks() * (inputs.gamma - 1)) / (1e6 * kboltz.in_mks())
        # Note the LTE temperature tables are limited to 1e5K
        # Non-LTE may require lower temperatures
        return (
            data[("PartType0", "InternalEnergy")]
            * const
            * (data[("PartType0", "InternalEnergy")] < yt.YTArray([300], "K") / const)
        )


# Definition of the target species field. Uses info from inputs
def _H2NumDensity(field, data):
    return (
        data[("PartType0", "Density")]
        * data[("PartType0", "MolecularMassFraction")]
        * data[("PartType0", "NeutralHydrogenAbundance")]
        * (1 - inputs.helium_mass_fraction)
        / (inputs.mol_hydrogen_ratio * mh)
    )


# Definition of the target species field. Uses info from inputs
def _MolecularNumDensity(field, data):
    return (
        data[("PartType0", "H2NumDensity")]
        * inputs.molecular_abundance
        * (data[("PartType0", "H2NumDensity")] > yt.YTArray([1e2], "cm**-3"))
        * (data[("PartType0", "gas_temperature")] < yt.YTArray([1e2], "K"))
    )


# Definition of the microturbulence at each point. Uses info from inputs
def _MicroTurb(field, data):
    turb = data["velocity_x"]
    turb[turb >= 0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")
    turb[turb <= 0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")

    return turb


# Definition of the dust temperature. Assumes same as gas temperature
def _DustTemperature(field, data):
    try:
        return yt.YTArray(data[("PartType0", "Dust_Temperature")], "K")
    except:
        return data["gas_temperature"].clip(0,20) # rough approximation - can do better using radiation fields if required


# Create a mask based on accreted particles
# Need to update if want to exclude particles not yet formed
def _Mask(field, data):
    mask = data[("PartType0", "Density")] * 0.0
    ids = np.unique(alllines[:, 1])
    # This is the slow way. See fast way below
    for id in ids:
        acclist = np.where(alllines[:, 1] == id)[0]  # All rows accreted by this star
        # Particle ids of the accreted particles
        ind = np.where(
            np.isin(data[("PartType0", "ParticleIDs")], alllines[acclist, 6]) == True
        )[0]
        mask[ind] = 1
    return mask


# Add all the fields

# yt.add_field(("PartType0", "gas_temperature"), function=_gas_temp, units="K", sampling_type='particle', force_override=True)

# yt.add_field(("PartType0", "dust_temperature"), function=_DustTemperature, units="K", sampling_type='particle', force_override=True)

yt.add_field(
    ("PartType0", "DustDensity"),
    function=_DustDensity,
    units="g/cm**3",
    sampling_type="particle",
    force_override=True,
)

yt.add_field(
    ("PartType0", "H2NumDensity"),
    function=_H2NumDensity,
    units="cm**-3",
    sampling_type="particle",
    force_override=True,
)

yt.add_field(
    ("PartType0", "MolecularNumDensity"),
    function=_MolecularNumDensity,
    units="cm**-3",
    sampling_type="particle",
    force_override=True,
)

yt.add_field(
    ("PartType0", "microturbulence_speed"),
    function=_MicroTurb,
    units="cm/s",
    sampling_type="particle",
    force_override=True,
)

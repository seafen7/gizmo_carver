"""
   main_gizmo_carver.py

   Purpose:
        Driver file for NH3-based RADMC carve routines. Call this file when you want
        to run the code. Should not ever need to edit this file.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022
        
        Modified from: main_CarveOut.py, written by:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.9, yt 4.0.2
"""

import yt
from writer_gizmo_carver import RadMC3DWriter_Gizmo
from globals_gizmo_carver import *
import inputs_gizmo_carver as inputs
from yt.units import *
import os
from datetime import datetime
import shutil
import numpy as np

# Definition of the dust field. Uses dust to gas ratio from inputs file
def _DustDensity(field, data):
    return inputs.dust_to_gas * data[('PartType0', 'Density')]
yt.add_field(('PartType0', 'DustDensity'), function=_DustDensity, units='g/cm**3', sampling_type='particle', force_override=True)

# Definition of the target species field. Uses info from inputs
def _MolecularNumDensity(field, data):
    return data[('PartType0', 'Density')]*data[('PartType0', 'MolecularMassFraction')]/(inputs.hydrogen_ratio*mh)*inputs.molecular_abundance
yt.add_field(('PartType0', 'MolecularNumDensity'), function=_MolecularNumDensity, units='cm**-3', sampling_type='particle', force_override=True)

# Definition of the microturbulence at each point. Uses info from inputs
def _MicroTurb(field, data):
    turb = data["velocity_x"]
    turb[turb>=0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")
    turb[turb<=0] = yt.YTQuantity(inputs.microturbulence_speed, "cm/s")
    #print(turb.min())
    #print(turb.max())
    return turb
yt.add_field(("PartType0", "microturbulence_speed"), function=_MicroTurb, units="cm/s", sampling_type='particle', force_override=True)

# Definition of the gas temperature. Uses info from inputs
def _gas_temp(field, data):
    y_helium = inputs.helium_mass_fraction/(4.0*(1-inputs.helium_mass_fraction))
    mu = (1+4.0*y_helium)/(1+y_helium) # Here electron abundance < 10^-5
    # print(‘const =’, mu, gamma, mh, kboltz, (mu*mh*(gamma-1))/kboltz)
    #mu = (1+4.0*y_helium)/(1+y_helium+data[(‘PartType0’,‘ElectronAbundance’)]
    #const = ((1.2*mh.in_mks()*(gamma-1))/kboltz.in_mks())
    return (data[('PartType0','InternalEnergy')]*((mu*mh.in_mks()*(inputs.gamma-1))/(1e6*kboltz.in_mks())))
yt.add_field(("PartType0", "gas_temperature"), function=_gas_temp, units="K", sampling_type='particle', force_override=True)

# Definition of the dust temperature. Assumes same as gas for now
def _DustTemperature(field, data):
    return data["gas_temperature"]
yt.add_field(("PartType0", "dust_temperature"), function=_DustTemperature, units="K", sampling_type='particle', force_override=True)


# Loads file into YT
ds = yt.load(inputs.hdf5_file, unit_base=inputs.unit_base)
try:
    print("Loaded file " + str(ds)) # using classic print() for try/except to work
except NameError:
    assert False, "Unable to properly load file into YT!"


now = datetime.now()
dt_string = now.strftime("%m.%d.%y_%H'%M'%S")

# Create working directory for this run
current_dir = inputs.output_filepath
working_dir_name = os.path.join(current_dir, 'RADMC_inputs_' + dt_string)
os.mkdir(working_dir_name)

box_left = np.add(inputs.box_center, -inputs.box_size)
box_right = np.add(inputs.box_center, inputs.box_size)
box_left_cgs = [Convert(x, inputs.box_units, 'cm', 'cm') for x in box_left]
box_right_cgs = [Convert(x, inputs.box_units, 'cm', 'cm') for x in box_right]

box_left_cgs = unyt_array(box_left, 'cm')
box_right_cgs = unyt_array(box_right, 'cm')

print("\nCarving between Left = " + str(unyt_array(box_left, inputs.box_units)))
print("            to Right = " + str(unyt_array(box_right, inputs.box_units)))
print("       w/ Resolution = " + str(abs(inputs.box_dim)) + " x " + str(abs(inputs.box_dim)) + "\n")

writer = RadMC3DWriter_Gizmo(ds, a_boxLeft=box_left_cgs, a_boxRight=box_right_cgs, a_boxDim=inputs.box_dim)

# Write the amr grid file (fast)
print("1/7: Writing amr grid file (fast!)")
writer.write_amr_grid(os.path.join(working_dir_name, inputs.out_afname))

# Write the number density file for species (slow)
print("2/7: Writing number density file (slow.)")
writer.write_line_file(('PartType0', 'AmmoniaNumDensity'), os.path.join(working_dir_name, inputs.out_nfname))

# Write the dust density file for dust (slow)
print("3/7: Writing dust density file (slow.)")
writer.write_dust_file(('PartType0', 'DustDensity'), os.path.join(working_dir_name, inputs.out_ddfname))

print("4.5/7: Writing microturbulence file (slow.)")
writer.write_line_file(("PartType0", "microturbulence_speed"), os.path.join(working_dir_name, inputs.out_mtfname))

# Write the temperature file for species or dust (slow)
print("5/7: Writing temperature file (slower..)")
writer.write_line_file(("PartType0", "gas_temperature"), os.path.join(working_dir_name, inputs.out_tfname))

# Assuming dust temperature is same as gas for now...
print("6/7: Writing dust temperature file (slower..)")
writer.write_dust_file(("PartType0", "dust_temperature"), os.path.join(working_dir_name, inputs.out_dtfname))

# Write the gas velocity file (slow x 3)
print("7/7: Writing velocity file (slowest...)")
velocity_fields = ["velocity_x", "velocity_y", "velocity_z"]
writer.write_line_file(velocity_fields, os.path.join(working_dir_name, inputs.out_vfname))

# Copy over existing files
print('Copying default files...')
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_molname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_wlmname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_cwlname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dksname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_dtpname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_linname), working_dir_name)
shutil.copy(os.path.join(inputs.existing_filepath, inputs.out_rmcname), working_dir_name)

print('Done! Output files generated at: \n\n' + working_dir_name)

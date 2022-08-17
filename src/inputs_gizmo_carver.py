"""
   inputs_gizmo_carver.py

   Purpose:
        Input file for RadMC carving routines. This is the only file you should
        edit. Follow the comments below to see what each variable represents.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022
        
        Modified from:
        inputs_CarveOut.py, written by
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.9, yt 4.0.2
"""

from yt.units import * 
import numpy as np

# Constants for calculating derived fields
dust_to_gas = 0.01
mol_hydrogen_ratio = 2.0
microturbulence_speed = 1e4 # cgs
gamma = 5.0/3.0 # Note gamma is not constant and this function is an approximation.
helium_mass_fraction = 0.284 # Default mass fraction in Gizmo

# number fraction of target species
molecular_abundance = 2*10**-8 # abundance of NH3 relative to H2

# Mask abundance based on accreted particles
mask_abundance = False

# Units of the below box values ('pc','cm','AU','ly' accepted)
box_units = 'pc'

# x, y, z coordinates for the center of the carved domain (e.g., location of a star core)
# The values should match the unit given by box_units
#box_center = [15.95957649, 15.54566532, 15.19446488] #M2e3, center = 15
box_center = [55.8,44.8, 57.2] # Core 1, 330  #M2e4 center = 50

# Routine will generate input files for a square area centered at box_center 
# extending to box_center += box_size on each side
# Use same units as box_units
box_size = 0.5 # pc (=L/2)

# Resolution of the resulting image (give as a complex number, e.g. for 
# box_dim = 64j, the resulting image will be 64x64)
box_dim = 256j

# Snapshot number
snap = '330'

# Name tag for output file directory
tag = 'sn'+snap+'_'+ np.str(np.int(np.imag(box_dim)))+'_'

# Filepath of the HDF5 file name to read in
# If the HDF5 file is located in the same directory as the script files, 
# you can just put the file name
hdf5_dir = '/scratch3/03532/mgrudic/STARFORGE_RT/production/M2e4_R10_S0_T1_B0.01_Res271_n2_sol0.5_42/output/'

# unit base to use for calculations
unit_base = {'UnitMagneticField_in_gauss':  1e+4,
             'UnitLength_in_cm'         : 3.08568e+18,
             'UnitMass_in_g'            :   1.989e+33,
             'UnitVelocity_in_cm_per_s' :      100}

# Filepath for directory containing input files that are not generated by the carver routine.
# These files are still necessary for running RADMC-3D. The files are:
# camera_wavelength_micron.inp
# dustkappa_silicate.inp
# dustopac.inp
# lines.inp
# molecule_nh3.inp (Or data file for other target species)
# radmc3d.inp
# wavelength_micron.inp
existing_filepath = '/home1/00653/tg458122/gizmo_carver/default_files'

# Filepath for storing output files. Routine will make a working directory within this
# output directory for each run.
output_filepath = '/work2/00653/tg458122/frontera/_gizmo_radmc/M2e4_fid_output_files/_cores'

#M2e3_R3_S0_T1_B0.01_Res126_n2_sol0.5_42/'
#/M2e4_fid_output_files/'#'./output_files' #Used for M2e3 tests

# Write a new line file rather than use defaults file
write_line_file = True
# velocity max/min for the wavelength file
vmax = 10 # km/s
# number of wavelengths in output file
nwav = 256
# Line rest frequency, see molecule_x.inp file
restfreq = 23.69449550E9 # 12CO (1-0): 115.2712018E9, (2-1) 230.5380000E9
#restfreq = 23.69449550E9 # NH3[1,1]

# Output file names for use in RADMC3D
out_afname = "amr_grid.inp"       # output file name for amr grid
out_nfname = "numberdens_nh3.inp" # output file name for target species above
out_nhname = "numberdens_h2.inp" # output file name for h2 number density
out_vfname = "gas_velocity.inp"   # output file name for velocity
out_tfname = "gas_temperature.inp"    # output file name for temperature
out_ddfname = "dust_density.inp" # output file name for dust density
out_dtfname = "dust_temperature.dat" # output for dust temperature (requires .dat)
out_mtfname = "microturbulence.inp" # output for microturbulence

# Names of existing files
out_molname = 'molecule_nh3.inp' # (Or data file for other target species)
out_wlmname = 'wavelength_micron.inp'
out_cwlname = 'camera_wavelength_micron.inp'
out_dksname = 'dustkappa_silicate.inp'
out_dtpname = 'dustopac.inp'
out_linname = 'lines.inp'
out_rmcname = 'radmc3d.inp'
out_execute = 'radmc3d'
out_subscript = 'submit_radmc.sh'
out_makeinput = 'input_info.txt' # Save the setup parameters and output

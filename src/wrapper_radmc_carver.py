"""
   wrapper_radmc_carver.py

   Purpose:
        Driver file for RADMC carve routines. Call this file when you want
        to run the code. Should not ever need to edit this file.

   Before running:
        1. Check variables defined in inputs_gizmo_carver.py
        2. Check submit_radmc.sh and consistency with (1) -- sizepc = box*2

   Author:
        Stella Offner

   Written/Tested with Python 3.9, yt 4.0.2
"""
from main_gizmo_carver import *
import subprocess
import os as os
import glob as glob

# Create file with all the batch submit commands
f = open("submit_all.sh", "w")

# Setup RADMC files
output_dir = main_gizmo_carver()  # Call with inputs_gizmo_carver defaults
print("Finished wrapper...", output_dir)

# Change the current working directory
os.chdir("/" + output_dir)

# Check result
subprocess.call(
    "python /home1/00653/tg458122/gizmo_carver/src/check_radmc_input.py", shell=True
)

# Change the current working directory
os.chdir("../../")

# Add the submission command to a script
f.write("sbatch %s/submit_radmc.sh" % output_dir)

f.close()

"""
   globals_gizmo_carver.py

   Purpose:
        Constants and simple functions for use in the RadMC carving routines.
        Should not need to edit this file.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022
        
        Modified from:
        globals_CarveOut.py, written by:
        Aaron T. Lee, aaron.t.lee@utexas.edu
        Spring 2018

   Written/Tested with Python 3.9, yt 4.0.2
"""

# Unit conversions
class UnitConv:
    cgsunits = ['cm','g','sec']
    lenunits = ['cm','au','ly','pc']
    massunits = ['sol','g']
    timeunits = ['sec','yr','Myr']
    au2cm = 1.496e13 # AU to cm
    pc2cm = 3.08567758128e18 # parsec to cm
    ly2cm = 9.463e18 # light-year to cm
    sol2g = 1.989e33 # solar mass to g
    sol2ergs = 3.83e33 # solar luminosity to erg/s
    sol2cm = 6.9550e10 # solar radius to cm
    yr2sec = 3.154e7   # year in seconds
    Myr2sec = 3.154e13 # Million years in seconds
    
# Unit conversion function
# Inputs: x = numerical value,
#         unit_in = unit of input x, unit_out = desired out unit
#         cgsunit = fundamental cgs unit (e.g., 'g' if mass, 'cm' if length)
def Convert(x,unit_in,unit_out,cgsunit):
    if cgsunit not in UnitConv.cgsunits:
        print("Unit conversion failure for cgs unit: " + str(cgsunit))

    if(cgsunit=='cm'):
        possible_units = UnitConv.lenunits
    elif(cgsunit=='g'):
        possible_units = UnitConv.massunits
    else:
        possible_units = UnitConv.timeunits

    units = [ unit_in.lower(), unit_out.lower()]
    if False in [x in possible_units for x in units]:
        print("Unit conversion failure for units: " + str(unit_in) + " and " + str(unit_out))

    # Convert first to CGS
    val1 = 1.0
    if units[0] not in UnitConv.cgsunits :
        val1 =eval( 'UnitConv.' + units[0] + '2' + cgsunit )

    # Convert to what you want
    val2 = 1.0
    if units[1] not in UnitConv.cgsunits :
        val2 = 1.0/eval( 'UnitConv.' + units[1] + '2' + cgsunit )

    return(x*val1*val2)
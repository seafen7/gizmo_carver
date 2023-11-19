#!/usr/bin/env python

import numpy
import math
import astropy.io.fits as pyfits
import warnings
import os
import sys

from astropy import coordinates as coord
from astropy.coordinates import FK5
from astropy import units as u

##################################################################################
# User-defined input variables

inputFileNameBase = "image"
outputFileNameBase = inputFileNameBase  # freq or_vel will be added

# intensity_unit = [0, 1, 2, 3]
intensity_unit = eval(sys.argv[2])
if isinstance(intensity_unit, int):
    intensity_unit = [intensity_unit]
#   0 = erg/cm**2/s/Hz/ster
#   1 = Jy/ster
#   2 = Jy/pixel
#   3 = K/ster -> K
# e.g. intensity_unit = [0, 2] returns two fits files

spectral_axis_flag = 0  # flag; sets the spectral axis unit in the fits file
#   0 = velocity
#   1 = frequency

if spectral_axis_flag:
    outputFileNameBase += "_freq"
else:
    outputFileNameBase += "_vel"


getBandwidthFromFile = 1  # flag; get bandwidth from file or from input keyword
#   0 = use input value
#   1 = use value determined from input file

# if getBandwidthFromFile flag is set to 1, this keyword is ignored.
delta_v_kms = 40.0 / 201.0

kB = 1.3807e-16  # Bolzmann's constant	 [erg/K]
c_ms = 2.99792458e8  # speed of light in m/s
velocityCenter_kms = 0.0  # center velocity [km/s]

# CO1-0
# restwavelength_mum = 2600.7575643
restwavelength_mum = numpy.float(sys.argv[1])

restfrequency_ghz = c_ms / (restwavelength_mum * 1e3)  # rest frequency [GHz]
print("Rest frequency [GHz]:", restfrequency_ghz)


#
centerCoordinates = "00:00:00.000  -00:00:00.00"  # assuming FK5
distance_pc = 1000.0  # distance from observer to the source [pc]

# General user input variables - they might not be changed very often

project = "SILCC zoomin: superfast zoomin"
software = "RADMC-3D (v0.41)"


# ================================================================================
#
# Starting the fits creation and header edition
#
# ================================================================================

counter = 0
intensity_ergSter = []
wavelength_micron = []

numberChannels = 0
numberPixels_x = 0
numberPixels_y = 0


# check if the input file exists, if not, abort

if os.path.isfile("%s.out" % inputFileNameBase) is True:
    print("")
    print("Reading in file: %s.out" % inputFileNameBase)
    print("")

    with open("%s.out" % inputFileNameBase) as image:
        counterw = 0
        counteri = 0
        for line in image:
            counter += 1
            line = line.rstrip("\n")
            # print status report to screen

            if counter % 1000000 == 0:
                numberLines = numberPixels_x * numberPixels_y * numberChannels
                print("	Reading in line   %s / %s" % (counter, numberLines))

            #
            if counter == 2:
                numberPixels_x = int(line.split()[0])
                numberPixels_y = int(line.split()[1])

                print("	Number of pixels in x-direction: ", numberPixels_x)
                print("	Number of pixels in y-direction: ", numberPixels_y)
                print("")

            #
            elif counter == 3:
                numberChannels = int(line.split()[0])
                wavelength_micron = numpy.zeros(numberChannels)

                print("	Number of channels: ", numberChannels)
                print("")

            #
            elif counter == 4:
                pixelSize_x_cm = float(line.split()[0])
                pixelSize_y_cm = float(line.split()[1])
                intensity_ergSter = numpy.zeros(
                    (numberChannels, numberPixels_y, numberPixels_x)
                )

                print("	Pixel size in x direction: %.4e cm" % pixelSize_x_cm)
                print("	Pixel size in y direction: %.4e cm" % pixelSize_y_cm)
                print("")

            #
            # get the wavelengths
            elif (
                (counter >= 5)
                and (counter < (5 + numberChannels))
                and (len(line.split()) > 0)
            ):
                wavelength_micron[counterw] = float(line)
                counterw += 1

            #
            # get the intensities

            elif (counter > (5 + numberChannels)) and (len(line.split()) > 0):
                intensity_ergSter[
                    numpy.int(counteri / numberPixels_x / numberPixels_y),
                    numpy.int(counteri / numberPixels_x) % numberPixels_y,
                    counteri % numberPixels_x,
                ] = line
                counteri += 1

    # defining some constants
    au_cm = 1.49597870700e13  # astronomical unit
    c_ms = 2.99792458e8  # speed of light

    # conversion factor from erg/cm**2/s/Hz to Jansky
    erg2jy = 1.0e23

    # conversion factor from arcseconds to degree, degree to radian
    arcsec2degree = 1.0 / 3600.0
    degree2radian = math.pi / 180.0

    # small angle approximation: tan(theta) ~ theta
    # definition of a parsec: 1 pc = 1 au / 1"
    # --> theta in arcsec = pixelSize in au / distance to the source in pc
    pixelSize_x_au = pixelSize_x_cm / au_cm
    pixelSize_y_au = pixelSize_y_cm / au_cm

    pixelSize_x_degree = pixelSize_x_au / distance_pc * arcsec2degree
    pixelSize_y_degree = pixelSize_y_au / distance_pc * arcsec2degree

    pixelSize_x_radian = pixelSize_x_degree * degree2radian
    pixelSize_y_radian = pixelSize_y_degree * degree2radian

    # conversion factor for steradian to pixel = area of a pixel in steradian
    ster2pixel = pixelSize_x_radian * pixelSize_y_radian

    # restfrequency in Hz
    restfrequency_hz = restfrequency_ghz * 1.0e9

    # calculate frequency resolution using the radio convention
    if getBandwidthFromFile == 0:
        delta_v_ms = -1.0 * delta_v_kms * 1.0e3
        delta_nu_hz = -1.0 * restfrequency_hz * delta_v_ms / c_ms

        if spectral_axis_flag == 0:
            print("")
            print("	Using a velocity bandwidth of %.10f m/s" % delta_v_ms)

        elif spectral_axis_flag == 1:
            print("")
            print("	Using a calculated frequency bandwidth of %.10f Hz" % delta_nu_hz)

    elif getBandwidthFromFile == 1:
        # calculating the bandwidth from the input file
        wavelength_m = wavelength_micron * 1.0e-6
        wavelength_cm = wavelength_micron * 1.0e-4
        nu_hz = c_ms / wavelength_m

        Delta_v_ms = []
        if numberChannels > 1:
            for i in range(len(wavelength_m) - 1):
                Delta_v_ms.append(c_ms * (nu_hz[i + 1] - nu_hz[i]) / restfrequency_hz)
        else:
            Delta_v_ms.append(0.0)

        delta_v_ms = numpy.mean(Delta_v_ms)
        delta_nu_hz = -1.0 * restfrequency_hz * delta_v_ms / c_ms

        if spectral_axis_flag == 0:
            print("")
            print(
                "	Calculated a velocity bandwidth of %.10f m/s from input file"
                % delta_v_ms
            )

        elif spectral_axis_flag == 1:
            print("")
            print(
                "	Calculated a frequency bandwidth of %.10f Hz from input file"
                % delta_nu_hz
            )

    # wavelength_cm = numpy.asarray(wavelength_m) * 100.	# cm: 1.e-2, mum: 1.e-6

    # create the different fits files (different intensity units, different spectral axis units)

    for intensity_unit_flag in intensity_unit:
        if intensity_unit_flag == 0:
            # character string describing the physical units in which the quantities in the array are expressed
            bunit = "erg/cm^2/s/Hz/ster"
            fileEnding = "ergSter.fits"

        elif intensity_unit_flag == 1:
            bunit = "Jy/ster"
            fileEnding = "jySter.fits"

        elif intensity_unit_flag == 2:
            bunit = "Jy/pixel"
            fileEnding = "jyPixel.fits"

        elif intensity_unit_flag == 3:
            bunit = "K"  # K/ster
            fileEnding = "KSter.fits"

        ###############################################

        print("")
        print("	Preparing fits file data array")

        if intensity_unit_flag == 0:
            hdu = pyfits.PrimaryHDU(numpy.array(intensity_ergSter))

        elif intensity_unit_flag == 1:
            hdu = pyfits.PrimaryHDU(numpy.array(intensity_ergSter * erg2jy))

        elif intensity_unit_flag == 2:
            hdu = pyfits.PrimaryHDU(
                numpy.array(intensity_ergSter * erg2jy * ster2pixel)
            )

        elif intensity_unit_flag == 3:
            hdu = pyfits.PrimaryHDU(
                numpy.array(intensity_ergSter.T * wavelength_cm**2).T / (2.0 * kB)
            )

        #
        # get the coordinates of the reference position

        refpos = coord.SkyCoord(centerCoordinates, frame=FK5, unit=(u.hour, u.degree))

        #
        # convert the reference position coordinates to degree

        ra_refpos_degree = refpos.ra.degree
        dec_refpos_degree = refpos.dec.degree

        #
        # General fits header information

        hdu.header["bitpix"] = -32

        hdu.header["software"] = software
        hdu.header["project"] = project
        hdu.header["observer"] = "ag-walch"

        hdu.header[
            "equinox"
        ] = 2000.0  # epoch keyword is deprecated, equinox should be used instead
        hdu.header["restfrq"] = restfrequency_hz  # physical unit must be Hertz!

        #
        # first axis (right-ascension)

        hdu.header[
            "ctype1"
        ] = "RA---SIN"  # type of the coodinate system; SIN = Slant orthographic
        hdu.header[
            "cdelt1"
        ] = (
            -pixelSize_x_degree
        )  # x coordinate increment; minus sign, since RA axis points from right to left
        hdu.header["crval1"] = ra_refpos_degree  # world coordinate at reference point
        hdu.header["crpix1"] = (
            numberPixels_x + 1.0
        ) / 2.0  # reference pixel x coordinate
        hdu.header["cunit1"] = "degree"

        #
        # second axis (declination)

        hdu.header[
            "ctype2"
        ] = "DEC--SIN"  # type of the coodinate system; SIN = Slant orthographic
        hdu.header["cdelt2"] = pixelSize_y_degree  # y coordinate increment
        hdu.header["crval2"] = dec_refpos_degree  # world coordinate at reference point
        hdu.header["crpix2"] = (
            numberPixels_y + 1.0
        ) / 2.0  # reference pixel y coordinate
        hdu.header["cunit2"] = "degree"

        #
        # third axis (spectral)

        if spectral_axis_flag == 0:
            hdu.header["ctype3"] = "VELO-LSR"
            hdu.header["cdelt3"] = delta_v_ms
            hdu.header["crval3"] = velocityCenter_kms * 1.0e3
            hdu.header["crpix3"] = (numberChannels + 1.0) / 2.0
            hdu.header["cunit3"] = "m/s"

        elif spectral_axis_flag == 1:
            hdu.header["ctype3"] = "FREQ"
            hdu.header["cdelt3"] = delta_nu_hz
            hdu.header["crval3"] = restfrequency_hz * (
                1 - velocityCenter_kms * 1.0e3 / c_ms
            )
            hdu.header["crpix3"] = numberChannels / 2.0
            hdu.header["cunit3"] = "Hz"

        hdu.header["bunit"] = bunit

        # fourth axis (stokes; not present here)

        ##############################################################################
        #

        print("	Writing out: %s_%s" % (outputFileNameBase, fileEnding))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            pyfits.writeto(
                "%s_%s" % (outputFileNameBase, fileEnding),
                numpy.float32(hdu.data),
                hdu.header,
                clobber=True,
            )

    print("DONE")

else:
    print("Input file %s.out not found" % inputFileNameBase)
    print("Aborting")

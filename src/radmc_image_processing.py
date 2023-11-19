"""
   radmc_image_processing.py

   Purpose:
        Contains functions for plotting RADMC-3D output image files using MatPlotLib.
        Edit this file to create custom plotting routines.

   Author:
        Sean Feng, feng.sean01@utexas.edu
        Spring 2022

   Written/Tested with Python 3.9, radmc3dpy 0.30.2
"""

from radmc3dPy import *
import matplotlib.pylab as plb
import numpy as np
from radmc3dPy.image import radmc3dImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

freq = 23.69449550  # NH3 transition in GHz
nu = freq * 1e9  # [Hz]
# wav = (c/nu)*1e4   # [micron]
dist = 300.0  # Region distance [pc]

plt.close("all")


def plotMomentMaps(
    moments=[0, 1, 2], imgPath=None, savePath=None, vclip=None, convolve=False
):
    """
    Parameters
    ----------

    moments: list of moments to be plotted (0, 1, 2)

    imgPath: name of directory containing image.out
           If not specified, function will look in current working directory.

    savePath: name of directory to save copy of moment map plot to
           Saves image to directory specified by imgPath and savePath (if any)

     vclip: [[min,max]] for each moment map

     Uses name of image directory to name image file
     Give path names as raw strings

    """

    for i in moments:
        plt.figure()
        m_image = radmc3dImage()

        if imgPath is not None:
            m_image.readImage(os.path.join(imgPath, "image.out"))
        else:
            m_image.readImage()

        if convolve:
            m_image.imConv(fwhm=[10, 10], pa=0, tdiam_prim=None, tdiam_sec=None)

        m_image.plotMomentMap(moment=i, nu0=nu, dpc=1, Tb=True, vclip=vclip[i])

        plt.title("Moment " + str(i) + " Map")

        if imgPath is not None:
            plt.savefig(
                os.path.join(
                    imgPath, os.path.basename(imgPath) + "moment_" + str(i) + ".png"
                )
            )
        if savePath is not None:
            plt.savefig(
                os.path.join(
                    savePath, os.path.basename(imgPath) + "moment_" + str(i) + ".png"
                )
            )


vclip = [[0.1, 20], [-2, 2], [0.1, 5]]
plotMomentMaps(moments=[0, 1, 2], imgPath=r"./", vclip=vclip)

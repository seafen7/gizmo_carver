from radmc3dPy import *
import matplotlib.pylab as plb
import numpy as np
from radmc3dPy.image import radmc3dImage
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

freq = 23.69449550
nu = freq*1e9 # [Hz]
wav = (3e10/nu)*1e4 # [micron]

plt.close('all')


# PARAMETERS
# moments: list of moments to be plotted
# imgPath: name of directory containing image.out (if not specified, function will look in current working directory)
# savePath: name of directory to save copy of moment map plot to
#
# Saves image to directory specified by imgPath and savePath (if any)
# Uses name of image directory to name image file
# Give path names as raw strings
def plotMomentMaps(moments=[0, 1, 2], imgPath=None, savePath=None):
    for i in moments:
        plt.figure()
        m_image = radmc3dImage()
        
        if imgPath is not None:
            m_image.readImage(os.path.join(imgPath, 'image.out'))
        else:
            m_image.readImage()
            
        m_image.plotMomentMap(moment=i, nu0=nu, wav0=wav)
        plt.title('moment ' + str(i) + ' map')
               
        if imgPath is not None:
            plt.savefig(os.path.join(imgPath, os.path.basename(imgPath) + '_moment ' + str(i) + '.png'))
        if savePath is not None:
            plt.savefig(os.path.join(savePath, os.path.basename(imgPath) + '_moment ' + str(i) + '.png'))

# plotMomentMaps(moments=[0], imgPath=(r"C:\Users\seany\Documents\Offner Research\RADMC_test\RADMC_inputs 05.12.22_21'54'18"),
#                savePath = r"C:\Users\seany\Documents\Offner Research\RADMC_test\StarParticles")

# plotMomentMaps(moments=[0], imgPath=r"C:\Users\seany\Documents\Offner Research\RADMC_test\RADMC_inputs 05.12.22_21'54'18")

# plotMomentMaps([0])

plotMomentMaps()
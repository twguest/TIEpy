# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

import numpy as np

from felpy.utils.vis_utils import Grids


from matplotlib import pyplot as plt
from scipy.constants import c,h,e

def ekev2wav(ekev):
    """
    convert energy to wavelength
    
    :params ekev: radiation energy in keV
    
    :returns: radiation wavelength in m
    """
    return (h*c)/(e*ekev*1000)

def get_required_distance(W, sigma_det, wav):

    """
    Calculate the propagation distance required to satisfy sampling conditions.
    
    :param W: approximate feature size [m]
    :param sigma_det: propagated plane (detector) pixel size [m]
    :param wav: source wavelength [m]
    
    :returns zreq: required distance to satisfy sampling conditions [m]
    """

    if type(sigma_det) == list:
        sigma_det = max(sigma_det)
    
    zreq = (W*sigma_det)/(wav)
    return zreq


def fresnel_criterion(W, z, wav):
   """
   determine the frensel number of a wavefield propagated over some distance z
   
   :param W: approx. beam/aperture size [m]
   :param z: propagation distance [m]
   :param wav: source wavelength [m]
   
   :returns F: Fresnel number
   """
   
   F = W**2/(z*wav)
   return F

def nfs_sampling(dz, feature_size, dx, wav):
    """
    
    :param dz: distance from the speckle mask to the detector
    :param feature_size: size of the phase-object for a nfs experiment
    :param dx: detector pixel-size
    :param wav: radiation wavelength
    
    :returns S: sampling criterion (dz/zreq)
    :return F: Fresnel Criterion
    """
    
    S = dz/get_required_distance(feature_size, dx, wav) ### required prop. dist for sampling 
    F = fresnel_criterion(feature_size, dz, wav) ### fresnel criterion
    
    return S,F

def angular_sensitivity(dx, dz):
    return np.arctan(dx/dz)


def theta(x,r, delta = 4e-07):
    return (2*r*delta)/np.sqrt(r**2-(x-r)**2)

def shift(theta, dx, z):
    return np.tan(theta)*z/dx
    
 
if __name__ == '__main__':
    
    ### usage
    
    ### Define the magnifcation from the source aperture at each of the imaging
    ### positions

    
    ### Define z1: the source to sample distance, and dz: the sample to det. distance
    z1 = 2 ### src to sample distance
    dz = 3
    
    zD = 0.5
    dzD = 4.5
    
    div = np.tan(13.5e-03/8)
    
    M01 = 1 #+ np.tan(div)*z1/1.5e-03
    M12 = 1 #+ np.tan(div)*dz/(M01*1.5e-03)
    M02 = 1 #+ np.tan(div)*(dz+z1)/1.5e-03
    
    M0D = 1 #+ np.tan(div)*zD/1.5e-03
    MD2 = 1 #+ np.tan(div)*dzD/(1.5e-03*M0D)
    
    r = 75e-03
    
    
    dx = 7.5e-06 ### detector pixel size
    ekev = 25.0 ### energy in keV
    
    dw = 15e-03 ### detector width 
    fov = np.linspace(dw/M02,r, 400)  ### detector field of view of source/cyl
    x = np.linspace(0,dw,400)

    
    ### plotting stuff here
    plots = Grids(global_aspect = 2, scale = 2)
    plots.create_grid(n = 1, m = 2, sharex = False, sharey = False)
    [ax1, ax2] = plots.axes
    
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("Phase Gradient ($\mu rad$)")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("Feature Displacement (pixels.)")
    ax1.legend()
    
    for r in [.75, 1.00, 1.50]:   

        
        t = theta(fov, r)
        t[np.where(np.isnan(t))] = 0

        ax1.plot(x*1e3, t*1e6, label = "{} mm".format(r*1e2))
        ax2.plot(x*1e3, shift(t, dx, dz+z1))

    for feature_size in [6e-06, 9e-06, 15e-06, 55e-06]:
        print("")
        print("Analyser")
        print("**** Feature Size: {}".format(feature_size))  
        
        S, F = nfs_sampling(dz, feature_size, dx, ekev2wav(ekev))
        print("Sampling Criterion (req. > 1): {}".format(S))
        
        
        print("Fresnel Criterion (req. > 1): {}".format(F))
        
        
        print("Feature Size @ Detector: {} (should be > dx) (pixels)".format(M12*feature_size/dx))
        
        #print("Angular Sensitivity: {}".format(angular_sensitivity(dx, dz)))
    

    for feature_size in [30e-06, 45e-06]:
        print("")
        print("Diffuser")
        print("**** Feature Size: {} um".format(feature_size*1e6))  
        
        S, F = nfs_sampling(dzD, feature_size, dx, ekev2wav(ekev))
        print("Sampling Criterion (req. > 1): {}".format(S))
        
        
        print("Fresnel Criterion (req. > 1): {}".format(F))
        
        
        print("Feature Size @ Detector: {} (should be > dx) (pixels)".format(MD2*feature_size/dx))
        
        #print("Angular Sensitivity: {}".format(angular_sensitivity(dx, dz)))
    
    ### plotting stuff here
    plots = Grids(global_aspect = 1.5, scale = 0.75)
    plots.create_grid(n = 2, m = 1, sharex = True, sharey = False)
    [ax1, ax2] = plots.axes
    
    plots.add_global_colorbar("Intensity (a.u.)", fontsize = 11, pad = 0.1, cmap = 'viridis', orientation = 'horizontal')

    y = np.ones([400,400])
    y *= shift(t, dx, dz+z1)
    
    ax1.imshow(y)
    ax1.set_aspect('auto')
    ax1.set_title("Simulated", pad = 2)
    ax1.set_ylabel("y(mm)")
    #ax1.colorbar(orientation="horizontal")
    ax2.set_title("Imaging and Medical Beamline", pad = 1)
    

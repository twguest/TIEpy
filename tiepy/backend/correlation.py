# -*- coding: utf-8 -*-

from scipy.signal import correlate2d as corr2d
import numpy as np

from tiepy.backend.processing import frame
"""
### TO-DO: function for finding point of max correlation between two images
### TO-DO: function for loading tiffs
### TO-DO: function for normalising intensity
"""





def get_coords(arr, shape):
    
    
    sh = arr.shape
    xmin, xmax = shape[0], sh[0]-shape[0]
    ymin, ymax = shape[1], sh[1]-shape[1]

    coords_x = np.arange(xmin, xmax)
    coords_y = np.arange(ymin, ymax)
    
    return coords_x, coords_y

def max_index(arr):
    return [np.argmax(arr, axis=1)[0], np.argmax(arr, axis=0)[0]]


def test_function(x, shift = 0):
    return (x+shift)*2

if __name__ == '__main__':
    
    from copy import copy
    from matplotlib import pyplot as plt
    
    global_shift = 1
    nx = ny = 150    
    
    i1 = np.ones([nx,ny])*np.arange(nx)
    i2 = np.ones([nx,ny])*np.arange(2,nx+2)
    
    i2 = frame(i2, c = [75,75], shape = [25, 25])
    
    #i1 = frame(i1, c = [75,75], shape = [50, 60])
 
    xc, yc = get_coords(i1, shape = [25,25])
    
    corr = corr2d(i1,i2, 'valid')
    plt.imshow(corr)
    
    mx, my = max_index(corr)
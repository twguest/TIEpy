# -*- coding: utf-8 -*-

import numpy 

def frame(arr, c, shape = [5,5]):
    """
    create a frame/window centered at c in arr.
    
    :param arr: numpy arry
    :param c: center of frame
    :param shape: size of the frame
    """
    sh = arr.shape
    
    return arr[c[0]-shape[0]//2:c[0]+shape[0]//2,c[1]-shape[1]//2:c[1]+shape[1]//2]
    
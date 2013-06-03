from bhmie import bhmie
import numpy as np


def calcq(optical,a):
    """
    optical is Nx3 numpy array
    first column: wavelength (um)
    second column: n (real optical const)
    third column: k (imaginary optical const)
    
    a is the grain size in um
    """
    qabs = np.zeros(optical.shape[0])
    for i,row in enumerate(optical):
        qext,qsca = bhmie(2*np.pi*a/row[0],row[1]+row[2]*1j,2) #2 is just a dummy number
        qabs[i] = qext-qsca
    return qabs

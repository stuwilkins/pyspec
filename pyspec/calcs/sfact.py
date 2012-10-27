"""Module to obtain and calculate x-ray scattering factors

This module provides access to the x-ray scattering factors.
It uses the DABAX files from the ESRF which are in type
of SPEC data format. Currently there is access to F0 and F1,F2.

The module uses datafiles which are assumed to be in the path
specified relative to the directory of the module. For example
by default the datafiles are expeced in ./data/<datafile>

When the module is loaded, the datafiles are read into global 
variables to cache their results for faster access. This is at
the detriment of speed of importing the modules.

Elements are specified by their name. If they are strings, this 
can also contain the valence (charge) in the case of F0

"""


from __future__ import with_statement
import scipy
import numpy as np
from scipy import interpolate, pi, inner, outer
import pickle
import os.path

# Do the caching part

F0 = {}
F1F2 = {}
AtomicConstants = {}

def getF0Params(element = None):
    """Returns parameters for caluclating F0

    f0 is approximated by by a function:

    f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
                i=1,4
 
    where k = sin(theta) / lambda and c, a_i and b_i
    are the so called Cromer-Mann coefficients.
    
    The function returns the coefficients as a tuple of
    the form (a, c, b)

    Parameters
    ----------

    element : string containing identifier of element

    Examples
    -------
    
    >>> sfact.getF0Params('He')


    """
    global F0
    thisElement = None
    with open(os.path.join(datadir, 'f0_InterTables.dat'),'r') as infile:
        for line in infile:
            if line.split()[0] == '#S':
                thisElement = line.split()[2]
            elif line[0] == '#':
                continue
            else:
                data = line.split()
                for i in range(len(data)):
                    data[i] = float(data[i])
                a = data[0:4]
                b = data[4]
                c = data[5:9]
                F0[thisElement] = (np.array(a), b, np.array(c))
                if thisElement == element:
                    infile.close()
                    return F0[thisElement]

    return None

def getF0(element, q = 0.):
    """Returns the scattering factor f_0

    We define the magnitude of Q as per physicists, with 
    |k| = 2\pi / \lambda

    Parameters
    ----------
    element : string 
        Element name  
    k : float or arraylike
        The magnitude of q to calculate the form factor at.
        This can be a single value or an array of values

    Examples
    --------

    >>> getF0('Fe', array([0.1, 0.2, 0.3]))
        Calculate f_0 for Fe at |Q| = 0.1, 0.2 and 0.3
    

    """
    
    global F0
    return calcF0(F0[element], q)

def getRealF0(element, x = 0):
    """Returns the electron density distribution

    Parameters
    ----------
    element : string
        Element name
    x : float or array like
        The position away from the atoms center in angstroms

    Examples
    --------

    >>> getRealF0('Fe', array([-0.1, 0, 0.1]))
        Calculate the electron density at x = -0.1, 0 and 0.1

    """
    global F0
    return calcRealF0(F0[element], x)
    
    
def calcRealF0(params, x = 0):
    """Calculate the electron density from Cromer-Mann coefficients

    The value returned is the electron density distribution of the 
    element, obtained by inverse fourier transforming the atomic
    form factor as supplied by the Cromer-Mann method.

    Parameters
    ----------
    params : tuple
        The Cromer-Mann coefficients
    x : float or arraylike
        Position in angstroms from the center of the atom

    """
                 
    r = scipy.exp(-1.0 * pow(pi, 2) * outer(1. / params[2], pow(x, 2)))
    r = ((params[0] / scipy.sqrt(params[2] / pi)) * r.T).T
    r = r.sum(0)
    return r

def calcF0(params, k = 0):
    """Calculate f_0 from Cromer-Mann coefficients

    We define the magnitude of Q as per physicists, with 
    |k| = 2\pi / \lambda

    """
    k = np.asarray(k)
    if k.ndim != 0:
        k = scipy.sqrt(pow(k, 2).sum(-1))
    stol = k / (4.0 * pi)
    f0 = (params[0] * scipy.exp(-1.0 * outer(params[2], pow(stol, 2))).T).T
    f0 = params[1] + f0.sum(0)
    return f0

def getF1F2(element):
    """Returns F1F2 Params"""
    return F1F2[element]

def getAtomicConstants(element = None):
    """Returns Atomic Constant for given element

    Array is in format:
    0 = Atomic Number
    1 = Atomic Radius [A]
    2 = CovalentRadius [A]
    3 = AtomicMass
    4 = BoilingPoint [K]
    5 = MeltingPoint [K]
    6 = Density [g/ccm]
    7 = Atomic Volume
    8 = CoherentScatteringLength [1E-12cm]
    9 = IncoherentX-section [barn]
    10 = Absorption@1.8A [barn]
    11 = DebyeTemperature [K]
    12 = ThermalConductivity [W/cmK]

    """
    return AtomicConstants[element]

def getAtomicConstantsParams(element = None):
    """Returns Atomic Constants

    Array is in format:
    0 = Atomic Number
    1 = Atomic Radius [A]
    2 = CovalentRadius [A]
    3 = AtomicMass
    4 = BoilingPoint [K]
    5 = MeltingPoint [K]
    6 = Density [g/ccm]
    7 = Atomic Volume
    8 = CoherentScatteringLength [1E-12cm]
    9 = IncoherentX-section [barn]
    10 = Absorption@1.8A [barn]
    11 = DebyeTemperature [K]
    12 = ThermalConductivity [W/cmK]

    """
    global AtomicConstants
    with open(os.path.join(datadir, 'AtomicConstants.dat'),'r') as infile:
        for line in infile:
            if line.split()[0] == '#S':
                s = line.split()
                thisElement = s[2]
                thisZ = np.array([s[1]]).astype('float32')
            elif line[0] == '#':
                continue
            else:
                data = np.array(line.split()).astype('float32')
                AtomicConstants[thisElement] = np.concatenate((thisZ, data))
                
def getF1F2Params(element = None):
    """Returns f1 and f2 scattering factors"""

    alldata = np.array([])
    global F1F2

    with open(os.path.join(datadir, 'f1f2_Henke.dat'),'r') as infile:
        for line in infile:
            if line.split()[0] == '#S':
                if len(alldata):
                    f1 = interpolate.interp1d(alldata[:,0], alldata[:,1] - thisZ)
                    f2 = interpolate.interp1d(alldata[:,0], alldata[:,2] * -1.)
                    F1F2[thisElement] = (f1, f2)
                    if thisElement == element:
                        infile.close()
                        return F1F2[element]

                s = line.split()
                thisElement = s[2]
                thisZ = int(s[1])

                alldata = np.array([])

            elif line[0] == '#':
                continue
            else:
                data = np.array(line.split()).astype('float32')
                if not len(alldata):
                    alldata = data
                else:
                    alldata = np.vstack((alldata, data))
    return alldata

# Now read into global variables

def _test():
    # Routune here to test the module and see if it works
    pass

if __name__ == "__main__":
    _test()
else:
    datadir = os.path.join(os.path.dirname(__file__), 'data')
    getF0Params()
    getF1F2Params()
    getAtomicConstantsParams()


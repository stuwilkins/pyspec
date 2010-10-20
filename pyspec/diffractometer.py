#
# diffractometer.py (c) Stuart B. Wilkins 2010, Sven Partzsch 2010
#
# $Id: fit.py 48 2010-10-16 18:05:14Z stuwilkins $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/fit.py $
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Part of the "pyspec" package
#
import numpy and np
import exceptions

"""
tardis lab frame:
Z up, Y along X-ray beam, X = Y x Z
sample rotations
mu    : along +Z     -> S'    (mu-frame)
theta : along +X'    -> S''   (theta-frame)
chi   : along +Y''   -> S'''  (chi-frame)
phi   : along +X'''  -> S'''' (phi-frame)
detector rotations
mu    : along +Z     -> S'    (mu-frame)
delta : along +X'    -> S*    (delta-frame)
gamma : along +Z*    -> S**   (gamma-frame)
"""

#
# Rotation Matricies 
#

def rotX(alpha):
    """Rotation matrix for a vector along X by angle alpha """
    rotMatrix = np.array([[ 1,          0,           0],
                          [ 0, cos(alpha), -sin(alpha)],
                          [ 0, sin(alpha),  cos(alpha)]])
    return rotMatrix

def rotY(alpha):
    """Rotation matrix for a vector along X by angle alpha """

    rotMatrix = np.array([[ cos(alpha), 0, sin(alpha)],
                          [          0, 1,          0],
                          [-sin(alpha), 0, cos(alpha)]])
    return rotMatrix

def rotZ(alpha):
    rotMatrix = np.array([[ cos(alpha), -sin(alpha), 0],
                          [ sin(alpha),  cos(alpha), 0],
                          [          0,           0, 1]])
    return rotMatrix

class Diffractometer():
    """Diffractometer class

    This class provides various functions to perform calculations to
    and from the sample and instrument frames"""

    def __init__(self, mode = 'sixc'):
        """Initialize the class
        'mode' defines the type of diffractometer"""
        self.mode = mode

    def setAllAngles(self, angles, mode = 'deg'):
        """Sets angles for calculation.
        Angles are expected in spec 'sixc' order:
        Delta     Theta       Chi       Phi        Mu     Gamma"""
        self.settingAngles = angles
        if mode == 'deg':
            self.settingAngles = self.settingAngles / 180.0 * pi

    def setAngles(self, delta = None, theta = None, chi = None,
                  phi = None, mu = None, gamma = None, mode = 'deg'):
        """Set the angles for calculation"""
        
        # Internally angles are stored in sixc (spec) order
        # Delta     Theta       Chi       Phi        Mu     Gamma

        maxlen = array([theta.len, delta.len, gamma.len, mu.len, phi.len, chi.len]).max()
        self.settingAngles = zeros((maxlen, 6))
        for aa, i in zip([delta, theta, chi, phi, mu, gamma], range(6)):
            if aa.len == maxlen:
                self.settingAngles[:,i] = aa
            elif aa.len == 1:
                self.settingAngles[:,i] = ones(maxlen) * aa
            elif aa is not None:
                raise exceptions.ValueError("%s must be numpy array or ")
            
        if mode == 'deg':
            self.settingAngles = self.settingAngles / 180.0 * pi
        

def thdelgam2Qxyz(theta, delta, gamma, waveLen = 2*pi,
                  mu = 0.0, mode = 'deg', shortForm = True, verbose = False):

    """
    transformation of set (theta, delta, gamma) into (Qx, Qy, Qz)
    in tardis sample frame (theta-frame)
    tardis sample frame for theta = 0:
    Qx, Qy, Qz along X, Y, Z, respectively

    angles as arrays or float
    all arrays have to have the same shape!
    """
        
    # for angels in degree
    if mode == 'deg':
        mu    = mu    / 180.0 * np.pi
        theta = theta / 180.0 * np.pi
        delta = delta / 180.0 * np.pi
        gamma = gamma / 180.0 * np.pi

    # wave vector length in 1/A, energy in eV
    kl = 2 * np.pi / waveLen 

    # wave vector for all diffractometer angles zero
    k0 = array([ 0, kl, 0]).T

    # initial wave vector in theta-frame
    ki = dot( rotX(-theta), dot( rotZ(-mu), k0 ) )

    # final   wave vector in theta-frame
    kf = dot( rotX(-theta+delta), dot( rotZ(gamma), k0 ) )

    #   scattering vector in theta-frame
    q  = kf - ki

    return q


def Qxyz2Qsam(Qxyz, chi=0.0, phi=0.0, mode = 'deg', shortForm = False, choise = 'tardis', verbose = False):

    """
    transformation of (Qx, Qy, Qz) from the theta frame into the sample (phi) frame
    in          'tardis' or 'sixc' choise
    up        :    x           z
    along beam:    y           y
    right hand:    z           x
    """

    # transformation from sixc to tardis
    Vts = array([ [ 0, 0, -1],
                  [ 0, 1,  0],
                  [ 1, 0,  0] ])

    if mode == 'deg':
        chi = chi / 180.0 * pi
        phi = phi / 180.0 * pi

    # identity
    if   type(chi) == numpy.ndarray:
        ident = chi * 0.0 + 1.0
    elif type(phi) == numpy.ndarray:
        ident = phi * 0.0 + 1.0
    else:
        ident = 1.0    

    if shortForm == False:
        Qphi = dot( rotMat(-phi, 'x'), dot( rotMat(-chi, 'y'), Qxyz ) )  
        if choise == 'sixc':
            Qphi = dot( Vts.T, Qphi )
    else:
        r11 = ident*          cos(chi)
        r12 = ident* 0.0
        r13 = ident*         -sin(chi)
        r21 = ident* sin(phi)*sin(chi)
        r22 = ident* cos(phi)
        r23 = ident* sin(phi)*cos(chi)
        r31 = ident* cos(phi)*sin(chi)
        r32 = ident*-sin(phi)
        r33 = ident* cos(phi)*cos(chi)

        Qphi = concatenate(( array([ r11*Qxyz[0] + r12*Qxyz[1] + r13*Qxyz[2] ]) ,
                             array([ r21*Qxyz[0] + r22*Qxyz[1] + r23*Qxyz[2] ]) ,
                             array([ r31*Qxyz[0] + r32*Qxyz[1] + r33*Qxyz[2] ]) ))

        if choise == 'sixc':
            Qphi = concatenate(( array([ Vts[0,0]*Qphi[0] + Vts[1,0]*Qphi[1] + Vts[2,0]*Qphi[2] ]) ,
                                 array([ Vts[0,1]*Qphi[0] + Vts[1,1]*Qphi[1] + Vts[2,1]*Qphi[2] ]) ,
                                 array([ Vts[0,2]*Qphi[0] + Vts[1,2]*Qphi[1] + Vts[2,2]*Qphi[2] ]) ))

    if verbose == True:
        if mode == 'deg':
            chi = chi * 180.0 / pi
            phi = phi * 180.0 / pi
        print '(chi, phi) = (%s, %s)' % (chi, phi)
        print '(Qx, Qy, Qz) = \n%s' % (Qxyz)
        print 'has been transformed into (Qx, Qy, Qz)_phi_%s = \n%s' % (choise, Qphi)

    return Qphi
    

def angle2Qsam(mu, theta, chi, phi, delta, gamma, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False):

    Qxyz = thdelgam2Qxyz(theta, delta, gamma, energy = energy, mu = mu, mode = mode, shortForm = shortForm, verbose = verbose)
    Qsam = Qxyz2Qsam(Qxyz, chi = chi, phi = phi, mode = 'deg', shortForm = False, choise = choise, verbose = False)

    return Qsam


####################################
#
# visualization for test
#
####################################

def k2Dshow(k, name = 'Q', conRoi = [1, 325, 1, 335]):

    kabs = sqrt( k[0]**2 + k[1]**2 + k[2]**2 )
    figure()
    ax = subplot(1, 4, 1)
    imshow(k[0])
    ax.set_ylabel(name+'x')
    ax = subplot(1, 4, 2)
    imshow(k[1])
    ax.set_ylabel(name+'y')
    ax = subplot(1, 4, 3)
    imshow(k[2])
    ax.set_ylabel(name+'z')
    ax = subplot(1, 4, 4)
    imshow( kabs )
    ax.set_ylabel('|'+name+'|')
   
    figure()
    cutAlpha = [0, 45, 90, 135]
    conCen   = [ 162.5, 167.5]
    cutXY    = ['x', 'x', 'y', 'x']
    for j in range(4):
            
        # xy-values, measured and fited values of j-th line cut
        linX, linY = getLineXY(conRoi, conCen, cutAlpha[j])
        linKx      = getSubSet(k[0],  conRoi, linX, linY, verbose = False)
        linKy      = getSubSet(k[1],  conRoi, linX, linY, verbose = False)
        linKz      = getSubSet(k[2],  conRoi, linX, linY, verbose = False)
        linKabs    = getSubSet(kabs,  conRoi, linX, linY, verbose = False)
       
        # new subplot
        ax = subplot( 2, 2, j+1)
        
        # plot value, and result of j-th line cut
        if cutXY[j] == 'x':
            ax.plot(linX, linKx  , '-bo', label=name+'x'  )
            ax.plot(linX, linKy  , '-r^', label=name+'y'  )
            ax.plot(linX, linKz  , '-gv', label=name+'z'  )
            ax.plot(linX, linKabs, '-c*', label='|'+name+'|' )
        else:
            ax.plot(linY, linKx  , '-bo', label=name+'x'  )
            ax.plot(linY, linKy  , '-r^', label=name+'y'  )
            ax.plot(linY, linKz  , '-gv', label=name+'z'  )
            ax.plot(linY, linKabs, '-c*', label='|'+name+'|' )
        
        # title and labels of j-th line cut
        yLabel = 'alpha = %d' % (cutAlpha[j])
        ax.set_ylabel(yLabel)
        ax.set_xlabel(cutXY[j])

        legend(loc='best')



####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    #rotM = rotMat(10, kind = 'x', mode = 'deg', verbose = True)

    q = thdelgam2Qxyz( 10, 20, 1, waveLen = 2*pi, mu = 5, mode = 'deg', shortForm = True, verbose = True)
    #q = thdelgam2Qxyz( array([5,10]), array([10,20]), array([0.5,1]), waveLen = 2*pi, mu = array([2.5,5]),
    #                   mode = 'deg', shortForm = True, verbose = True)



    ###########
    # comparison with sixc simulation mode, checked all six angles
    ###########

    """
    4.SIXC> setlat 6.28319 6.28319 6.28319 90 90 90
    (UB recalculated from orientation reflections and lattice.)
    """
    #14.SIXC> p UB
    #UB["0"] = 0.999999253114992
    #UB["1"] = -1.53360282755493e-16
    #UB["2"] = -1.60814167124132e-16
    #UB["3"] = -7.4538843686392e-18
    #UB["4"] = 0.999999253114992
    #UB["5"] = -6.12302719591125e-17
    #UB["6"] = 0
    #UB["7"] = 0
    #UB["8"] = 0.999999253114992
    """
    24.SIXC> LAMBDA = hc_over_e / 640

    34.SIXC> pa

    Six-Circle Geometry, Omega fixed (four circle, Mu = Gamma = 0) (mode 0)
    Sector 0

    Primary Reflection (at lambda 1.54):
     del th chi phi mu gam = 60 30 0 0 0 0 
                     H K L = 1 0 0

    Secondary Reflection (at lambda 1.54):
     del th chi phi mu gam = 60 30 0 90 0 0 
                     H K L = 0 1 0

    Lattice Constants (lengths / angles):
              real space = 6.283 6.283 6.283 / 90 90 90
        reciprocal space = 1 1 1 / 90 90 90

    Azimuthal Reference:
                   H K L = 0 0 1
               sigma tau = 0 0

    Monochromator:
                  Lambda = 19.3725 

    Cut Points:
      del   th  chi  phi   mu  gam
     -180 -180 -180 -180 -180 -180

    """

    #Qsam = angle2Qsam(mu, theta, chi, phi, delta, gamma, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    """
    Qsam = angle2Qsam(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    python: In [81]: Qsam
            Out[81]: array([ 0.,  0.,  0.])
    38.SIXC> wh

    H K L =  0  0  0
    Alpha = 0  Beta = 0  Azimuth = 180
    Two Theta = 0  Omega = 0  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000

    
    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 30.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    python: In [83]: Qsam
            Out[83]: array([ 0.16788546,  0.        ,  0.        ])
    40.SIXC> wh

    H K L =  0.16789  3.9558e-09  0
    Alpha = 0  Beta = 0  Azimuth = -90
    Two Theta = 30  Omega = 1.35e-06  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    30.0000   15.0000    0.0000    0.0000    0.0000    0.0000         


    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 40.0, 0.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [85]: Qsam
    Out[85]: array([ 0.22101043, -0.01933591,  0.        ])
    42.SIXC> wh

    H K L =  0.22101  -0.019336  0
    Alpha = 0  Beta = 0  Azimuth = -90
    Two Theta = 40  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000    0.0000    0.0000

    
    Qsam = angle2Qsam(0.0, 15.0, 0.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [87]: Qsam
    Out[87]: array([ 0.22048884, -0.02045445,  0.0282672 ])
    44.SIXC> wh

    H K L =  0.22049  -0.020455  0.028268
    Alpha = 0  Beta = 5  Azimuth = -87.318
    Two Theta = 40.259  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000    0.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 0.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [89]: Qsam
    Out[89]: array([ 0.21921356, -0.01569504,  0.08458648])
    46.SIXC> wh

    H K L =  0.21922  -0.015695  0.084588
    Alpha = 10  Beta = 5  Azimuth = -92.851
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000    0.0000    0.0000   10.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 30.0, 0.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [91]: Qsam
    Out[91]: array([ 0.14755127, -0.01569504,  0.18286083])
    48.SIXC> wh

    H K L =  0.14755  -0.015695  0.18286
    Alpha = 16.131  Beta = 16.618  Azimuth = -89.602
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000   30.0000    0.0000   10.0000    5.0000

    
    Qsam = angle2Qsam(10.0, 15.0, 30.0, 25.0, 40.0, 5.0, energy = 640.0, mode = 'deg', shortForm = True, choise = 'sixc', verbose = False)
    In [93]: Qsam
    Out[93]: array([ 0.14035988,  0.04813332,  0.18286083])
    50.SIXC> wh

    H K L =  0.14036  0.048134  0.18286
    Alpha = 16.131  Beta = 16.618  Azimuth = -89.602
    Two Theta = 42.574  Omega = -4.9999  Lambda = 19.3725

    Delta     Theta       Chi       Phi        Mu     Gamma 
    39.9999   15.0000   30.0000   25.0000   10.0000    5.0000
    """

#
# transformations.py (c) Stuart B. Wilkins 2010 and (c) Sven Partzsch 2010
#
# $Id$
# $HeadURL$
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

import os
import numpy as np
import matplotlib.pyplot as plt
from   pyspec import spec, fit, fitfuncs
from   pyspec.diffractometer import Diffractometer
from   pyspec.ccd.PrincetonSPE import *
from   pyspec.ccd.plotter import PlotGrid
import gridder
from ConfigParser import ConfigParser

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

class MyConfigParser(ConfigParser):
    def readAllLocations(self, filename):
        _f = os.path.split(filename)

        if _f[0] == '':
        
            locations = []
            if os.name is 'posix':
                if os.environ.has_key('HOME'):
                    locations.append(os.environ['HOME'] + os.path.sep + ".pyspec")
                locations.append('/usr/local/pyspec/etc')
                locations.append('/etc/pyspec')
            
            for l in locations:
                if os.path.isfile(l):
                    self.read(l)
                    break
        else:
            self.read(l)
            
    def getWithDefault(self,section, option, default):
        return Config.get(section, option, vars = { option : default })
    def _getWithConvert(self,_conv_fcn, section, option, default):
        try:
            val = self.getWithDefault(section, option, default)
        except:
            raise Exception("Unable to read option %s from config file." % (option))
        try:
            val = _conv_fcn(val)
        except:
            raise Exception("Unable to convert option %s to correct datatype." % (option))
        return val
    def getFloat(self, *args, **kwargs):
        self._getWithConvert(float, *args, **kwargs)
    def getInt(self, *args, **kwargs):
        self._getWithConvert(int, *args, **kwargs)
    def getFloat(self, *args, **kwargs):
        self._getWithConvert(float, *args, **kwargs)
            
#
# global help functions
#

def getAreaSet(image, roi):

    """
    selects a region of interest (ROI) from a CCD image
    
    image : array of the image values in the full window ( [   1, 325,   1, 335 ] )
    roi   : region of interest [xmin, dx, ymin, dy], e.g. [137, 51, 142, 51] for center

    output
    cut   : [ [val00, val01, ..., val0m],
              [val10, val11, ..., val1m],
              ...,
              [valn0, valn1, ..., valnm] ]
    valij is value at (x[i], y[j]) 
    """

    # data points in region of interest at (y, x)
    cut = image[roi[2]-1:roi[2]+roi[3]-1, roi[0]-1:roi[0]+roi[1]-1 ]

    return cut


def getAreaXY(roi):

    """
    calculates (x, y) = (x[0], x[1]) set for each point in region of interest

    output
    z  : [ [ [x00, x10, ..., xn0], [x01, x11, ..., xn1], ..., [x0m, x1m, ..., xnm] ],
           [ [y00, y10, ..., yn0], [y01, y01, ..., yn1], ..., [y0m, y1m, ..., ynm] ] ]
    """

    z  = np.array(np.meshgrid(np.arange(roi[0], roi[0] + roi[1]), np.arange(roi[2], roi[2] + roi[3])))
            
    return z 

#
# image processor class
#

class ImageProcessor():
    """Image Processor class

    This class provides the processing of single and sets of CCD-images
    in more detail: each pixel is transformed into reciprocal space
    the set of reciprocal vectors and intensities is gridded on a regular cuboid
    the needed informations can be provided by a spec scan from the pyspec package"""

    def __init__(self, configfile = None, ccdname = 'CCD'):
        # set parameters to configure the CCD setup
        # detector distance 30cm and detector pixel size 20um

        self.ccdName = ccdname
        #self.readConfigFile(configfile)
        
        self.detDis      = 300    # in mm
        self.detPixSizeX = 0.020  # in mm
        self.detPixSizeY = 0.020  # in mm
        # detector size in pixel
        self.detSizeX    = 1300
        self.detSizeY    = 1340
        self.detX0       = 1300 / 2.0
        self.detY0       = 1340 / 2.0
        self.detAngle    = 0.0
        # image treatment
        self.setName     = 'Set #'
        self.setNum      = 1
        self.conRoi      = None
        self.setFrameMode(1)
        # gridder options
        self.setGridOptions()
        # plot options
        self.plotFlag2D  = 7
        self.plotFlag1D  = 7
        self.logFlag1D   = 0
        self.logFlag2D   = 0
        self.fit1D       = False
        self.fitType     = 'lor2a'
        self.histBin     = 50
        # output information
        self.opTitle     = 'Image Processing'
        self.projName    = ''
        self._makeSetInfo()
        self._makeModeInfo()
        self.opProcInfo  = ''

        # Figure Size
        self._defaultFigureSize = (11, 8.5)

    def readConfigFile(self, filename):
        config = MyConfigParser()
        config.read(filename)
        config.get(self.ccdname, '', vars = {})
        
    #
    # set and get part
    #

    def setProjName(self, projName):
        """Set the name of the current project

        projName: name of the current project"""

        self.projName = projName

    def getProjName(self, projName):
        """Get the name of the current project

        projName: name of the current project"""

        return self.projName

    def setDetectorProp(self, detPixSizeX, detPixSizeY, detSizeX, detSizeY, detX0, detY0):
        """Set properties of used detector

        detPixSizeX : detector pixel size in detector X-direction (float in mm)
        detPixSizeY : detector pixel size in detector Y-direction (float in mm)
        detSizeX    : detector no. of pixels (size) in detector X-direction (integer)
        detSizeY    : detector no. of pixels (size) in detector Y-direction (integer)
        detX0       : detector X-coordinate of center for reference gamma-value (float)
        detY0       : detector Y-coordinate of center for reference delta-value (float)"""
        
        self.detPixSizeX = detPixSizeX  
        self.detPixSizeY = detPixSizeY
        self.detSizeX    = detSizeX
        self.detSizeY    = detSizeY
        self.detX0       = detX0
        self.detY0       = detY0

    def getDetectorProp(self):
        """Get properties of used detector

        detPixSizeX : detector pixel size in detector X-direction (float in mm)
        detPixSizeY : detector pixel size in detector Y-direction (float in mm)
        detSizeX    : detector no. of pixels (size) in detector X-direction (integer)
        detSizeY    : detector no. of pixels (size) in detector Y-direction (integer)
        detX0       : detector X-coordinate of center for reference gamma-value (float)
        detY0       : detector Y-coordinate of center for reference delta-value (float)"""

        return self.detPixSizeX, self.detPixSizeY, self.detSizeX, self.detSizeY, self.detX0, self.detY0

    def setDetectorDis(self, detDis):
        """Set the detector distance

        detDis : detector distance (float in mm)"""
        
        self.detDis = detDis

    def getDetectorDis(self):
        """Get the detector distance

        detDis : detector distance (float in mm)"""

        return self.detDis

    def setDetectorAngle(self, detAng):
        """Set the detector miss alignement angle in deg

        detAng : detector miss alignement angle in deg"""

        self.detAngle = detAng

    def getDetectorAngle(self, detAng):
        """Get the detector miss alignement angle in deg

        detAng : detector miss alignement angle in deg"""

        return self.detAngle

    def setBins(self, binX, binY):
        """Set no. of bins. Takes them into acount for pixel size, detector size and detector center
        
        binX : no. of pixels along detector X-direction which are bined
        binY : no. of pixels along detector Y-direction which are bined"""
        
        self.binX = binX
        self.binY = binY
        self._applyBins()

    def getBins(self, binX, binY):
        """Set no. of bins. Takes them into acount for pixel size, detector size and detector center
        
        binX : no. of pixels along detector X-direction which are bined
        binY : no. of pixels along detector Y-direction which are bined"""

        return self.binX, self.binY

    def setSetSettings(self, waveLen, imFilesNames, settingAngles, intentNorm, UBmat, setName, setNum, setSize):
        """Set the settings for the set 

        The set settings are:
        waveLen       : wavelength of the X-rays (float in Angstrom)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81"""

        self.waveLen       = conScan.wavelength
        self.imFileNames   = conScan.getCCDFilenames()
        self.darkFileNames = conScan.getCCDFilenames(dark = 1)
        self.settingAngles = conScan.getSIXCAngles()
        self.intentNorm    = conScan.Ring
        self.UBmat         = conScan.UB
        self.setName       = 'Scan #'
        self.setNum        = conScan.scan
        self.setSize       = self.settingAngles.shape[0]

    def getSetSettings(self, waveLen, imFilesNames, settingAngles, intentNorm, UBmat, setName, setNum, setSize):
        """Set the settings for the set 

        The set settings are:
        waveLen       : wavelength of the X-rays (float in Angstrom)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81"""

        return self.waveLen, self.imFilesNames, self.settingAngles, self.intentNorm, self.UBmat, self.setName, self.setNum, self.setSize

    def setSpecScan(self, conScan):
        """Set the settings for the set from the considered pyspec scan object

        conScan : pyspec scan object which contains the needed information

        The set settings are:
        temperature   : temperature of the sample (float in Kelvin)
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81"""

        self.conScan = conScan
        self._setBySpec()
   
    def getSpecScan(self):
        """Get the pyspec scan object which was used for the set settings
        
        conScan : pyspec scan object which contains the needed information"""

        return self.conScan

    def setConRoi(self, conRoi):
        """Set the considered region of interest 

        conRoi : [xMin, xStep, yMin, yStep]"""
        
        self.conRoi = conRoi

    def getConRoi(self):
        """Get the considered region of interest 

        conRoi : [xMin, xStep, yMin, yStep]"""
        
        return self.conRoi

    def setFrameMode(self, mode):
        """Set the mode of the output frame for (Qx, Qy, Qz)

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame)"""

        self.frameMode = mode
        self._preAxesLabels()
        self._makeModeInfo()

    def getFrameMode(self):
        """Get the mode of the output frame for (Qx, Qy, Qz)

        mode: : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame)"""

        return self.frameMode

    def setProcIm(self, procImSelect = None):
        """Set the selection of the images for processing

        procImSelect : list with the images which will be processed, all if None"""
                
        if procImSelect == None:
            self.procImSelect = range(self.setSize)
        else:
            self.procImSelect = procImSelect
        
    def getProcIm(self, procImSelect = None):
        """Get the selection of the images for processing

        procImSelect : list with the images which will be processed, all if None"""
                
        return self.procImSelect

    def setGridOptions(self, Qmin = None, Qmax = None, dQN = None):
        """Set the options for the gridding of the dataset

        Qmin : minimum values of the cuboid [Qx, Qy, Qz]_min
        Qmax : maximum values of the cuboid [Qx, Qy, Qz]_max
        dQN  : no. of grid parts (bins)     [Nqx, Nqy, Nqz]"""

        self.Qmin = Qmin
        self.Qmax = Qmax
        self.dQN  = dQN

    def getGridOptions(self):
        """Get the options for the gridding of the dataset

        Qmin : minimum values of the cuboid [Qx, Qy, Qz]_min
        Qmax : maximum values of the cuboid [Qx, Qy, Qz]_max
        dQN  : no. of grid parts (bins)     [Nqx, Nqy, Nqz]"""

        return self.Qmin, self.Qmax, self.dQN

    def getGridVectors(self):
        """Get the values for the underlying grid vectors

        qVal : list of the Qx, Qy, and Qz values"""

        return self.qVal

    #
    # set and get functions for plot settings
    #

    def setAxesLabels(self, axesLabels):
        """Set the plotting labels for the axes

        axesLabels : labels for the axes [xLabel, yLabel, zLabel]"""
        
        self.axesLabels = axesLabels

    def getAxesLabels(self):
        """Get the plotting labels for the axes

        axesLabels : labels for the axes [xLabel, yLabel, zLabel]"""
        
        return self.axesLabels

    def setPlotFlags(self, flag1D = 7, flag2D = 7):
        """Set the ploting flags for 1D and 2D plots

        flag1D : flag to select 2D plots
        flag2D : flag to select 1D plots

        binary code, flag & 1: intensity, flag & 2: occupation numbers of the grid parts (bins),
        flag & 4: histogram of occupation of the grid parts (bins)"""
        
        self.plotFlag1D = flag1D
        self.plotFlag2D = flag2D

    def getPlotFlags(self):
        """Get the ploting flags for 1D and 2D plots

        flag1D : flag to select 1D plots
        flag2D : flag to select 2D plots

        binary code, flag & 1: intensity, flag & 2: occupation numbers of the grid parts (bins),
        flag & 4: histogram of occupation of the grid parts (bins)"""
        
        return self.plotFlag2D, self.plotFlag1D

    def setLogFlags(self, flag1D = 0, flag2D = 0):
        """Set whether data are plotted on linear (0) or logarithmic (1) scale

        flag1D : flag to select scale of 1D plots
        flag2D : flag to select scale of 2D plots

        binary code, flag & 1: intensity, flag & 2: missing grid parts"""
        
        self.logFlag1D = flag1D
        self.logFlag2D = flag2D
        
    def getLogFlags(self):
        """Get whether data are plotted on linear (0) or logarithmic (1) scale

        flag1D : flag to select scale of 1D plots
        flag2D : flag to select scale of 2D plots

        binary code, flag & 1: intensity, flag & 2: missing grid parts"""
        
        return self.logFlag1D, self.logFlag2D

    def setFit1D(self, fit1D = 0, fitType = 'lor2a'):
        """Set whether 1D lines get fitted (1)

        fit1D   : 1 (on), 0 (off) fitting of the 1D data
        fitType : type of the peak function from pyspec fitfuncs, e.g. 'lor2a'"""
        
        self.fit1D     = fit1D
        self.fit1DType = fitType

    def getFit1D(self):
        """Get whether 1D lines get fitted (1)

        fit1D   : 1 (on), 0 (off) fitting of the 1D data
        fitType : type of the peak function from pyspec fitfuncs, e.g. 'lor2a'"""
        
        return self.fit1D, self.fit1DType

    def setHistBin(self, histBin = 50):
        """Set the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        self.histBin = histBin

    def getHistBin(self, histBin = 50):
        """Get the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        return self.histBin

    def setPlotIm(self, plotImSelect = None, plotImHor = 4, plotImVer = 3):
        """Set the options for ploting the raw images

        plotImSelect : list with the raw images which will be plotted, all if None
        plotImHor    : no. of horizontal images per window, e.g. 4
        plotImVer    : no. of vertical   images per window, e.g. 3"""
        
        if plotImSelect == None:
            self.plotImSelect = range(self.setSize)
        else:
            self.plotImSelect = plotImSelect
        self.plotImHor    = plotImHor
        self.plotImVer    = plotImVer
        
    def getPlotIm(self):
        """Get the options for ploting the raw images

        plotImSelect : list with the raw images which will be plotted, all if None
        plotImHor    : no. of horizontal images per window, e.g. 4
        plotImVer    : no. of vertical   images per window, e.g. 3"""
        
        return self.plotImSelect, self.plotImHor, self.plotImVer

    #
    # get set functions for input output
    #

    def setInfoFile(self, infoFile):
        """Set the path and file name for the info file about the current processing

        infoFile : path and file name for the info file about the current processing"""

        self.infoFile = infoFile

    def getInfoFile(self):
        """Get the path and file name for the info file about the current processing

        infoFile : path and file name for the info file about the current processing"""

        return self.infoFile

    #
    # help function part
    #

    def _applyBins(self):
        """Takes binning into acount for pixel size, detector size and detector center"""
        self.detPixSizeX *= self.binX
        self.detPixSizeY *= self.binY
        self.detSizeX    /= self.binX
        self.detSizeY    /= self.binY
        self.detX0       /= self.binX
        self.detY0       /= self.binY

    def _setBySpec(self):
        """Set the settings for the set from the considered pyspec scan object

        The set settings are:
        temperature   : temperature of the sample (float in Kelvin)
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81"""

        self.temperature   = self.conScan.Tsam.mean() # in Kelvin
        self.waveLen       = self.conScan.wavelength  # in Angstrom
        self.energy        = Diffractometer.hc_over_e / self.conScan.wavelength # in eV
        self.imFileNames   = self.conScan.getCCDFilenames()
        self.darkFileNames = self.conScan.getCCDFilenames(dark = 1)
        self.settingAngles = self.conScan.getSIXCAngles()
        self.intentNorm    = self.conScan.Ring
        self.UBmat         = self.conScan.UB
        self.setName       = 'Scan #'
        self.setNum        = self.conScan.scan
        self.setSize       = self.settingAngles.shape[0]
        self.procImSelect  = range(self.setSize)
        self._makeSetInfo()

    def _preAxesLabels(self, mode = None):
        """Prepare the labels of the axes regarding the frame mode

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        if mode == None:
            mode = self.frameMode

        if mode != 4:
            self.qLabel = ['Qx', 'Qy', 'Qz']
            self.setEntLabel = '(Qx, Qy, Qz, I)'
            self.setAxesLabels([ur"Qx (\u00c5$^{-1}$)", ur"Qy (\u00c5$^{-1}$)", ur"Qz (\u00c5$^{-1}$)"])
        else:
            self.qLabel = ['H', 'K', 'L']
            self.setEntLabel = '(H, K, L, I)'
            self.setAxesLabels(['H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'])

    def _readImage(self, imNum):
        """Read in the considered region of interest of the image
        dark image subtraction and normalization by ring current"""

        if imNum % 10 == 0:
            print '%s%d image #%03d: read image' % (self.setName, self.setNum, imNum)

        # get dark image (first CCD image)
        fileName = self.darkFileNames[0]
        darkVal  = PrincetonSPEFile(fileName)[0].astype(numpy.float64)
        
        # get file name and read image of data point
        pointMon = self.intentNorm[imNum]
        fileName = self.imFileNames[imNum]
        pointVal = PrincetonSPEFile(fileName)[0].astype(numpy.float64)

        # considered region of interest
        if self.conRoi == None:
            self.conRoi   = [1, len(pointVal[0]), 1, len(pointVal)]

        # get considered part of the images
        conDark  =  getAreaSet(darkVal,  self.conRoi)
        conIm    = (getAreaSet(pointVal, self.conRoi) - conDark) / pointMon

        return conIm

    def _XYCorrect(self, xVal, yVal):
        """Correct the miss alignement of the CCD camera

        xVal : measured X-values of the detector
        yVal : measured Y-values of the detector

        return
        xNew : corrected X-values of the detector
        yNew : corrected Y-values of the detector"""

        xNew = np.cos(self.detAngle) * xVal - np.sin(self.detAngle) * yVal
        yNew = np.sin(self.detAngle) * xVal + np.cos(self.detAngle) * yVal

        return xNew, yNew

    def _XY2delgam(self, del0, gam0):

        """
        calculate (delta, gamma) in deg from (x, y), flattend arrays are used
        x, y       : given (x, y)coordinates on CCD
        del0, gam0 : given (delta, gamma) in deg of point0, e.g. center
        x0, y0     : given (x, y) of point0, e.g. center
        """

        # (x, y)-values
        conZ = getAreaXY(self.conRoi)
        x    = np.ravel(conZ[0])
        y    = np.ravel(conZ[1])

        x, y = self._XYCorrect(x, y)

        # detector distance
        detDis = self.detDis
        # pixel distance binded 4x4 -> 80um (micro meter)
        pixDisX = self.detPixSizeX
        pixDisY = self.detPixSizeY
        # (x, y) for point0, take center of CCD
        x0 = self.detX0
        y0 = self.detY0

        # work in grad
        delta = del0 - np.arctan( (y-y0)*pixDisY/detDis )/np.pi*180.0
        gamma = gam0 - np.arctan( (x-x0)*pixDisX/detDis )/np.pi*180.0
    
        return delta, gamma

    def _processOneImage(self, imNum, mode = None):
        """Process one image to (Qx, Qy, Qz, I)

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        # used mode
        if mode == None:
            mode = self.frameMode
        # angle alias
        delta = self.settingAngles[imNum, 0]
        theta = self.settingAngles[imNum, 1]
        chi   = self.settingAngles[imNum, 2]
        phi   = self.settingAngles[imNum, 3]
        mu    = self.settingAngles[imNum, 4]
        gamma = self.settingAngles[imNum, 5]
        
        # intensities of considered image part
        intent = np.ravel( self._readImage(imNum) )

        # (delta, gamma)-values at each pixel
        delPix, gamPix = self._XY2delgam(delta, gamma)        
  
        # diffractometer for angle to q calculations
        scanDiff = Diffractometer()
        scanDiff.setLambda(self.waveLen)
        scanDiff.setAngles(delta = delPix, theta = theta, chi   = chi   ,
                           phi   = phi   , mu    = mu   , gamma = gamPix)
        scanDiff.calc()
        scanDiff.setUbMatrix(self.UBmat)
        if mode == 1:
            Qxyz = scanDiff.getQTheta()
        elif mode == 2:
            Qxyz = scanDiff.getQPhi()
        elif mode == 3:
            Qxyz = scanDiff.getQCart()
        elif mode == 4:
            Qxyz = scanDiff.getQHKL()
        else:
            print 'mode = %s is no proper mode for calculation of (Qx, Qy, Qz)!'
            print 'choose  theta- (1), phi- (2), cartesian- (3) or hkl-frame (4)'

        # out put (Qx, Qy, Qz, I)
        totIm = np.zeros((intent.size, 4))
        totIm[:,:3] = Qxyz
        totIm[:,3]  = intent

        return totIm

    def _calcVecDataSet(self):
        """Calculats the vector data set for the grid points"""

        # aliases
        minBox = self.Qmin
        maxBox = self.Qmax
        dVec   = (maxBox - minBox) / self.dQN
        
        # vector data set of the center of each grid part
        qxVal = np.arange(minBox[0], maxBox[0] - dVec[0]/2, dVec[0]) + dVec[0]/2
        qyVal = np.arange(minBox[1], maxBox[1] - dVec[1]/2, dVec[1]) + dVec[1]/2
        qzVal = np.arange(minBox[2], maxBox[2] - dVec[2]/2, dVec[2]) + dVec[2]/2
        self.qVal = [qxVal, qyVal, qzVal]

    def _calcMax(self):
        """Calculates the position of the maximum as indicies"""
        maxN = self.gridData.argmax()
        ind2 = maxN % self.dQN[2]
        ind1 = maxN / self.dQN[2] % self.dQN[1]
        ind0 = maxN / self.dQN[2] / self.dQN[1]
        self.maxInd = np.array([ind0, ind1, ind2])    

    def _fit1DData(self, xVal, yVal, fitType = None, infoDes = ''):
        """Fit a 1D data set

        xVal    : x-values
        yVal    : y-values
        fitType : peak shape to fit from pyspec fitfuncs, use object default if None, e.g. 'lor2a'
        infoDes : description of the current fit try for output, e.g. 'Qx Line cut of Scan #244'

        returns
        fitRes  : results of the fit, [a1, b1, cen1, width1, area1], [0, 0, 0, 0, 0] if unsuccessful fit"""

        if fitType == None:
            fitType = self.fit1DType

        fitRes = np.zeros(5)
        try:
            f = fit.fit(x=xVal, y=yVal, funcs = [fitfuncs.linear, getattr(fitfuncs, fitType)])
            f.go()
            fitRes = f.result
        except:
            print 'WARNING : %s could not be fitted to %s!' % (infoDes, fitType)

        return fitRes

    def _add1DFits(self, xVals, yVals, axes, fitType = None, infoDes = ''):
        """Tries fits of the 1D data and adds them to the data plot

        xVals   : list with the x-values
        yVals   : list with the y-values
        axes    : list of the axes where the data was plotted to add the fits
        fitType : peak shape to fit from pyspec fitfuncs, use object default if None, e.g. 'lor2a'
        infoDes : description of the current fit try for output, e.g. 'Qx Line cut of Scan #244'

        returns
        allRes  : all results of the fits, [[a1, b1, cen1, width1, area1],...], [0, 0, 0, 0, 0] if unsuccessful fit"""

        allRes = np.zeros((len(xVals), 5))

        # go through the list of 1D data sets
        for i in range(len(xVals)):
            # try to fit the 1D data
            infoDesI = infoDes + ' %d ' % (i)
            fitRes   = self._fit1DData(xVals[i], yVals[i], fitType = fitType, infoDes = infoDesI)
            allRes[i,:] = fitRes
            # if the fit was successful plot it to the data
            if not np.all(fitRes == 0):
                yFit = fitfuncs.linear(xVals[i], fitRes[:2]) + getattr(fitfuncs, fitType)(xVals[i], fitRes[2:])
                axes[i].plot(xVals[i], yFit, '-r')

        return allRes

    #
    # help functions for input / output
    #

    def _makeSetInfo(self):
        """Create the information about the set of images"""
        
        self.opSetInfo = '%s%s' % (self.setName, self.setNum)
        try:
            self.opSetInfo += ' Images: %s\n' % (self.procImSelect)
        except:
            self.opSetInfo += '\n'
        try:
            self.opSetInfo += 'Sample Temperature: %.2f K\t' % (self.temperature)
        except:
           self.opSetInfo += '' 
        try:
            self.opSetInfo += 'Photon Wavelength: %.2f Angst.\t' % (self.waveLen)
        except:
           self.opSetInfo += ''
        try:
            self.opSetInfo += 'Photon Energy: %.2f eV' % (self.energy)
        except:
           self.opSetInfo += ''
    
    def _makeModeInfo(self, mode = None):
        """Create info description of the used frame mode

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        if mode == None:
            mode = self.frameMode

        if mode == 1:
            self.opModeInfo = 'Frame Mode 1: (Qx, Qy, Qz) in theta-frame and (1/Angstrom)' 
        if mode == 2:
            self.opModeInfo = 'Frame Mode 2: (Qx, Qy, Qz) in phi-frame and (1/Angstrom)'
        if mode == 3:
            self.opModeInfo = 'Frame Mode 3: (Qx, Qy, Qz) in cartesian-frame and (1/Angstrom)'
        if mode == 4:
            self.opModeInfo = 'Frame Mode 4: (H, K, L) in hkl-frame and (reciprocal lattice units)'

    def _makeHeaderInfo(self):
        """Create header for output files"""

        if self.projName == '':
            opGap = ''
        else:
            opGap = ' - '
        self.opHeader = '%s%s%s\n%s\n%s' % (self.opTitle, opGap, self.projName, self.opSetInfo, self.opModeInfo)

    def _makeFitInfo1D(self, allRes, fitType = None, fitTitle = None, fitNames = None):
        """Create information output for the fit results in 1D

        allRes   : all results of the fittings, [[a, b, cen, width, area],...]
        fitType  : tpe of the fitting, e.g 'lor2a'
        fitTitle : title for the current fitting process, e.g. '1D Line cuts'
        fitNames : name for each fitting, e.g. ['Qx', 'Qy', 'Qz']"""

        # prepare information
        if fitTitle == None:
            fitTitle = ''
        else:
            fitTitle = fitTitle + ' '
        if fitType == None:
            fitType = self.fitType
        if fitNames == None:
            fitNames = []
            for i in range(allRes.shape[0]):
                fitNames.append('Fit %02d' % i)

        fitInfo = '%s%s\n\t a \t\t b \t\t cen \t\t width \t\t area' % (fitTitle, fitType)
        line    = '\n%s \t ' + 4*'%.5e\t ' + '%.5e'
        for i in range(allRes.shape[0]):
            fitInfo += line % (fitNames[i], allRes[i,0], allRes[i,1], allRes[i,2], allRes[i,3], allRes[i,4])

        return fitInfo

    def _makeGridInfo(self, gData = None, gOccu = None, gOut = None, gMin = None, gMax = None, dg = None):
        """Create information about the grid

        gridData : values at the grid parts (bins)
        gridOccu : occupation no. of the grid parts (bins)
        gridOut  : no. of data points outside the grid
        Qmin     : minimal Q-values of the grid
        Qmax     : maximal Q-values of the grid
        dQN      : no. of grid parts (bins) in the three principal directions"""

        # prepare information, if not given, take it from the object
        if gData == None:
            gData = self.gridData
        if gOccu == None:
            gOccu = self.gridOccu
        if gOut  == None:
            gOut  = self.gridOut
        if gMin == None:
            gMin = self.Qmin
        if gMax == None:
            gMax = self.Qmax
        if dg == None:
            dg = self.dQN
        emptNb = (gOccu == 0).sum()
        
        gridInfo = '\n\n%s sets processed to grid\n' % (self.setEntLabel) + \
            'Grid dimension = %s \t No. of bins in the grid %.2e\n' % (dg, gData.size) + \
            'Data points outside the grid : %.2e \t Bins with zero information : %.2e' % (gOut, emptNb)
        gridInfo += '\n\t min \t\t max \t\t step'
        line = '\n%s \t %.2e \t %.2e \t %.2e'
        for i in range(gMin.size):
            gridInfo += line % (self.qLabel[i], gMin[i], gMax[i], dg[i])

        return gridInfo

    #
    # process part
    #

    def processOneSet(self, procSelect = None, mode = None):
        """Process selcted images of the full set into (Qx, Qy, Qz, I)

        procSelect : list with the images which will be processed, take object default if None
        mode       : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        print '\n%s%d: process to %s' % (self.setName, self.setNum, self.setEntLabel)

        if procSelect == None:
            procSelect = self.procImSelect
        if mode == None:
            mode = self.frameMode

        # prepare size of full dataset and get data of first scan
        procSize  = len(procSelect)
        totFirst = self._processOneImage(procSelect[0])
        imSize   = totFirst.shape[0]
        npts     = procSize * imSize
        totSet   = np.zeros((npts, 4))

        # go through all selected images and get data sets
        j = 0
        k = imSize
        totSet[j:k,:] = totFirst
        for i in procSelect[1:]:
            j = j + imSize
            k = k + imSize
            totSet[j:k, :] = self._processOneImage(i, mode = mode)

        # for info file
        self.opProcInfo += '\n\nImage Set processed to %.2e %s sets' % (totSet.shape[0], self.setEntLabel)

        return totSet

    def makeGridData(self, procSelect = None, mode = None):
        """Grid the data set into a cuboid
        Size of the cuboid : Qmin, Qmax = [Qx, Qy, Qz]_min, max
        Number of parts    : dQN = [Nqx, Nqy, Nqz]

        procSelect : list with the images which will be processed, take object default if None
        mode       : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None

        returns
        gridData : values at the grid parts (bins)
        gridOccu : occupation no. of the grid parts (bins)
        gridOut  : no. of the data points outside the grid"""

        if procSelect == None:
            procSelect = self.procImSelect
        if mode == None:
            mode = self.frameMode

        totData = self.processOneSet(procSelect = procSelect, mode = mode)

        # prepare min, max,...
        if self.Qmin == None:
            self.Qmin = np.array([ totData[:,0].min(), totData[:,1].min(), totData[:,2].min() ])
        if self.Qmax == None:
            self.Qmax = np.array([ totData[:,0].max(), totData[:,1].max(), totData[:,2].max() ])
        if self.dQN  == None:
            self.dQN = [100, 100, 100]

        # use alias for grid options
        Qmin = self.Qmin
        Qmax = self.Qmax
        dQN  = self.dQN

        # 3D grid of the data set 
        gridData, gridOccu, gridOut = gridder.grid3d(totData,Qmin, Qmax, dQN, norm = 1)
        emptNb = (gridOccu == 0).sum()
        if gridOut != 0:
            print "Warning : There are %.2e points outside the grid (%.2e bins in the grid)" % (gridOut, gridData.size)
        if emptNb:
            print "Warning : There are %.2e values zero in the grid" % emptNb

        # mask the gridded data set
        #gridData = np.ma.array(gridData / gridOccu, mask = (gridOccu == 0))
        # store intensity, occupation and no. of outside data points of the grid
        self.gridData = gridData
        self.gridOccu = gridOccu
        self.gridOut  = gridOut

        # calculated the corresponding vectors and maximum intensity position of the grid
        self._calcVecDataSet()
        self._calcMax()
        
        # for info file
        self.opProcInfo += self._makeGridInfo()

        return gridData, gridOccu, gridOut

    def get1DSum(self):
        """1D Lines of the grid data and occupations by summing in the other directions

        returns
        gridData1DSum : intensity set  in the order Qx, Qy, Qz as list
        gridOccu1DSum : occupation no. in the order Qx, Qy, Qz as list"""

        gridData1DSum = [self.gridData.sum(1).sum(1),
                         self.gridData.sum(0).sum(1),
                         self.gridData.sum(0).sum(0)]
        gridOccu1DSum = [self.gridOccu.sum(1).sum(1),
                         self.gridOccu.sum(0).sum(1),
                         self.gridOccu.sum(0).sum(0)]

        return gridData1DSum, gridOccu1DSum

    def get2DSum(self):
        """2D Areas of the grid data and occupations by summing in the other direction

        returns
        gridData2DSum : intensity set  in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list
        gridOccu2DSum : occupation no. in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        gridData2DSum = [self.gridData.sum(0), self.gridData.sum(1), self.gridData.sum(2)]
        gridOccu2DSum = [self.gridOccu.sum(0), self.gridOccu.sum(1), self.gridOccu.sum(2)]

        return gridData2DSum, gridOccu2DSum
    
    def get1DCut(self):
        """1D Lines of the grid data and occupations at the position of the maximum intensity

        returns
        gridData1DCut : intensity set  in the order Qx, Qy, Qz as list
        gridOccu1DCut : occupation no. in the order Qx, Qy, Qz as list"""

        gridData1DCut = [self.gridData[:,self.maxInd[1],self.maxInd[2]],
                         self.gridData[self.maxInd[0],:,self.maxInd[2]],
                         self.gridData[self.maxInd[0],self.maxInd[1],:]]
        gridOccu1DCut = [self.gridOccu[:,self.maxInd[1],self.maxInd[2]],
                         self.gridOccu[self.maxInd[0],:,self.maxInd[2]],
                         self.gridOccu[self.maxInd[0],self.maxInd[1],:]]

        return gridData1DCut, gridOccu1DCut
    
    def get2DCut(self):
        """2D Areas of the grid data and occupations at the position of the maximum intensity

        returns
        gridData2DCut : intensity set  in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list
        gridOccu2DCut : occupation no. in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        gridData2DCut = [self.gridData[self.maxInd[0],:,:],
                         self.gridData[:,self.maxInd[1],:],
                         self.gridData[:,:,self.maxInd[2]]]
        gridOccu2DCut = [self.gridOccu[self.maxInd[0],:,:],
                         self.gridOccu[:,self.maxInd[1],:],
                         self.gridOccu[:,:,self.maxInd[2]]]

        return gridData2DCut, gridOccu2DCut

    def get1DCutAv(self):
        """1D averaged Lines of the grid data and occupations at the position of the maximum
        intensity and its eight neighbored lines 

        returns
        gridData1DCutAv : intensity set  in the order Qx, Qy, Qz as list
        gridOccu1DCutAv : occupation no. in the order Qx, Qy, Qz as list"""

        # initialize with correct size as zeros
        gridData1DCutAv = [np.zeros(self.dQN[0]),np.zeros(self.dQN[1]),np.zeros(self.dQN[2])]
        gridOccu1DCutAv = [np.zeros(self.dQN[0]),np.zeros(self.dQN[1]),np.zeros(self.dQN[2])]

        #print self.dQN[0]
        # go through the neighbors
        for i in range(3):
            for j in range(3):
                gridData1DCutAv[0] += self.gridData[:,self.maxInd[1]+i-1,self.maxInd[2]+j-1]/9.0
                gridData1DCutAv[1] += self.gridData[self.maxInd[0]+i-1,:,self.maxInd[2]+j-1]/9.0
                gridData1DCutAv[2] += self.gridData[self.maxInd[0]+i-1,self.maxInd[1]+j-1,:]/9.0
                gridOccu1DCutAv[0] += self.gridOccu[:,self.maxInd[1]+i-1,self.maxInd[2]+j-1]/9.0
                gridOccu1DCutAv[1] += self.gridOccu[self.maxInd[0]+i-1,:,self.maxInd[2]+j-1]/9.0
                gridOccu1DCutAv[2] += self.gridOccu[self.maxInd[0]+i-1,self.maxInd[1]+j-1,:]/9.0
        
        return gridData1DCutAv, gridOccu1DCutAv

    def get2DCutAv(self):
        """2D average Areas of the grid data and occupations at the position of the maximum
        intensity and the their two neighbors

        returns
        gridData2DCutAv : intensity set  in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list
        gridOccu2DCutAv : occupation no. in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        # initialize with correct size as zeros
        gridData2DCutAv = [np.array(np.meshgrid(np.zeros(self.dQN[2]),np.zeros(self.dQN[1])))[0],
                           np.array(np.meshgrid(np.zeros(self.dQN[2]),np.zeros(self.dQN[0])))[0],
                           np.array(np.meshgrid(np.zeros(self.dQN[1]),np.zeros(self.dQN[0])))[0]]
        gridOccu2DCutAv = [np.array(np.meshgrid(np.zeros(self.dQN[2]),np.zeros(self.dQN[1])))[0],
                           np.array(np.meshgrid(np.zeros(self.dQN[2]),np.zeros(self.dQN[0])))[0],
                           np.array(np.meshgrid(np.zeros(self.dQN[1]),np.zeros(self.dQN[0])))[0]]
        print gridData2DCutAv[0].shape
        # go through the neighbors
        for i in range(3):
            gridData2DCutAv[0] += self.gridData[self.maxInd[0]+i-1,:,:]/3.0
            gridData2DCutAv[1] += self.gridData[:,self.maxInd[1]+i-1,:]/3.0
            gridData2DCutAv[2] += self.gridData[:,:,self.maxInd[2]+i-1]/3.0
            gridOccu2DCutAv[0] += self.gridOccu[self.maxInd[0]+i-1,:,:]/3.0
            gridOccu2DCutAv[1] += self.gridOccu[:,self.maxInd[1]+i-1,:]/3.0
            gridOccu2DCutAv[2] += self.gridOccu[:,:,self.maxInd[2]+i-1]/3.0
        return gridData2DCutAv, gridOccu2DCutAv

    def getIntIntensity(self):
        """Get the integrated intensity of the peak by summing over all grid parts (bins)

        returns
        intInten : integrated intensity"""

        intInten = self.gridData.sum()

        return intInten
    
    #
    # plot part
    #

    def plotGrid1DSum(self):
        """Plots the 1D Lines of the data grid summed over the other dimensions

        retrurns
        fig1   : plt.figure object of the plotting window
        allax1 : list of plt.axes objects which carry the figures
        allRes : all results of the fits, [[a1, b1, cen1, width1, area1],...], [0, 0, 0, 0, 0] if unsuccessful fit"""

        # results for fit of 1D data, None if no fitting
        allRes = None

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag1D = 7)
        gridPlot.setLogFlags(flag1D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot1DAxes(self.qVal, self.axesLabels)
        gridData1DSum, gridOccu1DSum = self.get1DSum()
        gridPlot.setPlot1DData(gridData1DSum, gridOccu1DSum,
                               plotTitle = '1D Lines, over other directions is summed')
        # plot, get figure and axes back
        fig1, allax1 = gridPlot.plot1DData()
        allRes = np.zeros((3,5))
        # try to fit the 1D data
        if self.fit1D:
            infoDes = '%s%s 1D Line summed' % (self.setName, self.setNum)
            allRes  = self._add1DFits(self.qVal, gridData1DSum, axes = allax1[:3], fitType = self.fit1DType, infoDes = infoDes)
            # for info file
            self.opProcInfo += '\n\n' + self._makeFitInfo1D(allRes, fitType = None, fitTitle = '1D Line summed', fitNames = self.qLabel)
        
        return fig1, allax1, allRes

    def plotGrid2DSum(self):
        """Plots the 2D Areas of the data grid summed over the other dimension

        retrurns
        fig2   : plt.figure object of the plotting window
        allax2 : list of plt.axes objects which carry the figures"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag2D = 7)
        gridPlot.setLogFlags(flag2D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot2DAxes([self.Qmin[2], self.Qmin[2], self.Qmin[1]], [self.Qmax[2], self.Qmax[2], self.Qmax[1]],
                               [self.Qmin[1], self.Qmin[0], self.Qmin[0]], [self.Qmax[1], self.Qmax[0], self.Qmax[0]],
                               [self.axesLabels[2], self.axesLabels[2], self.axesLabels[1]],
                               [self.axesLabels[1], self.axesLabels[0], self.axesLabels[0]])
        gridData2DSum, gridOccu2DSum = self.get2DSum()
        for i in range(3):
            gridData2DSum[i] = np.ma.array(gridData2DSum[i], mask = (gridOccu2DSum[i] == 0))
        gridPlot.setPlot2DData(gridData2DSum, gridOccu2DSum,
                               plotTitle = '2D Areas, over other direction is summed')
        # plot, get figure and axes back
        fig2, allax2 = gridPlot.plot2DData()

        return fig2, allax2

    def plotGrid1DCut(self):
        """Plots the 1D Lines of the data grid summed over the other dimensions

        retrurns
        fig1   : plt.figure object of the plotting window
        allax1 : list of plt.axes objects which carry the figures
        allRes : all results of the fits, [[a1, b1, cen1, width1, area1],...], [0, 0, 0, 0, 0] if unsuccessful fit"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag1D = 7)
        gridPlot.setLogFlags(flag1D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot1DAxes(self.qVal, self.axesLabels)
        gridData1DCut, gridOccu1DCut = self.get1DCut()
        gridPlot.setPlot1DData(gridData1DCut, gridOccu1DCut,
                               plotTitle = '1D Line Cuts at Maximum Position')
        # plot, get figure and axes back
        fig1, allax1 = gridPlot.plot1DData()
        allRes = np.zeros((3,5))
        # try to fit the 1D data
        if self.fit1D:
            infoDes = '%s%s 1D Line cut' % (self.setName, self.setNum)
            allRes  = self._add1DFits(self.qVal, gridData1DCut, axes = allax1[:3], fitType = self.fit1DType, infoDes = infoDes)
            # for info file
            self.opProcInfo += '\n\n' + self._makeFitInfo1D(allRes, fitType = None, fitTitle = '1D Line cut', fitNames = self.qLabel)

        return fig1, allax1, allRes

    def plotGrid2DCut(self):
        """Plots the 2D Areas of the data grid summed over the other dimension

        retrurns
        fig2   : plt.figure object of the plotting window
        allax2 : list of plt.axes objects which carry the figures"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag2D = 7)
        gridPlot.setLogFlags(flag2D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot2DAxes([self.Qmin[2], self.Qmin[2], self.Qmin[1]], [self.Qmax[2], self.Qmax[2], self.Qmax[1]],
                               [self.Qmin[1], self.Qmin[0], self.Qmin[0]], [self.Qmax[1], self.Qmax[0], self.Qmax[0]],
                               [self.axesLabels[2], self.axesLabels[2], self.axesLabels[1]],
                               [self.axesLabels[1], self.axesLabels[0], self.axesLabels[0]])
        gridData2DCut, gridOccu2DCut = self.get2DCut()
        for i in range(3):
            gridData2DCut[i] = np.ma.array(gridData2DCut[i], mask = (gridOccu2DCut[i] == 0))
        gridPlot.setPlot2DData(gridData2DCut, gridOccu2DCut,
                               plotTitle = '2D Area Cuts at Maximum Position')
        
        # plot, get figure and axes back
        fig2, allax2 = gridPlot.plot2DData()

        return fig2, allax2

    def plotGrid1DCutAv(self):
        """Plots the 1D Lines of the data grid summed over the other dimensions

        retrurns
        fig1   : plt.figure object of the plotting window
        allax1 : list of plt.axes objects which carry the figures
        allRes : all results of the fits, [[a1, b1, cen1, width1, area1],...], [0, 0, 0, 0, 0] if unsuccessful fit"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag1D = 7)
        gridPlot.setLogFlags(flag1D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot1DAxes(self.qVal, self.axesLabels)
        gridData1DCutAv, gridOccu1DCutAv = self.get1DCutAv()
        gridPlot.setPlot1DData(gridData1DCutAv, gridOccu1DCutAv,
                               plotTitle = '1D Line Cuts at Maximum Position and 8 Neighbors Averaged')
        # plot, get figure and axes back
        fig1, allax1 = gridPlot.plot1DData()
        allRes = np.zeros((3,5))
        # try to fit the 1D data
        if self.fit1D:
            infoDes = '%s%s 1D Line cuts average' % (self.setName, self.setNum)
            allRes  = self._add1DFits(self.qVal, gridData1DCutAv, axes = allax1[:3], fitType = self.fit1DType, infoDes = infoDes)
            # for info file
            self.opProcInfo += '\n\n' + self._makeFitInfo1D(allRes, fitType = None, fitTitle = '1D Line cut average', fitNames = self.qLabel)

        return fig1, allax1, allRes

    def plotGrid2DCutAv(self):
        """Plots the 2D Areas of the data grid summed over the other dimension

        retrurns
        fig2   : plt.figure object of the plotting window
        allax2 : list of plt.axes objects which carry the figures"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag2D = 7)
        gridPlot.setLogFlags(flag2D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot2DAxes([self.Qmin[2], self.Qmin[2], self.Qmin[1]], [self.Qmax[2], self.Qmax[2], self.Qmax[1]],
                               [self.Qmin[1], self.Qmin[0], self.Qmin[0]], [self.Qmax[1], self.Qmax[0], self.Qmax[0]],
                               [self.axesLabels[2], self.axesLabels[2], self.axesLabels[1]],
                               [self.axesLabels[1], self.axesLabels[0], self.axesLabels[0]])
        gridData2DCutAv, gridOccu2DCutAv = self.get2DCutAv()
        for i in range(3):
            gridData2DCutAv[i] = np.ma.array(gridData2DCutAv[i], mask = (gridOccu2DCutAv[i] == 0))
        gridPlot.setPlot2DData(gridData2DCutAv, gridOccu2DCutAv,
                               plotTitle = '2D Area Cuts at Maximum Position and 2 Neighbors Averaged')
        
        # plot, get figure and axes back
        fig2, allax2 = gridPlot.plot2DData()

        return fig2, allax2
    
    def plotImages(self, images = None):
        """Plots the selcted images

        images : selction of the images for plotting, use object default if None

        retrurns
        allfig : list of plt.figure objects of the plotting windows
        allax  : list of plt.axes objects which carry the figures"""
        
        # prepare plots
        plotImNum = self.plotImHor * self.plotImVer
        j = 0
        allfig = []
        allax  = []
        plotImTitle  = '%s%s' % (self.setName, self.setNum)
        # read in first image to get conRoi if not set
        if self.conRoi == None:
            self._readImage(0)
        plotImExtent = [self.conRoi[0], self.conRoi[0] + self.conRoi[1],
                        self.conRoi[2], self.conRoi[2] + self.conRoi[3]]

        if images is None:
            images = self.plotImSelect

        # go through images numbers which should be plotted
        for i in images:

            # label for y-axis
            yLabel = 'image # %d' % (i)

            if j%plotImNum == 0:
                # prepare plot window
                fig = plt.figure(figsize = self._defaultFigureSize)
                fig.suptitle(plotImTitle, fontsize = 24)
                allfig.append(fig)

            # new subplot
            ax = plt.subplot(self.plotImVer, self.plotImHor, j%plotImNum+1)
            plt.subplots_adjust(hspace = 0.4)
            allax.append(ax)
            ax.imshow(self._readImage(i), extent = plotImExtent)

            # show the image with number as y-label    
            ax.set_ylabel(yLabel, fontsize = 18)
                     
            # increment the plot image counter
            j += 1

        return allfig, allax

    def plotAll(self):
        """Plots 1D/2D sums and cuts"""
        self.plotGrid1DSum()
        self.plotGrid2DSum()
        self.plotGrid1DCut()
        self.plotGrid2DCut()
        self.plotGrid1DCutAv()
        self.plotGrid2DCutAv()

    #
    # input / output part
    #

    def makeInfo(self):
        """Create the information about the current processing"""

        self._makeHeaderInfo()
        curInfo = '%s%s' % (self.opHeader, self.opProcInfo)

        return curInfo

    def writeInfoFile(self, outFile = None):
        """Write information about the current processing into a file

        outFile : path and file name of the output file"""

        if outFile == None:
            outFile = self.infoFile

        out = file(outFile, 'w')
        out.write(self.makeInfo())
        out.close()
        
####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    sf   = spec.SpecDataFile('/home/tardis/spartzsch/2010_09_X1A2/ymn2o5_sep10_1', 
			     ccdpath = '/mounts/davros/nasshare/images/sept10')
    scan = sf[244]

    testData = ImageProcessor()

    testData.setDetectorAngle(-1.24)
    testData.setBins(4, 4)
    testData.setSpecScan(scan)
    #testData.setConRoi([1, 325, 1, 335])
    testData.setFrameMode(1)
    testData.setGridOptions(Qmin = None, Qmax = None, dQN = [90, 160, 30])
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [200, 400, 100])
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [100, 100, 100])

    #testData.setPlotIm(plotImSelect = [40], plotImHor = 4, plotImVer = 3)
    #testData.plotImages()
    
    #imSet  = testData.processOneSet(procSelect = [40])
    #totSet = testData.processOneSet()
    testData.makeGridData(procSelect = [40])
    #testData.makeGridData()

    #print 'Peak integrated intensity : %.2e' % (testData.getIntIntensity())
    #lineData, lineOccu = testData.get1DCut()
    #areaData, areaOccu = testData.get2DCut()

    # plot options
    #testData.setAxesLabels([ur"Qx (\u00c5$^{-1}$)", ur"Qy (\u00c5$^{-1}$)", ur"Qz (\u00c5$^{-1}$)"])
    #testData.setAxesLabels(['H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'])
    testData.setPlotFlags(7, 7)
    testData.setLogFlags(0, 3)
    testData.setFit1D(False)
    testData.setHistBin(50)
 
    testData.plotGrid1DSum()
    #testData.plotGrid2DSum()
    testData.plotGrid1DCut()
    #testData.plotGrid2DCut()
    testData.plotGrid1DCutAv()
    #testData.plotGrid2DCutAv()
    #testData.plotAll()
 
    # test of input output file handling

    print '\n\n'
    print testData.makeInfo()
    #testData.writeInfoFile(outFile = 'infoFile.dat')

    plt.show()

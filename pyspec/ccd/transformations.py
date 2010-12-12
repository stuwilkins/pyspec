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
import gc
#gc.set_debug(gc.DEBUG_LEAK | gc.DEBUG_STATS)
import time
import sys
import numpy as np
from   scipy.optimize import leastsq
import matplotlib.pyplot as plt
from   pyspec import fit, fitfuncs
from   pyspec.diffractometer import Diffractometer
try:
    from   pyspec.ccd.files import *
except:
    passfrom   pyspec.ccd.plotter import PlotImages, PlotGrid
try:
    import pyspec.ccd.ctrans as ctrans
except:
    try:
        import ctrans
    except:
        pass

from ConfigParser import ConfigParser

try:
    import princeton
except:
    print "No princeton"
    pass

#gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_LEAK)
gc.enable()

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>, Sven Partzsch <SvenPartzsch@gmx.de>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

class CCDParamsConfigParser(ConfigParser):
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

def get3DMesh(xVal, yVal, zVal):
    """make a 3D meshgrid from the given values

    xVal : x-values, [x0, x1, x(Nx-1)] (list or np.array)
    yVal : y-values, [y0, y1, y(Ny-1)] (list or np.array)
    zVal : z-values, [z0, z1, z(Nz-1)] (list or np.array)

    returns xGrid, yGrid, zGrid as Nz x Ny x Nx np.arrays"""

    # work first with lists
    xVal = list(xVal)
    yVal = list(yVal)
    zVal = list(zVal)

    # dimensions
    Nx = len(xVal)
    Ny = len(yVal)
    Nz = len(zVal)

    # grid results
    xGrid = np.array(xVal * Ny * Nz).reshape(Nz, Ny, Nx)
    yGrid = np.array(yVal * Nz * Nx).reshape(Nz, Nx, Ny).transpose(0,2,1)
    zGrid = np.array(zVal * Nx * Ny).reshape(Nx, Ny, Nz).transpose(2,1,0)

    return xGrid, yGrid, zGrid

#
# FileProcessor Class
#

class FileProcessor():
    """FileProcessor Class

    This class processes CCD files and returns a numpy array of 
    float values for all the images. Images can be saved to disk
    for faster processing in the future."""

    def __init__(self, filenames = None, 
                 darkfilenames = None, 
                 norm = None, format = 'SPE',
                 spec = None, mon = 'Monitor'):
        """Initialize the class

        filenames      : list of filenames for all images
        darkfilenames  : list of filenames for the dark images
        spec           : SpecScan object from which to obtain data"""
        self._format = format
        self.filenames = filenames
        self.darkfilenames = darkfilenames
        self.normData = None

        if spec is not None:
            self.setFromSpec(spec)
        if norm is not None:
            self.normData = norm
        self.bgndParams = np.array([])

    def setFromSpec(self, scan, mon = 'Monitor'):
        """Set the filenames from a SpecScan instance

        scan    : SpecScan instance
        norm    : spec counter to normalize against"""
        self.filenames = scan.getCCDFilenames()
        self.darkfilenames = scan.getCCDFilenames(dark = 1)
        self.normData = scan.values[mon]

    def _processBgnd(self, maskroi = None, mask = None):

        bgndfunc = self._residualsLinear

        x, y = np.meshgrid(range(self.images.shape[2]), 
                           range(self.images.shape[1]))

        if mask is None:
            mask = np.ravel(np.ones(self.images.shape[1:]) == 1)
            if maskroi is not None:
                for m in maskroi:
                    xmask = (x >= m[0]) & (x <= (m[0] + m[2]))
                    ymask = (y >= m[1]) & (y <= (m[1] + m[3]))
                    mask = mask & (np.ravel((xmask & ymask)) == False)

        _x = np.ravel(x)[mask]
        _y = np.ravel(y)[mask]
        allplsq = np.array([])
        for i in range(self.images.shape[0]):
            guess = [1e-3, 1e-3, self.images[i].mean()]
            z = np.ravel(self.images[i])[mask]
            plsq = leastsq(bgndfunc, guess, args = (z, _y, _x))
            self.images[i] = self.images[i] - x*plsq[0][0] - y*plsq[0][1] - plsq[0][2]
            allplsq = np.concatenate((allplsq, plsq[0]))
        allplsq = allplsq.reshape(-1, 3)
        self.bgndParams = allplsq            

    def _residualsLinear(self, p, z, x, y):
        mx, my, c = p
        err = z - (mx * x) - (my * y) - c
        return err

    def processBgnd(self, maskroi = None, mask = None):
        """Process the background by fitting each image

        maskroi : list of 4 element tuples or lists for ROI to mask
        mask    : (n x m) array for the background mask"""
        print "---- Subtracting background from images"
        self._processBgnd(maskroi = maskroi, mask = mask)
        print "---- Done."

    def process(self, dark = True, norm = True, dtype = np.float):
        """Read in images and process into 3D array.
        
        dark : if True, subtract the dark images from the data"""

        images = []
        darkimages = []

        if norm:
            if self.normData is None:
                normData = np.ones(len(self.filenames))
                print "XXXX No normalization data found"
            else:
                normData = self.normData
                print "---- Normalizing data."
        else:
            normData = np.ones(len(self.filenames))

        if dark:
            print "---- Correcting for dark current"

        for i, (iname, diname, normVal) in enumerate(zip(self.filenames, self.darkfilenames, normData)):
            if type(iname) == list:
                _images = []
                _darkimages = []
                for j, (_in, _din) in enumerate(zip(iname, diname)):
                    if (j % 50) == 0:
                        print "---- Reading image %-3d of %-3d (sub image %-3d of %-3d)     \r" % (i + 1, len(self.filenames),
                                                                                                   j + 1, len(iname)),
                        sys.stdout.flush()
                    image = self._getRawImage(_in).astype(dtype)
                    _images.append(image)
                    if os.path.exists(_din):
                        darkimage =  self._getRawImage(_din).astype(dtype)
                        _darkimages.append(darkimage)
                        if (j % 50) == 0:
                            print "---- Reading dark image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.filenames),
                                                                                                       j + 1, len(iname)),
                            sys.stdout.flush()
                image = np.array(_images).sum(0)
                if len(_darkimages):
                    darkimage = np.array(_darkimages).sum(0)
                    darkimages.append(darkimage)
            else:
                image = self._getRawImage(iname).astype(dtype)
                if os.path.exists(diname):
                    darkimage =  self._getRawImage(diname).astype(dtype)
                    darkimages.append(darkimage)
                print "---- Reading image %-3d of %-3d\r" % (i, len(self.filenames)),
            if dark:
                if len(darkimages):
                    image = image - darkimages[-1]
                else:
                    print "XXXX Unable to dark currect correct. No Image found"
            
            if norm:
                image = image / normVal
            images.append(image)
        print "---- Processed %d images (%d dark images)" % (len(images), len(darkimages))

        self.images = np.array(images)

        print "---- Done. Array size %s (%.3f Mb)." % (str(self.images.shape), 
                                                       self.images.nbytes / 1024**2) 

    def _getRawImage(self, iname):
        """Read raw image"""
        if self._format == 'SPE':
            return PrincetonSPEFile(iname).getBinnedData()
        else:
            raise Exception("Unknown file format \"%s.\"" % self._format)

    def getImage(self, n = None):
        """Return the image data"""
        if n is None:
            return self.images
        else:
            return self.images[n]

    def save(self, filename, compressed = False):
        """Save an image sequance to a numpy binary file
        
        filename   : filename to save to.
        compressed : if True save as compressed file"""
        
        obj = {'images'   : self.images,
               'normData' : self.normData}

        np.savez(filename, **obj)
        print "*** Image saved to %s" % filename

    def load(self, filename):
        """Load images from previous save operation

        filename : filename of images"""
        
        print "Loading image from %s" % filename
        print np.load(filename)

#
# image processor class
#

class ImageProcessor():
    """Image Processor class

    This class provides the processing of single and sets of CCD-images
    in more detail: each pixel is transformed into reciprocal space
    the set of reciprocal vectors and intensities is gridded on a regular cuboid
    the needed informations can be provided by a spec scan from the pyspec package"""

    def __init__(self, fP, configfile = None, ccdname = 'CCD',
                 spec = None):
        """Initialize the image processor

        fP         : file processor object for getting the images
        configfile : 
        ccdname    :
        spec       : """
        
        # file processor to provied the needed images
        self.fileProcessor = fP

        self.ccdName = ccdname
        #self.readConfigFile(configfile)
        
        # set parameters to configure the CCD setup
        # detector distance 30cm and detector pixel size 20um
        self.setDetectorProp(0.020, 0.020, 1300, 1340, 650.0, 670.0)
        self.setDetectorPos(300, 0.0)
        # image treatment
        self.setName     = 'Set #'
        self.setNum      = 1
        self.conRoi      = None
        self.setFrameMode(1)
        # gridder options
        self.setGridOptions()
        # grid background subtraction
        self.setGridBackOptions()
        # cut indicies and position
        self.cutInd = None
        self.setCutPos(cutMode = 'max')
        # 1D fit options
        self.fit1D       = False
        self.fitType     = 'lor2a'
        # output information
        self.opTitle     = 'Image Processing'
        self._makeSetInfo()
        self._makeModeInfo()
        self.opProcInfo  = ''
        # Define variables
        #self.fileProcessor = None
        self.processMode = 'fast'
        self.totSet = None

    def readConfigFile(self, filename):
        config = CCDParamsConfigParser()
        config.read(filename)
        config.get(self.ccdname, '', vars = {})
        
    #
    # set and get part
    #

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
        """Get properties of used detector returned as a tuple

        detPixSizeX : detector pixel size in detector X-direction (float in mm)
        detPixSizeY : detector pixel size in detector Y-direction (float in mm)
        detSizeX    : detector no. of pixels (size) in detector X-direction (integer)
        detSizeY    : detector no. of pixels (size) in detector Y-direction (integer)
        detX0       : detector X-coordinate of center for reference gamma-value (float)
        detY0       : detector Y-coordinate of center for reference delta-value (float)"""

        return self.detPixSizeX, self.detPixSizeY, self.detSizeX, self.detSizeY, self.detX0, self.detY0

    def setDetectorPos(self, detDis = 300.0, detAng = 0.0):
        """Set the detector position

        detDis : detector distance (float in mm)
        detAng : detector miss alignement angle (float in deg)"""
        
        self.detDis = detDis
        self.detAngle = detAng

    def getDetectorPos(self, detDis = 30.0, detAng = 0.0):
        """Get the detector position

        detDis : detector distance (float in mm)
        detAng : detector miss alignement angle (float in deg)"""
        
        return self.detDis, self.detAngle    

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

    def setSetSettings(self, waveLen, imFilesNames, darkFileNames, settingAngles, intentNorm, UBmat, setName, setNum, setSize):
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

        self.waveLen       = waveLen
        self.imFileNames   = imFilesNames
        self.darkFileNames = darkFileNames
        self.settingAngles = settingAngles
        self.intentNorm    = intentNorm
        self.UBmat         = UBmat
        self.setName       = setName
        self.setNum        = setNum
        self.setSize       = setSize

    def getSetSettings(self):
        """Get the settings for the set 

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

        return self.waveLen, self.imFilesNames, self.darkFilesNames, self.settingAngles, self.intentNorm, self.UBmat, self.setName, self.setNum, self.setSize

    def setSpecScan(self, conScan):
        """Set the settings for the set from the considered pyspec scan object

        conScan : pyspec scan object which contains the needed information

        The set settings are:
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        imFileNames   : filenames for each image
        darkFileNames : filenames of the dark images
        settingAngles : setting angles of the diffractometer at each image (data point)
        intentNorm    : normalization factor (division) for each image
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        setName       : name of the considered set, e.g. 'Scan #' in the spec case
        setNum        : no. to determine the set, e.g. 244 in the spec case
        setSize       : no. of images in the set, e.g. 81

        The file Processor processes the corresponding images"""

        self.conScan = conScan
        self._setBySpec()
        self.fileProcessor.setFromSpec(conScan)
        self.fileProcessor.process()
   
    def getSpecScan(self):
        """Get the pyspec scan object which was used for the set settings
        
        returns pyspec scan object which contains the needed information"""

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

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame)
        mode : 'theta', 'phi', 'cart', 'hkl'"""

        if isinstance(mode, str):
            if mode == 'theta':
                self.frameMode = 1
            elif mode == 'phi':
                self.frameMode = 2
            elif mode == 'cart':
                self.frameMode = 3
            elif mode == 'hkl':
                self.frameMode = 4
            else:
                raise Exception("Unknown mode %s." % mode)

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
        
        self.Qmin = np.array(Qmin)
        self.Qmax = np.array(Qmax)
        self.dQN  = np.array(dQN)

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

    def getGridMesh(self):
        """Return the grid vectors as a mesh.

        Returns the grid x,y,z coordinates as 3 3D arrays of values"""
        
        grid = np.mgrid[0:self.dQN[0], 0:self.dQN[1], 0:self.dQN[2]]
        r = (self.Qmax - self.Qmin) / self.dQN

        X = grid[0] * r[0] + self.Qmin[0]
        Y = grid[1] * r[1] + self.Qmin[1]
        Z = grid[2] * r[2] + self.Qmin[2]
        
        return X, Y, Z

    def setGridBackOptions(self, backSub = False, backMaskBox = None, backType = 'threeDLin'):
        """Set the masked box for the background subtraction of the grid

        backSub     : substract background from grid if True
        backMaskBox : [xMin, yMin, zMin, xMax, yMax, zMax]
        backType    : type of background subtraction, 'mean', 'threeDLin'"""

        self.backSub     = backSub
        self.backMaskBox = backMaskBox
        self.backType    = backType

    def getGridBackOptions(self):
        """Get the masked box for the background subtraction of the grid

        backSub     : substract background from grid if True
        backMaskBox : [xMin, yMin, zMin, xMax, yMax, zMax]
        backType    : type of background subtraction, 'mean', 'threeDLin'"""

        return self.backSub, self.backMaskBox, self.backType

    def setIntegrationOptions(self, intMaskBox = None):
        """Set the options for the integration of intensity and background

        intMaskBox : [xmin, ymin, zmin, xmax, ymax, zmax] as np.array"""

        self.intMaskBox = intMaskBox

    def getIntegrationOptions(self):
        """Get the options for the integration of intensity and background                                      
        
        intMaskBox : [xmin, ymin, zmin, xmax, ymax, zmax] as np.array"""

        return self.intMaskBox

    def setCutInd(self, cutInd):
        """Set the cut indicies

        cutInd : cut indices, [nx, ny, nz]
                 if None, indicies of maximum is taken"""

        self.cutInd = cutInd

    def getCutInd(self):
        """Get the cut indicies

        cutInd : cut indices, [nx, ny, nz]
                 if None, indicies of maximum is taken"""

        return self.cutInd
    
    def setCutPos(self, cutPos = None, cutMode = 'fix'):
        """Set the cut position and cut mode

        cutPos  : cut position, [qx, qy, qz]
        cutMode : 'fix' for fixed setting of 'max' for cut at positon of maximum"""

        self.cutPos  = np.array(cutPos)
        self.cutMode = cutMode

    def getCutPos(self):
        """Get the cut position and cut mode

        cutPos  : cut position, [qx, qy, qz]
        cutMode : 'fix' for fixed setting of 'max' for cut at positon of maximum"""

        return self.cutPos, self.cutMode

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
   
    #
    # get set functions for input output
    #

    def setFileProcessor(self, fp = None):
        self.fileProcessor = fp

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
        else:
            self.qLabel = ['H', 'K', 'L']
            self.setEntLabel = '(H, K, L, I)'

    def _readImage(self, imNum):
        """Read in the considered region of interest of the image
        dark image subtraction and normalization by ring current"""

        #if imNum % 10 == 0:
        #    print '%s%d image #%03d: read image' % (self.setName, self.setNum, imNum)

        # get dark image (first CCD image)
        #if self.darkVal is None:
        #    self.darkVal  = PrincetonSPEFile(self.darkFileNames[0])[0].astype(numpy.float)
        
        # get file name and read image of data point
        #fileName = self.imFileNames[imNum]
        #pointVal = PrincetonSPEFile(fileName)[0].astype(numpy.float) - self.darkVal

        pointVal = self.fileProcessor.getImage()[imNum] 

        # considered region of interest
        if self.conRoi == None:
            self.conRoi   = [1, pointVal.shape[1], 1, pointVal.shape[0]]

        # get considered part of the image
        pointMon = self.intentNorm[imNum]
        conIm    = getAreaSet(pointVal, self.conRoi) / pointMon

        return conIm

    def _XYCorrect(self, xVal, yVal):
        """Correct the miss alignement of the CCD camera

        xVal : measured X-values of the detector
        yVal : measured Y-values of the detector

        return
        xNew : corrected X-values of the detector
        yNew : corrected Y-values of the detector"""

        # detetoc angle in rad
        detAn = self.detAngle /180.0 * np.pi

        xNew = np.cos(detAn) * xVal - np.sin(detAn) * yVal
        yNew = np.sin(detAn) * xVal + np.cos(detAn) * yVal

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

    def _processOneImage(self, outArray, imNum, mode = None):
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
        outArray[:,:3] = Qxyz
        outArray[:,3]  = intent
        #del Qxyz 
        #del intent

    def _calcVecDataSet(self):
        """Calculats the vector data set for the grid points

        stores:
        qVal   : list of Qx-, Qy- and Qz-values
        dVec   : np.array of difference step in Qx, Qy and Qz
        qxGrid : 3D grid of Qx-values
        qyGrid : 3D grid of Qy-values
        qzGrid : 3D grid of Qz-values"""

        # aliases
        minBox = self.Qmin
        maxBox = self.Qmax
        dVec   = (maxBox - minBox) / self.dQN
        
        # vector data set of the center of each grid part (bin)
        qxVal = np.arange(minBox[0], maxBox[0] - dVec[0]/2, dVec[0]) + dVec[0]/2
        qyVal = np.arange(minBox[1], maxBox[1] - dVec[1]/2, dVec[1]) + dVec[1]/2
        qzVal = np.arange(minBox[2], maxBox[2] - dVec[2]/2, dVec[2]) + dVec[2]/2
        self.qVal = [qxVal, qyVal, qzVal]
        self.dVec = dVec
        
        # qx, qy, qz as gird in order z,y,x
        qx, qy, qz = get3DMesh(self.qVal[0], self.qVal[1], self.qVal[2])
        # qx, qy, qz as gird in order x,y,z like data grid and everything else
        self.qxGrid = qx.transpose(2,1,0)
        self.qyGrid = qy.transpose(2,1,0)
        self.qzGrid = qz.transpose(2,1,0)

    def _q2n(self, q):
        """Transform q-vector in grid indices

        q : q-vector as np.array

        return
        n : grid indices as np.array"""
        n = list((q - self.Qmin)/(self.Qmax - self.Qmin)*self.dQN).astype(int)
        if (n > np.array(self.dQN)).all():
            print '\n\nXXXX q-vector %s gives grid indices %s,' % (q, n)
            print '---- which are bigger than grid size %s' % (self.dQN)
            n = np.array(self.dQN)/2
            print '---- indices set to %s' % (n)
        return n

    def _box2mask(self, maskBox):
        """Transform box-values to mask grid region

        maskBox  : [xmin, ymin, zmin, xmax, ymax, zmax] as np.array

        return
        maskGrid : masking grid, True if not in considered region"""

        # qx, qy, qz as gird in order x,y,z                                                                     
        qx, qy, qz = self.qxGrid, self.qyGrid, self.qzGrid

        xMask = (qx >= maskBox[0]) & (qx <= maskBox[3])
        yMask = (qy >= maskBox[1]) & (qy <= maskBox[4])
        zMask = (qz >= maskBox[2]) & (qz <= maskBox[5])

        maskGrid = xMask & yMask & zMask
        return maskGrid

    def _calcMax(self):
        """Calculates the position of the maximum as indicies

        stores
        maxInd : indices of intensity maximum as np.array"""
        maxN = self.gridData.argmax()
        ind2 = maxN % self.dQN[2]
        ind1 = maxN / self.dQN[2] % self.dQN[1]
        ind0 = maxN / self.dQN[2] / self.dQN[1]
        self.maxInd = np.array([ind0, ind1, ind2])

    def _calcCutInd(self):
        """Calculates the cutting indices regarding the cutting modes

        'fix' : from the given position
        'max' : from the position of maximum"""
        
        if self.cutMode == 'fix':
            self.cutInd = self._q2n(self.cutPos)
        elif self.cutMode == 'max':
            self.cutInd = self.maxInd

        # for info file
        self.opProcInfo += self._makeCutIndInfo()

    def _processBgnd(self, maskBox = None):
        """Background subtraction of grid

        maskBox  : masked box [xMin, yMin, zMin, xMax, yMax, zMax]

        stores
        pBack    : parameters of background [mx, my, mz, d]
        gridBack : grid of the background
        maskBack : True if point in maskBox (False if part of background)
        maskFit  : False if point used for fit"""

        # function to discribe background
        bgndfunc = self._threeDLinRes

        # qx, qy, qz as gird in order x,y,z
        qx, qy, qz = self.qxGrid, self.qyGrid, self.qzGrid
                        
        # background and fit mask
        if maskBox == None:
            maskBox = self.backMaskBox
        if maskBox == None:
            self.maskBack = np.ones(self.gridData.shape) == 0
            self.maskFit  = self.maskOccu
        else:
            self.maskBack = self._box2mask(maskBox)
            self.maskFit  = self.maskOccu | self.maskBack

        # considered positions for background fit
        conMask = (self.maskFit == False)
        _qx = np.ravel(qx[conMask])
        _qy = np.ravel(qy[conMask])
        _qz = np.ravel(qz[conMask])

        # background subtraction
        if self.backType == 'mean':
            backMean      = self.gridData[conMask].mean()
            self.pBack    = [backMean]
            self.gridBack = np.ones(qx.shape) * backMean
        elif self.backType == 'threeDLin':
            bgndFunc = self._threeDLin
            bgndRes  = self._threeDLinRes
            guess    = [1e-6, 1e-6, 1e-6, self.gridData.mean()]
            fMes     = np.ravel(self.gridData[conMask])
            plsq     = leastsq(bgndRes, guess, args = (fMes, _qx, _qy, _qz))
            self.pBack    = plsq[0]
            self.gridBack = bgndFunc(plsq[0], qx, qy, qz)
        self.gridBack[self.maskOccu] = 0.0
        self.gridData = self.gridData - self.gridBack     

        # for info file
        self.opProcInfo += self._makeGridBackInfo(maskBox)
        
    def _threeDLin(self, p, x, y, z):
        """3D linear function

        p       : parameter [mx, my, mz, d]
        x, y, z : coordinates"""
        
        mx, my, mz, d = p
        f = (mx * x) + (my * y) + (mz * z) + d
        return f

    def _threeDLinRes(self, p, f, x, y, z):
        """3D linear function residual

        p       : parameter [mx, my, mz, d]
        f       : measured value
        x, y, z : coordinates"""
        
        err = f - self._threeDLin(p, x, y, z)
        return err


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

    def _calcIntMask(self, maskBox):
        """Calculate the mask for the integration region

        maskBox : [xmin, ymin, zmin, xmax, ymax, zmax] as np.array

        stores
        intMask : masking grid, True if not in integration region"""

        nMin = self._q2n(maskBox[:3])
        nMax = self._q2n(maskBox[3:])

    #
    # help functions for input / output
    #

    def _makeSetInfo(self):
        """Create the information about the set of images"""
        
        self.opSetInfo = '**** %s%s' % (self.setName, self.setNum)
        try:
            self.opSetInfo += '\nImages: %s' % (self.procImSelect)
        except:
            pass
        try:
            self.opSetInfo += '\nPhoton Wavelength: \t %.2f Angst.' % (self.waveLen)
        except:
            pass
        try:
            self.opSetInfo += '\nPhoton Energy: \t\t %.2f eV' % (self.energy)
        except:
            pass
    
    def _makeModeInfo(self, mode = None):
        """Create info description of the used frame mode

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        self.opModeInfo = '\n\n**** '

        if mode == None:
            mode = self.frameMode

        if mode == 1:
            self.opModeInfo += 'Frame Mode 1: (Qx, Qy, Qz) in theta-frame and (1/Angstrom)' 
        if mode == 2:
            self.opModeInfo += 'Frame Mode 2: (Qx, Qy, Qz) in phi-frame and (1/Angstrom)'
        if mode == 3:
            self.opModeInfo += 'Frame Mode 3: (Qx, Qy, Qz) in cartesian-frame and (1/Angstrom)'
        if mode == 4:
            self.opModeInfo += 'Frame Mode 4: (H, K, L) in hkl-frame and (reciprocal lattice units)'
        

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

        fitInfo  = '\n\n**** %s%s' % (fitTitle, fitType)
        fitInfo += '\n\t a \t\t b \t\t cen \t\t width \t\t width/step \t area' 
        line    = '\n%s \t ' + 4*'%.5e\t ' + '%.2f\t\t ' + '%.5e'
        for i in range(allRes.shape[0]):
            fitInfo += line % (fitNames[i], allRes[i,0], allRes[i,1], allRes[i,2], allRes[i,3], allRes[i,3]/self.dVec[i], allRes[i,4])

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
            dg = (gMax - gMin) / self.dQN
        emptNb = (gOccu == 0).sum()
        
        gridInfo = '\n\n**** %s sets processed to grid\n'   % (self.setEntLabel) + \
            'No. of bins in the grid : \t %.2e\n'      % (gData.size) + \
            'Data points outside the grid : \t %.2e\n' % (gOut) + \
            'Bins with zero information : \t %.2e'     % (emptNb)
        gridInfo += '\n\t min \t\t max \t\t step \t\t width'
        line = '\n%s' + 4*' \t %.2e'
        for i in range(gMin.size):
            gridInfo += line % (self.qLabel[i], gMin[i], gMax[i], dg[i], gMax[i]-gMin[i])

        return gridInfo

    def _makeGridBackInfo(self, maskBox):
        """Create information about background subtraction of the grid"""

        intInten, backInten = self.getIntIntensity()
    
        gridBackInfo  = '\n\n**** Background of grid fited linear'
        gridBackInfo += '\nTotal Intensity : \t %.4e' % (intInten + backInten)
        gridBackInfo += '\nGrid  Intensity : \t %.4e' % (intInten)
        gridBackInfo += '\nBackground Intensity : \t %.4e' % (backInten) 
        if maskBox == None:
            gridBackInfo += '\nNo masked region'
        else:
            gridBackInfo += '\nMasked region:'
            gridBackInfo += '\n\t min \t\t max'
            line = '\n%s' + 2*' \t %.2e'
            for i in range(3):
                gridBackInfo += line % (self.qLabel[i], self.backMaskBox[i], self.backMaskBox[i+3])

        gridBackInfo += '\nFitted parameters:'
        if self.backType == 'mean':
            gridBackInfo += '\nmean'
            gridBackInfo += ('\n%.2e') % self.pBack[0]
        elif self.backType == 'threeDLin':
            gridBackInfo += '\nmx \t\t my \t\t mz \t\t d'
            gridBackInfo += ('\n%.2e' + 3*' \t %.2e') % tuple(self.pBack)

            gridBackInfo += '\nBackground values at corners'
            cornList = [['min', 'min', 'min'], ['max', 'min', 'min'], ['max', 'max', 'min'],
                        ['min', 'max', 'min'], ['min', 'min', 'max'], ['max', 'min', 'max'], 
                        ['max', 'max', 'max'], ['min', 'max', 'max']]
            cornLayo = ['\n(%s_%s, ', '%s_%s, ', '%s_%s) : \t ']
            for i in range(len(cornList)):
                cornVal = self.pBack[3]
                for j in range(3):
                    gridBackInfo += cornLayo[j] % (self.qLabel[j], cornList[i][j])
                    cornVal      += getattr(self, 'Q' + cornList[i][j])[j]*self.pBack[j]
                    gridBackInfo += '%.2e' % (cornVal)
        
        return gridBackInfo

    def _makeCutIndInfo(self):
        """Create information about the cutting position of the grid"""

        cutIndInfo = '\n\n**** Cutting the grids'
        if self.cutMode == 'fix':
            cutIndInfo += ' at fixed position'
        elif self.cutMode == 'max':
            cutIndInfo += ' at maximum postion'
        cutIndInfo += '\n\t index \t value'
        line = '\n%s' + ' \t %d' + ' \t %.2e'
        for i in range(3):
            cutIndInfo += line % (self.qLabel[i], self.cutInd[i], self.qVal[i][self.cutInd[i]])

        return cutIndInfo

    #
    # process part
    #
    
    def processOneSet(self, procSelect = None, mode = None):
        """Process selcted images of the full set into (Qx, Qy, Qz, I)

        procSelect : list with the images which will be processed, take object default if None
        mode       : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        
        if self.processMode == 'builtin':

            print '\n%s%d: process to %s' % (self.setName, self.setNum, self.setEntLabel)

            if procSelect == None:
                procSelect = self.procImSelect
            if mode == None:
                mode = self.frameMode
            
            imSize = self.fileProcessor.getImage()[0].size
            npts = imSize * len(procSelect)

            if self.totSet is None:
                self.totSet = np.zeros((npts, 4))

            j = 0
            k = imSize
            for i in procSelect:
                print "**** Processing image %d" % i
                self._processOneImage(self.totSet[j:k,:], i, mode = mode)
                j = j + imSize
                k = k + imSize
        
        else:
            ccdToQkwArgs = {}
            if self.totSet is not None:
                del self.totSet
                gc.collect()
                #ccdToQkwArgs['outarray'] = self.totSet
                
            print "\n**** Converting to Q"
            t1 = time.time()
            self.totSet = ctrans.ccdToQ(mode        = self.frameMode,
                                        angles      = self.settingAngles * np.pi / 180.0, 
                                        ccd_size    = (self.detSizeX, self.detSizeY),
                                        ccd_pixsize = (self.detPixSizeX, self.detPixSizeY),
                                        ccd_cen     = (self.detX0, self.detY0),
                                        dist        = self.detDis,
                                        wavelength  = self.waveLen,
                                        UBinv       = np.matrix(self.UBmat).I,
                                        **ccdToQkwArgs)
            t2 = time.time()
            print "---- DONE (Processed in %f seconds)" % (t2 - t1)
            self.totSet[:,3] = np.ravel(self.fileProcessor.getImage())
            #print gc.get_referents(totSet)

        # for info file
        self.opProcInfo += '\n\n**** Image Set processed to %.2e %s sets' % (self.totSet.shape[0], self.setEntLabel)

    def makeGridData(self, procSelect = None, mode = None, delete = False, backSub = None):
        """Grid the data set into a cuboid
        Size of the cuboid : Qmin, Qmax = [Qx, Qy, Qz]_min, max
        Number of parts    : dQN = [Nqx, Nqy, Nqz]

        procSelect : list with the images which will be processed, take object default if None
        mode       : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None
        backSub    : subtract linear background if True; if None take class default

        stores
        gridData : values at the grid parts (bins)
        gridOccu : occupation no. of the grid parts (bins)
        gridBack : background values if backSub == True
        gridOut  : no. of the data points outside the grid
        maskData : mask True if gridData == 0
        maskOccu : mask True if gridOccu == 0
        maskBack : mask True if point not for background fit, only if backSub == True"""

        if procSelect == None:
            procSelect = self.procImSelect
        if mode == None:
            mode = self.frameMode

        self.processOneSet(procSelect = procSelect, mode = mode)
        print "---- Total data is %f MBytes\n" % (self.totSet.nbytes / 1024.0**2)
        
        # prepare min, max,...
        if (self.Qmin == np.array(None)).all():
            self.Qmin = np.array([ self.totSet[:,0].min(), self.totSet[:,1].min(), self.totSet[:,2].min() ])
        if (self.Qmax == np.array(None)).all():
            self.Qmax = np.array([ self.totSet[:,0].max(), self.totSet[:,1].max(), self.totSet[:,2].max() ])
        if (self.dQN  == np.array(None)).all():
            self.dQN = [100, 100, 100]

        # 3D grid of the data set 
        print "**** Gridding Data."
        t1 = time.time()
        gridData, gridOccu, gridOut = ctrans.grid3d(self.totSet, self.Qmin, self.Qmax, self.dQN, norm = 1)
        t2 = time.time()
        print "---- DONE (Processed in %f seconds)" % (t2 - t1)
        emptNb = (gridOccu == 0).sum()
        if gridOut != 0:
            print "---- Warning : There are %.2e points outside the grid (%.2e bins in the grid)" % (gridOut, gridData.size)
        if emptNb:
            print "---- Warning : There are %.2e values zero in the grid" % emptNb

        # mask the gridded data set
        #gridData = np.ma.array(gridData / gridOccu, mask = (gridOccu == 0))
        # store intensity, occupation and no. of outside data points of the grid
        self.gridData = gridData
        self.gridOccu = gridOccu
        self.gridOut  = gridOut

        # masks for the data and occupation no.
        self.maskData = (gridData == 0)
        self.maskOccu = (gridOccu == 0)

        #del self.totSet
        #gc.collect()

        # calculated the corresponding vectors and maximum intensity position of the grid
        self._calcVecDataSet()

        # for info file
        self.opProcInfo += self._makeGridInfo()

        # background subtraction
        if backSub == None:
            backSub = self.backSub
        if backSub == True:
            self._processBgnd()

        # position of maximum and cutting of the data grid
        self._calcMax()
        self._calcCutInd()

    def get1DSum(self, selType = 'gridData'):
        """1D Lines of the selected grid by summing in the other directions

        selType : select type of summed grid, e.g. 'gridData'

        returns
        oneDSum : set in the order Qx, Qy, Qz as list"""

        try:
            selGrid = getattr(self, selType)
            oneDSum = [selGrid.sum(1).sum(1),
                       selGrid.sum(0).sum(1),
                       selGrid.sum(0).sum(0)]
        except:
            print 'xxxx %s is not a corecct type of a grid from ImageProcessor!' % (selType)
        
        return one1DSum

    def get2DSum(self, selType = 'gridData'):
        """2D Areas of the selected grid by summing in the other direction

        selType : select type of summed grid, e.g. 'gridData'

        returns
        twoDSum : set in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        try:
            selGrid = getattr(self, selType)
            twoDSum = [selGrid.sum(0), selGrid.sum(1), selGrid.sum(2)]
        except:
            print 'xxxx %s is not a corecct type of a grid from ImageProcessor!' % (selType)

        return twoDSum
    
    def get1DCut(self, selType = 'gridData', cutInd = None):
        """1D Lines of the selected grid at the cut position

        selType : select type of cutted grid, e.g. 'gridData'
        cutInd  : cut indices, [nx, ny, nz], if None, default
                  if default None, indicies of maximum is taken

        returns
        oneDCut : set in the order Qx, Qy, Qz as list"""
        
        if cutInd == None:
            cutInd = self.cutInd
        if cutInd == None:
            cutInd = self.maxInd

        try:
            selGrid = getattr(self, selType)
            oneDCut = [selGrid[:, cutInd[1], cutInd[2]],
                       selGrid[cutInd[0], :, cutInd[2]],
                       selGrid[cutInd[0], cutInd[1], :]]
        except:
            print '\nxxxx %s is not a corecct type of a grid from ImageProcessor!\n' % (selType)
           
        return oneDCut
    
    def get2DCut(self, selType = 'gridData', cutInd = None):
        """2D Areas of the selected grid at the cut position

        selType : select type of cutted grid, e.g. 'gridData'
        cutInd  : cut indices, [nx, ny, nz], if None, default
                  if default None, indicies of maximum is taken

        returns
        twoDCut : set in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        if cutInd == None:
            cutInd = self.cutInd
        if cutInd == None:
            cutInd = self.maxInd 

        try:
            selGrid = getattr(self, selType)
            twoDCut = [selGrid[cutInd[0], :, :],
                       selGrid[:, cutInd[1], :],
                       selGrid[:, :, cutInd[2]]]
        except:
            print '\nxxxx %s is not a corecct type of a grid from ImageProcessor!\n' % (selType)
        
        return twoDCut

    def get1DCutAv(self):
        """1D averaged Lines of the grid data and occupations at the position of the maximum
        intensity and its eight neighbored lines 

        cutInd : cut indices, [nx, ny, nz], if None, default
                 if default None, indicies of maximum is taken

        returns
        gridData1DCutAv : intensity set  in the order Qx, Qy, Qz as list
        gridOccu1DCutAv : occupation no. in the order Qx, Qy, Qz as list"""

        if cutInd == None:
            cutInd = self.cutInd
        if cutInd == None:
            cutInd = self.maxInd 

        # initialize with correct size as zeros
        gridData1DCutAv = [np.zeros(self.dQN[0]),np.zeros(self.dQN[1]),np.zeros(self.dQN[2])]
        gridOccu1DCutAv = [np.zeros(self.dQN[0]),np.zeros(self.dQN[1]),np.zeros(self.dQN[2])]

        #print self.dQN[0]
        # go through the neighbors
        for i in range(3):
            for j in range(3):
                gridData1DCutAv[0] += self.gridData[:,self.cutInd[1]+i-1,self.cutInd[2]+j-1]/9.0
                gridData1DCutAv[1] += self.gridData[self.cutInd[0]+i-1,:,self.cutInd[2]+j-1]/9.0
                gridData1DCutAv[2] += self.gridData[self.cutInd[0]+i-1,self.cutInd[1]+j-1,:]/9.0
                gridOccu1DCutAv[0] += self.gridOccu[:,self.cutInd[1]+i-1,self.cutInd[2]+j-1]/9.0
                gridOccu1DCutAv[1] += self.gridOccu[self.cutInd[0]+i-1,:,self.cutInd[2]+j-1]/9.0
                gridOccu1DCutAv[2] += self.gridOccu[self.cutInd[0]+i-1,self.cutInd[1]+j-1,:]/9.0
        
        return gridData1DCutAv, gridOccu1DCutAv

    def get2DCutAv(self):
        """2D average Areas of the grid data and occupations at the position of the maximum
        intensity and the their two neighbors

        cutInd : cut indices, [nx, ny, nz], if None, default
                 if default None, indicies of maximum is taken

        returns
        gridData2DCutAv : intensity set  in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list
        gridOccu2DCutAv : occupation no. in the order (Qy, Qz), (Qx, Qz), (Qx, Qy) as list"""

        if cutInd == None:
            cutInd = self.cutInd
        if cutInd == None:
            cutInd = self.maxInd 

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
            gridData2DCutAv[0] += self.gridData[self.cutInd[0]+i-1,:,:]/3.0
            gridData2DCutAv[1] += self.gridData[:,self.cutInd[1]+i-1,:]/3.0
            gridData2DCutAv[2] += self.gridData[:,:,self.cutInd[2]+i-1]/3.0
            gridOccu2DCutAv[0] += self.gridOccu[self.cutInd[0]+i-1,:,:]/3.0
            gridOccu2DCutAv[1] += self.gridOccu[:,self.cutInd[1]+i-1,:]/3.0
            gridOccu2DCutAv[2] += self.gridOccu[:,:,self.cutInd[2]+i-1]/3.0
        return gridData2DCutAv, gridOccu2DCutAv

    def getIntIntensity(self, maskBox = None):
        """Get the integrated intensity of the peak by summing over all grid parts (bins)

        maskBox   : [xmin, ymin, zmin, xmax, ymax, zmax] as np.array

        returns
        intInten  : integrated intensity
        backInten : background intensity"""

        # prepare considered region with mask
        if maskBox == None:
            try:
                maskBox = self.intMaskBox
            except:
                pass
        if maskBox == None:
            conMask = np.ones(self.gridData.shape).astype(int)
        else:
            conMask = (self._box2mask(maskBox) == False)

        # integrated intensity and background in considered region
        intInten = self.gridData[conMask].sum()
        if self.backSub == True:
            backInten = self.gridBack[conMask].sum()
        else:
            backInten = 0.0
        
        return intInten, backInten

    #
    # fit part
    #

    def get1DFit(self, xVal, yVal, fitType = None, infoDes = ''):
        """Fit a 1D data set

        xVal    : x-values as list of arrays
        yVal    : y-values as list of arrays
        fitType : peak shape to fit from pyspec fitfuncs, use object default if None, e.g. 'lor2a'
        infoDes : description of the current fit try for output, e.g. 'Line cut of Scan #244'

        returns
        yFit    : y-values of the fit
        fitRes  : results of the fit, [a1, b1, cen1, width1, area1], [0, 0, 0, 0, 0] if unsuccessful fit"""

        if fitType == None:
            fitType = self.fitType

        yFit   = []
        for i in range(len(xVal)):
            yFit.append(np.zeros(len(xVal[i])))
        fitRes = np.zeros((len(xVal),5))
        for i in range(len(xVal)):
            try:
                f = fit.fit(x=xVal[i], y=yVal[i], funcs = [fitfuncs.linear, getattr(fitfuncs, fitType)])
                f.go()
                yFit[i]   = fitfuncs.linear(xVal[i], f.result[:2]) + getattr(fitfuncs, fitType)(xVal[i], f.result[2:])
                fitRes[i] = f.result
            except:
                print 'WARNING : %s %s could not be fitted to %s!' % (self.qLabel[i], infoDes, fitType)

        # for info file
        self.opProcInfo += self._makeFitInfo1D(fitRes, fitType = None, fitTitle = infoDes, fitNames = self.qLabel)

        return yFit, fitRes
    
    #
    # input / output part
    #

    def makeInfo(self):
        """Create the information about the current processing"""

        self._makeSetInfo()
        curInfo = '%s%s%s' % (self.opSetInfo, self.opModeInfo, self.opProcInfo)

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

    from pyspec import spec
    #sf   = spec.SpecDataFile('/home/tardis/swilkins/data/lbco5/lbco125.01', 
    #			     ccdpath = '/mounts/davros/nasshare/images/nov10')
    #scan = sf[91]

    sf   = spec.SpecDataFile('/home/tardis/spartzsch/data/ymn2o5_oct10/ymn2o5_oct10_1', 
    			     ccdpath = '/mounts/davros/nasshare/images/oct10')
    scan = sf[360]

    ###
    # file processor
    ###

    fp = FileProcessor()
    #fp.setFromSpec(scan)
    #fp.process()
    #fp.processBgnd(maskroi =[[100, 100, 100, 100]])
    #fp.save('/mounts/timelord/storage/tmpimage.npz')
    #fp.load('/mounts/timelord/storage/tmpimage.npz')

    ###
    # image processor and options
    ###

    testData = ImageProcessor(fp)
    #testData.processMode = 'builtin'
    testData.setDetectorPos(detAng = -1.24)
    testData.setBins(4, 4)
    #testData.setFileProcessor(fp)
    testData.setSpecScan(scan)
    #fp.processBgnd(maskroi =[[63, 100, 67, 100]])
    #testData.setConRoi([1, 325, 1, 335])
    testData.setFrameMode(4)
    testData.setGridOptions(Qmin = None, Qmax = None, dQN = [90, 80, 80])
    # ICM, HKL
    #testData.setCutPos([0.4750, 0.0, 0.2856])
    # CM, HKL
    #testData.setCutPos([0.4924, 0.0, 0.246])
    #testData.processMode = 'builtin'
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [200, 400, 100])
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [100, 100, 100])

    ###
    # grid of image set
    ###    
    
    #testData.makeGridData(procSelect = [40])
    testData.setGridBackOptions(backSub = True, backMaskBox = [0.475, -0.009, 0.235, 0.508, 0.009, 0.255], backType = 'mean')
    testData.setIntegrationOptions(intMaskBox = [0.475, -0.009, 0.235, 0.508, 0.009, 0.255])
    testData.makeGridData()

    #testData._processBgnd(m))


    #print 'Peak integrated intensity : %.2e' % (testData.getIntIntensity())
    #lineData, lineOccu = testData.get1DCut()
    #areaData, areaOccu = testData.get2DCut()

    ###
    # plot of raw images
    ###

    #plIm = PlotImages(fp, testData)
    #plIm.setPlotContent(plotSelect = range(0,121,10), plotType = 'norm')
    #plIm.setPlotLayout(, plotImHor = 4, plotImVer = 3)
    #plIm.plotImages()

    ###
    # plot grid data
    ###

    testPlotter = PlotGrid(testData)
    #testPlotter.setPlotFlags(7, 7)
    testPlotter.setLogFlags(7, 7)
    #testPlotter.setFit1D(True)
    #testPlotter.setHistBin(50)
    testPlotter.setPlot1DFit(True)

    # plot jobs

    #testPlotter.plotGrid1D('sum')
    #testPlotter.plotGrid2D('sum')
    testPlotter.plotGrid1D('cut')
    testPlotter.plotGrid2D('cut')
    #testPlotter.plotMask1D('cut')
    testPlotter.plotMask2D('cut')

    ###
    # processing information
    ###

    info1 = testData.makeInfo()
    #testData.writeInfoFile(outFile = 'infoFile.dat')

    ###
    # second run with same objects
    ###

    scan2 = sf[361]
    #testData.setSpecScan(scan2)
    #testData.makeGridData()
    #testPlotter.plotGrid1D('cut')
    #testPlotter.plotGrid2D('cut')

    #info2 = testData.makeInfo()

    print '\n\n'
    print info1
    #print info2

    plt.show()

    #raw_input('Waiting...')

    #del testData
    #del fp
    #gc.collect()
    #raw_input('Waiting...')
     

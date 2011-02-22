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
    pass

from   pyspec.ccd.plotter import PlotImages, PlotGrid

try:
    import pyspec.ccd.ctrans as ctrans
except:
    try:
        import ctrans
    except:
        pass

from ConfigParser import ConfigParser

#gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_LEAK)
gc.enable()

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>, Sven Partzsch <SvenPartzsch@gmx.de>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

#
# FileProcessor Class
#

class FileProcessor():
    """FileProcessor Class

    This class processes CCD files and returns a numpy array of 
    float values for all the images. Images can be saved to disk
    for faster processing in the future.

    The basic building block of this concept is a list of all filenames
    to process. This should be provided through the helper functions or
    can be set through a specfile object.

    """

    def __init__(self, filenames = None, 
                 darkfilenames = None, 
                 norm = None, format = 'SPE',
                 spec = None):
        """Initialize the class

        filenames      : list of filenames for all images (2d list will bin images)
        darkfilenames  : list of filenames for the dark images (2d list will bin images)
        spec           : SpecScan object from which to obtain data file names
        norm           : Normalize individual images to this numpy array.
        format         : Data file format. 'SPE' princeton format.
        """
        self._format = format
        self.filenames = filenames
        self.darkfilenames = darkfilenames
        self.normData = None
        self.meanMonitor = False

        if spec is not None:
            self.setFromSpec(spec)
        if norm is not None:
            self.normData = norm
        self.bgndParams = np.array([])

    def setFilenames(self, filenames = None, darkfilenames = None):
        """Set the list of filenames and darkfilenames"""
        self.filenames = filenames
        self.darkfilenames = dar

    def setMeanMonitor(self, b):
        """Set if the images are normalized by the mean of the monitor

        If True then normalize all images by the mean of the monitor counts present"""
        
        self.meanMonitor = b

    def setFromSpec(self, scan, mon = 'Monitor'):
        """Set the filenames from a SpecScan instance

        scan    : SpecScan instance
        norm    : spec counter to normalize against"""
        
        self.filenames     = scan.getCCDFilenames()
        # alternative because dark image not any more provided by pyspec Scan
        self.darkfilenames = len(self.filenames) * [[self.filenames[0][0][:-9] + '-DARK_0000.spe']]
        self.normData      = scan.values[mon]

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

        maskroi : list of tuples
           list of 4 element tuples or lists for ROI to mask
        mask    : ndarray
           (n x m) array for the background mask"""

        print "---- Subtracting background from images"
        self._processBgnd(maskroi = maskroi, mask = mask)
        print "---- Done."

    def process(self, dark = True, norm = True, dtype = np.float,
                quiet = False):
        """Read in images and process into 3D array.
        
        dark  : bool
           If True, subtract the dark images from the data
        norm  : bool
           If True, normalize by monitor.
        dtype : datatype 
           numpy datatype of processed array.
        quiet : bool
           If True, dont write to screen when reading images."""

        images = []
        darkimages = []

        if norm:
            if self.normData is None:
                normData = np.ones(len(self.filenames))
                print "XXXX No normalization data found"
            else:
                if self.meanMonitor:
                    normData = np.ones(len(self.filenames)) * self.normData.mean()
                    print "---- Normalizing data (using mean value)."
                else:
                    normData = self.normData
                    print "---- Normalizing data."
        else:
            normData = np.ones(len(self.filenames))

        if dark:
            print "---- Correcting for dark current"

        if quiet:
            print "---- Reading Images"

        for i, (iname, diname, normVal) in enumerate(zip(self.filenames, self.darkfilenames, normData)):
            if type(iname) == list:
                _images = None
                _darkimages = None
                
                #Start reading the light images
                for j, _in in enumerate(iname):
                    if os.path.exists(_in):
                        image = self._getRawImage(_in).astype(dtype)
                        if _images is not None:
                            _images = _images + image
                        else:
                            _images = image
                        #_images.append(image)
                        if not quiet:
                            print "---- Reading image %-3d of %-3d (sub image %-3d of %-3d)     \r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()
                    else:
                        if not quiet:
                            print "---- Missing image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()

                for j, _din in enumerate(diname):
                    if os.path.exists(_din):
                        darkimage =  self._getRawImage(_din).astype(dtype)
                        if _darkimages is not None:
                            _darkimages = _darkimages + darkimage
                        else:
                            _darkimages = darkimage
                        #_darkimages.append(darkimage)
                        if not quiet:
                            print "---- Reading dark image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.darkfilenames), j + 1, len(diname)),
                            sys.stdout.flush()
                    else:
                        if not quiet:
                            print "---- Missing dark image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.darkfilenames), j + 1, len(diname)),
                            sys.stdout.flush()

                image = _images
                if _darkimages is not None:
                    darkimage = _darkimages
                    darkimages.append(darkimage)
                #image = np.array(_images).sum(0)
                #if len(_darkimages):
                #    darkimage = np.array(_darkimages).sum(0)
                #    darkimages.append(darkimage)

            else:
                # Process only single image pair
                image = self._getRawImage(iname).astype(dtype)
                if os.path.exists(diname):
                    darkimage =  self._getRawImage(diname).astype(dtype)
                    darkimages.append(darkimage)
                if not quiet:
                    print "---- Reading image %-3d of %-3d\r" % (i, len(self.filenames)),

            if dark:
                if len(darkimages):
                    image = image - darkimages[-1]
                else:
                    print "XXXX Unable to dark currect correct. No Image found"
            
            if norm:
                image = image / normVal
            images.append(image)

        print "\n---- Processed %d images (%d dark images)" % (len(images), len(darkimages))

        self.images = np.array(images)

        print "---- Done. Array size %s (%.3f Mb)." % (str(self.images.shape), 
                                                       self.images.nbytes / 1024**2) 

    def _getRawImage(self, iname):
        """Read raw image"""
        if self._format == 'SPE':
            return PrincetonSPEFile(iname).getBinnedData()
        else:
            raise Exception("Unknown file format \"%s.\"" % self._format)

    def _computeMeanImage(self):
        N = self.images.shape[0]
        self.mean = self.images.sum(0) / N
        stdev = (self.images - self.mean)**2

    def getMeanImage(self):
        """Return the mean image after processing.

        Calculates the mean of all images processed. Returns a tuple
        of the mean image and the standard error"""
        _computeMeanImage()
        return self.mean, self.stderr

    def getImage(self, n = None):
        """Return the image data

        n : int or None
           Image number to return. If None, return all images."""

        if n is None:
            return self.images
        else:
            return self.images[n]

    def saveImage(self, filename, inum = None, dtype = np.float32):
        """Save image to binary file of specified datatype"""
        if inum is None:
            self.images.astype(dtype).tofile(filename)
        else:
            self.images[n].astype(dtype).tofile(filename)
                               

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
        obj = np.load(filename)
        self.images = obj['images']
        self.normData = obj['images']

#
# image processor class
#

class ImageProcessor():
    """Image Processor class

    This class provides the processing of single and sets of CCD-images
    in more detail: each pixel is transformed into reciprocal space
    the set of reciprocal vectors and intensities is gridded on a regular cuboid
    the needed informations can be provided by a spec scan from the pyspec package"""

    def __init__(self, fP = None, configfile = None, ccdname = 'CCD',
                 spec = None):
        """Initialize the image processor

        fP         : file processor object for getting the images
        configfile : file defining config of CCD
        ccdname    : name of CCD (used in config file)
        spec       : SpecScan object used to obtain information"""
        
        # file processor to provied the needed images
        self.fileProcessor = fP

        self.ccdName = ccdname
        
        # set parameters to configure the CCD setup
        # detector distance 30cm and detector pixel size 20um
        self.setDetectorProp(0.020, 0.020, 1300, 1340, 650.0, 670.0)
        self.setDetectorPos(300, 0.0)
        # Overide if a config file is used.
        if configfile:
            self.readConfigFile(configfile)
        
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

    def getImageData(self):
        """Get the totat image data, transformed into Q"""
        return self.totSet

    def getGridData(self):
        """Get the data on the grid"""
        return self.gridData

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

    def setSetSettings(self, waveLen, imFilesNames, darkFileNames,
                       settingAngles, intentNorm, UBmat, setName,
                       setNum, setSize):
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
        # Sven : Don't do this. Always wait for the user to iniciate this!
        #self.fileProcessor.process() # NO NO NO NO NO NO !!!!!
   
    def getSpecScan(self):
        """Get the pyspec scan object which was used for the set settings
        
        returns pyspec scan object which contains the needed information"""

        return self.conScan

    def setConRoi(self, conRoi):
        """Set the region of interest used to process the images

        conRoi : [xMin, xStep, yMin, yStep]"""
        
        self.conRoi = conRoi

    def getConRoi(self):
        """Get the region of interest used to process all the images

        conRoi : [xMin, xStep, yMin, yStep]"""
        
        return self.conRoi

    def setFrameMode(self, mode):
        """Set the mode of the output frame for (Qx, Qy, Qz)

        The image processor uses a number of modes which defile the
        coordinates of the final grid. These are:

        mode 1 : 'theta'    : Theta axis frame.  
        mode 2 : 'phi'      : Phi axis frame.
        mode 3 : 'cart'     : Crystal cartesian frame.
        mode 4 : 'hkl'      : Reciproal lattice units frame."""

        if mode == 'theta':
            self.frameMode = 1
        elif mode == 'phi':
            self.frameMode = 2
        elif mode == 'cart':
            self.frameMode = 3
        elif mode == 'hkl':
            self.frameMode = 4
        else:
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

        This function returns the X, Y and Z coordinates of the grid as 3d
        arrays.

        Example:

        X, Y, Z = ip.getGridMesh()

        These values can be used for obtaining the coordinates of each voxel.
        For instance, the position of the (0,0,0) voxel is given by

        x = X[0,0,0]
        y = Y[0,0,0]
        z = Z[0,0,0]"""
        
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
   
    #
    # get set functions for input output
    #

    def setFileProcessor(self, fp = None):
        self.fileProcessor = fp

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
        # alternative because dark image not any more provided by pyspec Scan
        self.darkFileNames = len(self.imFileNames) * [[self.imFileNames[0][0][:-9] + '-DARK_0000.spe']]
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
        n = np.array((q - self.Qmin)/(self.Qmax - self.Qmin)*self.dQN).astype(int)
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
        maskGrid : masking grid, True if in maskBox"""

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
    
    #
    # help functions for input / output
    #

    def _makeSetInfo(self):
        """Create the information about the set of images"""
        
        self.opSetInfo = '**** %s%s\n' % (self.setName, self.setNum)
        try:
            self.opSetInfo += '\nImages: %s\n' % (self.procImSelect)
        except:
            pass
        try:
            self.opSetInfo += '\nPhoton Wavelength: \t %.2f Angst.\n' % (self.waveLen)
        except:
            pass
        try:
            self.opSetInfo += '\nPhoton Energy: \t\t %.2f eV\n' % (self.energy)
        except:
            pass
    
    def _makeModeInfo(self, mode = None):
        """Create info description of the used frame mode

        mode : 1 (theta-) , 2 (phi-), 3 (cartesian-) or 4 (hkl-frame), take object default if None"""

        self.opModeInfo = '\n\n**** '

        if mode == None:
            mode = self.frameMode

        if mode == 1:
            self.opModeInfo += 'Frame Mode 1: (Qx, Qy, Qz) in theta-frame and (1/Angstrom)\n' 
        if mode == 2:
            self.opModeInfo += 'Frame Mode 2: (Qx, Qy, Qz) in phi-frame and (1/Angstrom)\n'
        if mode == 3:
            self.opModeInfo += 'Frame Mode 3: (Qx, Qy, Qz) in cartesian-frame and (1/Angstrom)\n'
        if mode == 4:
            self.opModeInfo += 'Frame Mode 4: (H, K, L) in hkl-frame and (reciprocal lattice units)\n'

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
    
        gridBackInfo  = '\n\n**** Background subtraction of grid'
        if self.backType == 'mean':
            gridBackInfo += ' by mean value'
        elif self.backType == 'threeDLin':
            gridBackInfo += ' by fitted 3D linear function'
        gridBackInfo += '\n---- Masked region'
        if maskBox == None:
            gridBackInfo += '\nNo'
        else:
            gridBackInfo += '\n\t min \t\t max'
            line = '\n%s' + 2*' \t %.2e'
            for i in range(3):
                gridBackInfo += line % (self.qLabel[i], maskBox[i], maskBox[i+3])

        gridBackInfo += '\n---- Fitted parameters:'
        if self.backType == 'mean':
            gridBackInfo += '\nmean'
            gridBackInfo += ('\n%.2e') % self.pBack[0]
        elif self.backType == 'threeDLin':
            gridBackInfo += '\nmx \t\t my \t\t mz \t\t d'
            gridBackInfo += ('\n%.2e' + 3*' \t %.2e') % tuple(self.pBack)

            gridBackInfo += '\n---- Background values at corners'
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

    def _makeIntIntensityInfo(self):
        """Create information about integrated intensity of the grid"""

        intInten, backInten = self.getIntIntensity()
    
        gridIntInfo  = '\n\n**** Integrated intensities of the grid'
        gridIntInfo += '\nTotal Intensity : \t %.4e' % (intInten + backInten)
        gridIntInfo += '\nGrid  Intensity : \t %.4e' % (intInten)
        gridIntInfo += '\nBackground Intensity : \t %.4e' % (backInten)
        gridIntInfo += '\n---- Considered region:'
        try:
            maskBox = self.intMaskBox
        except:
            maskBox = None
        if maskBox == None:
            gridIntInfo += '\nAll'
        else:
            gridIntInfo += '\n\t min \t\t max'
            line = '\n%s' + 2*' \t %.2e'
            for i in range(3):
                gridIntInfo += line % (self.qLabel[i], maskBox[i], maskBox[i+3])

        return gridIntInfo

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
        ccdToQkwArgs = {}
        if self.totSet is not None:
            del self.totSet
            gc.collect()

        # if images not yet processed, do it
        try:
            getattr(self.fileProcessor, 'images')
        except:
            self.fileProcessor.process()
                    
        print "\n**** Converting to Q"
        t1 = time.time()
        self.totSet = ctrans.ccdToQ(angles      = self.settingAngles * np.pi / 180.0, 
                                    mode        = self.frameMode,
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
            
        # for info file
        self.opProcInfo += '\n\n**** Image Set processed to %.2e %s sets (Took %f seconds)' % (self.totSet.shape[0], self.setEntLabel, (t2 - t1))

    def makeGridData(self, procSelect = None, mode = None, delete = False, backSub = None):
        """Grid the data set into a cuboid.
        
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
        gridRes = ctrans.grid3d(self.totSet, self.Qmin, self.Qmax, self.dQN, norm = 1)
        #print gridRes
        print len(gridRes)
        for i in range(2):
            print gridRes[i].shape
        print gridRes[2]
        #gridData, gridOccu, gridStdErr, gridOut = ctrans.grid3d(self.totSet, self.Qmin, self.Qmax, self.dQN, norm = 1)
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
        self.gridData   = gridData
        self.gridOccu   = gridOccu
        self.gridOut    = gridOut
        self.gridStdErr = gridStdErr

        # masks for the data and occupation no.
        self.maskData = (gridData == 0)
        self.maskOccu = (gridOccu == 0)

        # calculated the corresponding vectors and maximum intensity position of the grid
        self._calcVecDataSet()

        # for info file about gridding
        self.opProcInfo += self._makeGridInfo()

        # background subtraction
        if backSub == None:
            backSub = self.backSub
        if backSub == True:
            self._processBgnd()

        # for info file about integrated intensities
        self.opProcInfo += self._makeIntIntensityInfo()

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
        
        return oneDSum

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
            conMask = self._box2mask(maskBox)

        # integrated intensity and background in considered region
        intInten = self.gridData[conMask].sum()
        if self.backSub == True:
            backInten = self.gridBack[conMask].sum()
        else:
            backInten = 0.0
        
        return intInten, backInten

    #
    # input / output part
    #

    def makeInfo(self):
        """Create and return the information about the current processing"""
        return str(self)

    def __str__(self):
        """Output the text about the current processing"""
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

class CCDParamsConfigParser(ConfigParser):
    """Class to read config file which defines all CCD parameters"""
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
    testData.setGridBackOptions(backSub = True,
                                backMaskBox = [0.475, -0.009, 0.235, 0.508, 0.009, 0.255],
                                backType = 'threeDLin')
    testData.setIntegrationOptions(intMaskBox = [0.475, -0.009, 0.235, 0.508, 0.009, 0.255])
    testData.makeGridData()

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
    #testPlotter.plotGrid1D('cut')
    #testPlotter.plotGrid2D('cut')
    #testPlotter.plotMask1D('cut')
    #testPlotter.plotMask2D('cut')

    

    ###
    # second run with same objects
    ###

    scan2 = sf[361]
    #testData.setSpecScan(scan2)
    #testData.makeGridData()
    #testPlotter.plotGrid1D('cut')
    #testPlotter.plotGrid2D('cut')

    ###
    # processing information
    ###

    info1 = testData.makeInfo()
    #testData.writeInfoFile(outFile = 'infoFile.dat')
    print '\n\n'
    print info1



    #info2 = testData.makeInfo()

    
    #print info2
    print '\n'

    plt.show()

    #raw_input('Waiting...')

    #del testData
    #del fp
    #gc.collect()
    #raw_input('Waiting...')
     

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
import struct
import numpy as np
import numpy.ma as ma
import operator
from   scipy.optimize import leastsq
import matplotlib.pyplot as plt
from   pyspec import fit, fitfuncs
from   pyspec.ccd import utils as ccdutils
from   pyspec.diffractometer import Diffractometer
try:
    import h5py
except:
    pass

try:
    from PIL import Image as PILImage
except:
    pass

try:
    from   pyspec.ccd.files import *
except:
    pass

try:
    from   pyspec.ccd.plotter import PlotImages, PlotGrid, PlotGrid2, PlotWindow
except:
    pass

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
                 spec = None, process = False,
                 spoolpath = None):
        """Initialize the class

        filenames      : list of strings.
           List of filenames for all images (2d list will bin images).
        darkfilenames  : list of strings.
           List of filenames fror all dark images. (2d list will bin images).
        spec           : SpecScan object.
           SpecScan object from which to obtain image file names.
        norm           : ndarray.
           Normalize individual images to this numpy array.
        format         : string.
           Data file format. 'SPE' princeton format.
        process        : Bool.
           If True, then process the images once class is initialized.
        spoolpath      : String
           Path of directories to spool data to upon processing.
        """
        self._format = format
        self.filenames = filenames
        self.darkfilenames = darkfilenames
        self.normData = None
        self.meanMonitor = False

        self.images = None
        self.mask = None

        self._cropOnRead = None
        self._binOnRead = None
        
        self.spoolFilename = None
        self.spoolfd = None
        self.spoolfdin = None

        if spec is not None:
            if spoolpath is not None:
                self.setFromSpec(spec, spoolpath = spoolpath, spool = True)
            else:
                self.setFromSpec(spec)
            self.setFromSpec(spec)

        if norm is not None:
            self.normData = norm

        self.bgndParams = np.array([])

        self.processed = False
        self.currentFrame = 0

        if process:
            self.process()

    def setCropOnRead(self, roi = None):
        """Sets the portion of the image to crop on loading each image.

        This function sets the region of interest which will be taken from
        each image after reading from disk and before placing in the
        image processor object. This can be used to reduce the size of an image
        stack.

        roi : None or tuple
            The ROI. If None use whole image. The tuple should be formated
            as (xstart, ystart, ystop, ystop). Please note, this uses python
            standard indexing. The stop values are NOT included in the crop.

        """

        self._cropOnRead = roi

    def setCropOnRead(self,bins = None):
        """Sets the portion of the image to bin on loading each image.

        This function will bin the images before loading them into the image processor
        
        bins : None or tuple
               The number of bins in the x and y directions. If None use (1x1) binning.
               The tuple should be formated as (xbins, ybins).
        """

        self._binOnRead = bins

    def setFilenames(self, filenames = None, darkfilenames = None):
        """Set the list of filenames and darkfilenames

        filenames     : list of strings
           Set the list of filenames to process for the light images
        darkfilenames : list of strings
           Set the list of filenames to process the dark images
        """
        self.filenames = filenames
        self.darkfilenames = darkfilenames

    def setMeanMonitor(self, b):
        """Set if the images are normalized by the mean of the monitor

        b : bool
           If True, normalize to mean of monitor counts."""
        
        self.meanMonitor = b

    def setFromSpec(self, scan, mon = 'Monitor', spool = False, spoolpath = ''):
        """Set the filenames from a SpecScan instance

        This function sets the following from a SpecScan instance:
           
           CCD Filenames
           CCD Dark Filenames
           Monitor normalization values.

        scan    : SpecScan instance
           SpecScan to use for data
        mon     : String
           Name of monitor in SpecScan instance
        spoolpath : String
           Path to use for spool file.

        Note: If scan is a list of spec scans, then the filenames and 
        normalization data will be concatenated to make a single 
        ImageProcessor module from all the data."""
        
        if type(scan) != list:
            scan = [scan]

        self.filenames = []
        self.darkfilenames = []
        self.normData = np.array([])
        for s in scan:
            self.filenames += s.ccdFilenames
            self.darkfilenames += s.ccdDarkFilenames
            self.normData = np.concatenate((self.normData, s.values[mon]))

        # set spool path

        if spool:
            sn = ''.join(['_%d' % s.scanno[0] for s in scan])
            self.setSpoolFile(spoolpath + os.sep + scan[0].datafile.filename.split(os.sep)[-1] + sn)
            self.setSpool()

    def setSpoolFile(self, filename = None):
        """Set the intermediary spool (process) file

        filename : filename (with path) of process file"""

        self.spoolFilename = filename
        

    def setSpoolOut(self, spool = True):
        """Set the fileprocessor to spool the images out to disk"""
        if spool:
            self.spoolfd = open(self.spoolFilename, 'wb')
        else:
            self.spoolfd = None

    def setSpoolIn(self, spool = True):
        """Set the fileprocessor to spool the images from disk"""
        if spool:
            self.spoolfdin = open(self.spoolFilename, 'rb')
        else:
            self.spoolfdin = None

    def setSpool(self, write = False):
        """Set the use of the spool file.

        If the spool file set by FileProcessor.setSpoolFile exists
        then use that to read when used with an iterator. If the
        file doesn't exist then the images are spooled to the given file
        
        Please note this is fragile. This is intended for use to speed up
        internal processing NOT for saving data.

        If write is True then writing is forced"""

        print "---- Spool Filename = %s" % self.spoolFilename
        if self.spoolFilename is not None:
            if os.path.exists(self.spoolFilename):
                self.setSpoolIn(True)
                self.setSpoolOut(False)
                print "---- Reading from spool file"
            else:
                self.setSpoolOut(True)
                self.setSpoolIn(False)
                print "---- Writing to spool file"

    def __iter__(self):
        return self

    def next(self):
        a = self.pop(quiet = True)
        if not a:
            raise StopIteration
        else:
            return self.images

    def pop(self, *args, **kwargs):
        """Process the next image in the sequence

        This routine will call process() on the next image
        in a sequence (stack) of images. It will return true if an image
        was processed and false if we are at the end of the stack.

        For arguments, see the process() function.
        """

        if self.spoolfdin is not None:
            # We have a valid spool file
            return self._readSpoolImage()

        if self.currentFrame == len(self.filenames):
            # We are at the end
            return False
        else:
            self.process(frames = self.currentFrame, *args, **kwargs)
            return True

    def __getitem__(self, index):
        """Convinience function to process and get an image"""
        if self.spoolfdin is not None:
            self._readSpoolImage(index)
        else:
            self.process(frames = index, quiet = True)

        return self.images

    def process(self, dark = True, norm = True, 
                dtype = np.float, quiet = False,
                crop = False, BG = False,
                frames = None,
                keepdark = False):
        """Read in images and process into 3D array.
        
        dark  : bool
           If True, subtract the dark images from the data
        norm  : bool
           If True, normalize by monitor.
        dtype : datatype 
           numpy datatype of processed array.
        quiet : bool
           If True, dont write to screen when reading images.
        frames : list/int
           If not none, only process the frames passed to the 
           variable frames.
        keepdark : bool
           If true keep an array of all dark images. If false
           only store the latest in an array of darkimages"""
   
        images = []
        stdevs = []
        darkstdevs = []

        if not self.processed:
            self.darkimages = []
            self.processed = True
        
        if norm:
            if self.normData is None:
                normData = np.ones(len(self.filenames))
                print "XXXX No normalization data found"
            else:
                if self.meanMonitor:
                    normData = np.ones(len(self.filenames)) * self.normData.mean()
                    if not quiet:
                        print "---- Normalizing data (using mean value)."
                else:
                    normData = self.normData
                    if not quiet:
                        print "---- Normalizing data."
        else:
            normData = np.ones(len(self.filenames))
            
        if dark and not quiet:
            print "---- Correcting for dark current"

        if not quiet:
            print "---- Reading Images"

        if len(self.darkfilenames) == 0:
            self.darkfilenames = self.filenames

        if frames is None:
            self.framesToProcess = range(0, len(self.filenames))
            if not quiet:
                print "---- Processing all frames (%d)" % len(self.framesToProcess)
        else:
            if type(frames) != list:
                frames = [frames]
            self.framesToProcess = frames
            if not quiet:
                print "---- Processing %d frames" % len(self.framesToProcess)

        normIterator = operator.itemgetter(*self.framesToProcess)(normData.tolist())
        if type(normIterator) != tuple:
            normIterator = [normIterator]

        # Loop over all frames to process

        for (i, iname, diname, normVal) in zip(self.framesToProcess,
                                               operator.itemgetter(*self.framesToProcess)(self.filenames), 
                                               operator.itemgetter(*self.framesToProcess)(self.darkfilenames), 
                                               normIterator):

            # Store current frame number

            self.currentFrame = i + 1

            # Now choose if we have multiple acquisitions

            if type(iname) == list:
                _images = None
                _darkimages = None 
                _darkstdev = None
                _stdev = None
                #Start reading the light images
                for j, _in in enumerate(iname): 
                    if os.path.exists(_in):
                        image = self._getRawImage(_in).astype(dtype)
                        _images, _stdev = self._binImageWithStdev(_images, _stdev, image)    
                        
                        if not quiet:
                            print "---- Reading image %-3d of %-3d (sub image %-3d of %-3d)     \r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()
                    else:
                        if not quiet:
                            print "---- Missing image %-3d of %-3d (sub image %-3d of %-3d)\r" % (i + 1, len(self.filenames), j + 1, len(iname)),
                            sys.stdout.flush()

                if dark:
                    for j, _din in enumerate(diname):
                        if os.path.exists(_din):
                            darkimage =  self._getRawImage(_din).astype(dtype)
                            _darkimages, _darkstdev = self._binImageWithStdev(_darkimages, _darkstdev, darkimage)
                        
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
                    if keepdark:
                        self.darkimages.append(darkimage)
                    else:
                        self.darkimages = [darkimage]
            else:
                _darkstdev = None
                _stdev = None
                # Process only single image pair
                image = self._getRawImage(iname).astype(dtype)

                if os.path.exists(diname):
                    darkimage =  self._getRawImage(diname).astype(dtype)
                    if keepdark:
                        self.darkimages.append(darkimage)
                    else:
                        self.darkimages = [darkimage]

                if not quiet:
                    print "---- Reading image %-3d of %-3d\r" % (i, len(self.filenames)),

            if dark:
                if len(self.darkimages):
                    image = image - self.darkimages[-1]
                else:
                    print "XXXX Unable to dark currect correct. No Image found"
            
            if norm:
                image = image / normVal

            images.append(image)

            # Write out the image to the spool file if present

            if self.spoolfd is not None:
                self._writeSpoolImage(image)

            if _stdev is not None:
                stdevs.append(np.sqrt(_stdev[1] / _stdev[2]))
            else:
                stdevs.append(zeros(image.shape))

            if _darkstdev is not None:
                darkstdevs.append(np.sqrt(_darkstdev[1] / _darkstdev[2]))
            else:
                darkstdevs.append(zeros(image.shape))

        if not quiet:
            print "\n---- Processed %d images (%d dark images)" % (len(images), len(self.darkimages))
            
        
        self.images = np.array(images)
        self.stdevs = np.array(stdevs)
        self.darkstdevs = np.array(darkstdevs)
        if not quiet:
            print "---- Done. Array size %s (%.3f Mb)." % (str(self.images.shape), 
                                                           self.images.nbytes / 1024**2) 
    
    def _writeSpoolImage(self, image):
        """Write out spool image"""
        if self.spoolfd is not None:
            if self.spoolfd.tell() == 0:
                # Start of file, write out image size
                self.spoolfd.write("%d,%d,%s\n" % (image.shape[0], image.shape[1], str(image.dtype)))
            scipy.io.numpyio.fwrite(self.spoolfd, image.size, image, 'f')

    def _readSpoolImage(self, imageno = None):
        """Read Spool Image"""
        if self.spoolfdin is not None:
            if self.spoolfdin.tell() == 0:
                # Read Parameters
                d = self.spoolfdin.readline().strip('\n')
                d = d.split(',')
                self.spoolImageSize = (int(d[0]), int(d[1]))
                self.spoolImageType = d[2]
            if imageno is not None:
                # Seek to find place in file
                self.spoolfdin.seek(4 * self.spoolImageSize[0] * self.spoolImageSize[1] * imageno)
            image = scipy.io.numpyio.fread(self.spoolfdin, 
                                           self.spoolImageSize[0] * self.spoolImageSize[1],
                                           'f')
            if len(image):
                image = image.reshape(1,self.spoolImageSize[0],self.spoolImageSize[1])
                self.images = image
                return True
            else:
                return False
        else:
            return False
            
    def _binImageWithStdev(self, images, stdev, newimage):
        """Bin the images and provide a stdev of all bins"""
        newstdev = [None,None,None]

        # Note stdev[0] is M_k and stdev[1] is Q_k
        if stdev is None:
            newstdev = [newimage, zeros(newimage.shape), 1]
        else:
            newstdev[2] = stdev[2] + 1
            newstdev[1] = stdev[1] + ((newstdev[2] - 1) * pow(newimage - stdev[0],2) / newstdev[2])
            newstdev[0] = stdev[0] + ((newimage - stdev[0]) / stdev[2])

        if images is not None:
            return (newimage + images), newstdev
        else:
            return newimage, newstdev

    def _getRawImage(self, iname):
        """Read raw image"""
        
        if self._format == 'SPE':
            img = PrincetonSPEFile(iname).getBinnedData()
        elif self._format == 'LCLS':
            img = LCLSdataformat(iname)
        elif self._format == 'TIFF':
            img = np.array(PILImage.open(iname)).sum(-1)
        else:
            raise Exception("Unknown file format \"%s\"" % self._format)

        if self._cropOnRead is not None:
            img = img[self._cropOnRead[0]:self._cropOnRead[2],
                      self._cropOnRead[1]:self._cropOnRead[3]]
        if self._binOnRead is not None:
            img = ccdutils.binArray(img, self._binOnRead)

        return img
    
    def _computeMeanImage(self):
        N = self.images.shape[0]
        self.mean = self.images.sum(0) / N
        self.stderr = (self.images - self.mean)**2

    def maskImageWithStdev(self, light = None, dark = None, quiet = False):
        """Mask the image with thresholds of the stdev

        light : float
                Threshold (stdev) for light images
        dark  : float
                Threshold (stdev) for dark images
        quiet : bool
                Don't output anything to console
                
        If dark or light are none then the mask will not be 
        computed for that threshold

        Note: The mask values are TRUE is they are BAD"""

        self.mask = np.zeros(self.stdevs.shape, dtype = np.bool)

        if light is not None:
            self.mask = self.mask | (self.stdevs > light)
        if dark is not None:
            self.mask = self.mask | (self.darkstdevs > dark)

        if not quiet:
            print "---- Masked data. (%d masked values)" % (self.mask.sum())

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

    def getMask(self, n = None):
        """Return the mask of the image data

        n : int or None
           Image mask to return. If None, return all masks."""

        if self.mask is None:
            self.mask = np.zeros(self.images.shape, dtype = np.bool)

        if n is None:
            return self.mask
        else:
            return self.mask[n]
        
        

    def saveImage(self, filename, inum = None, dtype = np.float32):
        """Save image to binary file of specified datatype

        filename : string
           Filename to save image stack to.
        inum     : Integer
           Image number to save
        dtype    : numpy data type
           Datatype for images

        """
        if inum is None:
            self.images.astype(dtype).tofile(filename)
        else:
            self.images[n].astype(dtype).tofile(filename)
                               

    def save(self, filename, compressed = False):
        """Save an image sequance to a numpy binary file
        
        filename   : string
           Filename to save to.
        compressed : bool
           If True save as compressed file"""
        
        obj = {'images'   : self.images,
               'normData' : self.normData}

        np.savez(filename, **obj)
        print "**** Image saved to %s" % filename

    def load(self, filename):
        """Load images from previous save operation

        filename : string
           Filename of images previously stored by save()"""
        
        print "**** Loading image from %s" % filename
        obj = np.load(filename)
        self.images = obj['images']
        self.normData = obj['normData']

    def __str__(self):
        """Return text string of FileProcessor stats"""
        

#
# HDF5 FileProcessor
#

class FileProcessorHDF5():
    def __init__(self, filename = None, quiet = False):
        """Initialize Class"""
        self.quiet = quiet
        if filename is not None:
            self.setFilename(filename)
            self._read()

    def setFilename(self, filename):
        """Set the filename of the HDF5 file"""
        if type(filename) != list:
            self.filename = [filename]
        else:
            self.filename = filename
    
    def process(self):
        """Process (read in) the data"""
        self._read()

    def _read(self):
        """Read the HDF5 Container"""
        if not self.quiet:
            print "**** Reading data from %s" % self.filename

        self._h5data = h5py.File(self.filename, 'r')
        self.images = array(self._h5data['entry']['instrument']['detector']['data'])
        
        n = self.images.shape[0]
        
        angles = self._h5data['entry']['instrument']['NDAttributes']
        self.sixcAngles = ones((n,6))
        self.sixcAngles[:,0] = array(angles['Delta'])[0:n]
        self.sixcAngles[:,1] = array(angles['Gamma'])[0:n]
        self.sixcAngles[:,2] = array(angles['Theta'])[0:n]
        self.sixcAngles[:,3] = array(angles['Chi'])[0:n]
        self.sixcAngles[:,4] = array(angles['Phi'])[0:n]
        self.sixcAngles[:,5] = array(angles['Mu'])[0:n]
        self.monitor = angles['Monitor'][:n]
        #self.images = np.inner(self.images, ( 1 / self.monitor ))
        if not self.quiet:
            print "---- Done. Array size %s (%.3f Mb)." % (str(self.images.shape), 
                                                           self.images.nbytes / 1024**2) 
    
    def getImage(self, n = None):
        """Get image data"""
        if n is None:
            return self.images
        else:
            return self.images[n]

    def getSIXCAngles(self):
        """Get the diffractometer setting angles"""
        return self.sixcAngles

#
# image processor class
#

class ImageProcessor():
    """Image Processor class

    This class provides the processing of single and sets of CCD-images
    in more detail: each pixel is transformed into reciprocal space
    the set of reciprocal vectors and intensities is gridded on a regular cuboid
    the needed informations can be provided by a spec scan from the pyspec package"""

    def __init__(self, fP = None, configfile = "ccd.cfg", ccdname = None, spec = None):
        """Initialize the image processor

        fP         : pyspec.ccd.transformations.FileProcessor instance
                     File processor object for getting the images
        configfile : string
                     File defining config of CCD
        ccdname    : string
                     Name of CCD as defined in config file. If none, then 
                     do not process config file
        spec       : pyspec.spec.SpecScan instance
                     Spec Scan used to obtain information on current set"""
        
        # file processor to provied the needed images
             
        # set parameters to configure the CCD setup
        # detector distance 30cm and detector pixel size 20um
        self.setDetectorProp(0.020, 0.020, 1300, 1340, 650, 670)
        self.setDetectorPos(300, 0.0)
        self.setDetectorMask()
        
        # Overide if a config file is used.
        #self.ccdname = ccdname
        #self.ccdConfigFile = configfile
        #if ccdname is not None:
        #    self.readConfigFile(configfile)
        
        # Default to frame mode 1
        self.setFrameMode(1)
        
        self.totSet   = None
        
        self.gridData   = None
        self.gridOccu   = None
        self.gridOut    = None
        self.gridStdErr = None

        self.Qmin  = None
        self.Qmax  = None
        self.dQN   = None

        self._mask = None

        if fP is not None:
            self.setFileProcessor(fP)

    def readConfigFile(self, filename):
        """Read config file and set detector parameters

        This routine will read the detector vital statistics from 
        a config file of name 'filename'. If the config file does not
        exitst in the path, then standard locations are searched.

        filename : string
           filename of config file.

        """
        config = CCDParamsConfigParser()
        fname = config.readAllLocations(filename)
        print "---- Reading config '%s' from %s" % (self.ccdName, fname)
        xsize = config.getInt(self.ccdName, 'xsize', 1300)
        ysize = config.getInt(self.ccdName, 'ysize', 1340)
        pixx  = config.getFloat(self.ccdName, 'pixel_xsize', 0.020)
        pixy  = config.getFloat(self.ccdName, 'pixel_ysize', 0.020)
        xcen  = config.getFloat(self.ccdName, 'xcen', 650.0)
        ycen  = config.getFloat(self.ccdName, 'ycen', 570.0)
        dist  = config.getFloat(self.ccdName, 'sample_det_distance', 300.0)
        rot   = config.getFloat(self.ccdName, 'detector_rotation', 0.0)
        self.setDetectorProp(pixx, pixy, xsize, ysize, xcen, ycen)
        self.setDetectorPos(dist, rot)
        
    #
    # set and get part
    #

    def getGrid(self):
        """Convinience function to return useful grid values

        This function returns the X, Y and Z coordinates of the grid
        along with the intensity and standard error grids. For example::

            >>>X, Y, Z, I, E, N = ip.getGrid()

        """
        X, Y, Z = self.getGridMesh()
        I = self.getGridData()
        E = self.getGridStdErr()
        N = self.getGridOccu()
        return X, Y, Z, I, E, N

    def getImageData(self):
        """Get the totat image data, transformed into Q"""
        return self.totSet

    def getGridData(self):
        """Return the intensity grid

        This function will mask the data if previously
        set to be masked"""
        # Here we check for mask
        if self._mask is not None:
            return ma.masked_array(self.gridData, mask = self._mask)
        else:
            return self.gridData

    def getGridStdErr(self):
        """Return the standard error grid"""
        if self._mask is not None:
            return ma.masked_array(self.gridStdErr, mask = self._mask)
        else:
            return self.gridStdErr

    def getGridStdDev(self):
        """Return the standard deviation grid"""
        st = self.gridStdErr * sqrt(self.gridOccu)
        if self._mask is not None:
            return ma.masked_array(st, mask = self._mask)
        else:
            return st
    
    def getGridOccu(self):
        """Return the occupation of the grid"""
        return self.gridOccu

    def setGridMaskOnOccu(self, thresh = 10):
        """Set grid mask from bin occupation

        thresh : int
            Threshold for grid occupation"""

        self._mask = self.gridOccu <= thresh

    def setDetectorProp(self, detPixSizeX, detPixSizeY, detSizeX, detSizeY, detX0, detY0):
        """Set properties of used detector

        detPixSizeX : float (mm)
             Detector pixel size in detector X-direction.
        detPixSizeY : float (mm)
             Detector pixel size in detector Y-direction.
        detSizeX    : integer
             Detector no. of pixels (size) in detector X-direction.
        detSizeY    : integer
             Detector no. of pixels (size) in detector Y-direction.
        detX0       : float
             Detector X-coordinate of center for reference.
        detY0       : float
             Detector Y-coordinate of center for reference."""
        
        self.detPixSizeX = detPixSizeX  
        self.detPixSizeY = detPixSizeY
        self.detSizeX    = detSizeX
        self.detSizeY    = detSizeY
        self.detX0       = detX0
        self.detY0       = detY0

    def getDetectorProp(self):
        """Get properties of used detector returned as a tuple

        detPixSizeX : float (mm)
             Detector pixel size in detector X-direction.
        detPixSizeY : float (mm)
             Detector pixel size in detector Y-direction.
        detSizeX    : integer
             Detector no. of pixels (size) in detector X-direction.
        detSizeY    : integer
             Detector no. of pixels (size) in detector Y-direction.
        detX0       : float
             Detector X-coordinate of center for reference.
        detY0       : float
             Detector Y-coordinate of center for reference."""

        return self.detPixSizeX, self.detPixSizeY, self.detSizeX, self.detSizeY, self.detX0, self.detY0

    def setDetectorPos(self, detDis = 300.0, detAng = 0.0):
        """Set the detector position

        detDis : float (mm)
           Detector distance from sample.
        detAng : float (deg)
           Detector miss alignement angle (rotation around incident beam)"""
        
        self.detDis = detDis
        self.detAngle = detAng

    def getDetectorPos(self, detDis = 30.0, detAng = 0.0):
        """Get the detector position

        detDis : float (mm)
           Detector distance from sample.
        detAng : float (deg)
           Detector miss alignement angle (rotation around incident beam)"""
        
        return self.detDis, self.detAngle

    def setDetectorMask(self, mask = None):
        """Set the detector mask, to exclude data points

        mask : ndarray
           Array of same form as CCD images, with 1's for valid data points

        """
        self.detMask = mask

    def setBins(self, binX = 1, binY = 1):
        """Set detector binning.

        This function sets the detector binning, and modifies the detector paremeters accordingly.

        binX : integer
           Binning in x-direction.
        binY : integer
           Binning in y-direction.
           
        """
    
        # try to get old bins
        try:
            oldBinX = self.binX
        except:
            oldBinX = 1
        try:
            oldBinY = self.binY
        except:
            oldBinY = 1

        # scaling ratios
        ratX = 1.0*binX/oldBinX
        ratY = 1.0*binY/oldBinY

        # apply changes to detector probperties
        self.detPixSizeX *= ratX
        self.detPixSizeY *= ratY
        self.detSizeX     = int(self.detSizeX / ratX)
        self.detSizeY     = int(self.detSizeY / ratY)
        self.detX0       /= ratX
        self.detY0       /= ratY

        self.binX = binX
        self.binY = binY
        
    def getBins(self, binX, binY):
        """Set no. of bins. Takes them into acount for pixel size, detector size and detector center
        
        binX : no. of pixels along detector X-direction which are bined
        binY : no. of pixels along detector Y-direction which are bined"""

        return self.binX, self.binY

    def setSetSettings(self, waveLen, settingAngles, UBmat):
        """Set the settings for the set 

        The set settings are:
        waveLen       : float
           Wavelength of the X-rays (Angstroms).
        settingAngles : setting angles of the diffractometer at each image (data point)
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        """

        self.waveLen       = waveLen
        self.settingAngles = settingAngles
        self.UBmat         = UBmat

    def setSpecScan(self, conScan):
        """Set the settings for the set from the considered pyspec scan object

        The set settings are:
        waveLen       : wavelength of the X-rays  (float in Angstrom)
        energy        : photon energy             (float in eV)
        settingAngles : setting angles of the diffractometer at each image (data point)
        UBmat         : UB matrix (orientation matrix) to transform the HKL-values into the sample-frame (phi-frame)
        
        Note: If the scans are a list, then they will be constructed from that list from all angles.
        The first scan will be used for the wavelength and UB matrix."""

        self.conScan = conScan
        self._setFromSpecScan()
   
    def getSpecScan(self):
        """Get the pyspec scan object which was used for the set settings
        
        returns pyspec scan object which contains the needed information"""

        return self.conScan

    def setFrameMode(self, mode):
        """Set the mode of the output frame for (Qx, Qy, Qz)

        mode : Integer
           Mode (see below)

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

    def getFrameMode(self):
        """Get the mode of the output frame for (Qx, Qy, Qz)

        mode 1 : 'theta'    : Theta axis frame.  
        mode 2 : 'phi'      : Phi axis frame.
        mode 3 : 'cart'     : Crystal cartesian frame.
        mode 4 : 'hkl'      : Reciproal lattice units frame."""

        return self.frameMode

    def setGridOptions(self, *args, **kwargs):
        """Depreciated function

        .. warning::
           This function is depreciated. See :func:setGridSize()
        """
        raise Exception("Depriciated Function, use : setGridSize()")

    def getGridOptions(self, *args, **kwargs):
        """Depreciated function

        .. warning::
           This function is depreciated. See :func:getGridSize()
        """
        raise Exception("Depriciated Function, use : getGridSize()")

    def setGridSize(self, Qmin = None, Qmax = None, dQN = None):
        """Set the options for the gridding of the dataset

        Qmin : ndarray
           Minimum values of the cuboid [Qx, Qy, Qz]_min
        Qmax : ndarray
           Maximum values of the cuboid [Qx, Qy, Qz]_max
        dQN  : ndarray
           No. of grid parts (bins)     [Nqx, Nqy, Nqz]"""

        if Qmin is not None:
            self.Qmin = np.array(Qmin)
        if Qmax is not None:
            self.Qmax = np.array(Qmax)
        if dQN is not None:
            self.dQN  = np.array(dQN)

    def getGridSize(self):
        """Get the options for the gridding of the dataset

        returns tuple of (Qmin, Qmax, dQN)"""

        return self.Qmin, self.Qmax, self.dQN

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

    #
    # get set functions for input output
    #

    def setFileProcessor(self, fp = None):
        """Set the FileProcessor object for treating the CCD images

        fp : FileProcessor instance
           FileProcessor to use to obtain CCD images"""
        
        self.fileProcessor = fp

    def getFileProcessor(self):
        """Get the FileProcessor object for treating the CCD images

        fp : FileProcessor object with the CCD images"""
        
        return self.fileProcessor

    #
    # help function part
    #

    def _calcBMatrix(self, angles):
        """Calculate the B matrix from reciprocal space angles

        anges: ndarray
           Array of real space values, [a, b, c, alpha, beta, gamma]

        returns B matrix as (3x3) ndarray.
        """
        B = np.ones((3,3))
        B[0,0] = angles[0]
        B[1,0] = 0
        B[2,0] = 0
        B[0,1] = angles[1] * cos(angles[5])
        B[1,1] = angles[1] * sin(angles[5])
        B[2,1] = 0
        B[0,2] = angles[2] * cos(angles[4])
        B[1,2] = -1.0 * angles[2] * sin(angles[4]) * cos(angles[3])
        B[2,2] = 2 * np.pi / angles[2]

        return B
    
    def _setFromSpecScan(self):
        """See setFromSpecScan"""

        if type(self.conScan) != list:
            scan = [self.conScan]
        else:
            scan = self.conScan

        self.waveLen       = scan[0].wavelength  # in Angstrom
        self.energy        = Diffractometer.hc_over_e / scan[0].wavelength # in eV
        self.settingAngles = scan[0].getSIXCAngles()
        for s in scan[1:]:
            self.settingAngles = np.concatenate((self.settingAngles,s.getSIXCAngles()))
        
        self.UBmat         = scan[0].UB
        

    def _setFromFileProcessor(self):
        """Set Values from the ImageProcessor"""
        self.waveLen = self.fileProcessor.getWavelength()
        self.energy = Diffractometer.hc_over_e / self.waveLen
        self.settingAngles = self.fileProcessor.getSIXCAngles()
        self.UBMat = self.fileProcessor.getUB()
        

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

    #
    # help functions for input / output
    #

    def __str__(self):
        """Create the information about the grid process"""
 
        _s = "Detector Parameters:\n\n"

        _s += "Detector size        :%d x %d [pixels]\n" % (self.detSizeX, self.detSizeY)
        _s += "Detector pixel size  :%f x %f [mm]\n" % (self.detPixSizeX, self.detPixSizeY)
        _s += "Zero (center) of CCD :(%f, %f) [pixels]\n" % (self.detX0, self.detY0)
        _s += "Sample Detector dist.:%f [mm]\n" % self.detDis
        _s += "Detector rot. ang.   :%f [deg]\n" % self.detAngle

        if self.totSet is not None:
            _s += "\n\nData Set:\n" 

            _s += "Number of Pixels : %.2e\n" % self.totSet.shape[0]
            _s += 'Q_min            : [%.3e %.3e %.3e]\n' % (self.totSet[:,0].min(), 
                                                             self.totSet[:,1].min(), 
                                                             self.totSet[:,2].min())
            _s += 'Q_max            : [%.3e %.3e %.3e]\n' % (self.totSet[:,0].max(), 
                                                             self.totSet[:,1].max(), 
                                                             self.totSet[:,2].max())
            _s += 'I_min            : %.3e\n' % self.totSet[:,3].min()
            _s += 'I_max            : %.3e\n' % self.totSet[:,3].max()
            _s += 'I_ave            : %.3e\n' % self.totSet[:,3].mean()

        _s += "\n\nGrid Parameters:\n\n"
        mode = self.frameMode
        if mode == 1:
            _s += 'Frame Mode 1: (Qx, Qy, Qz) in theta-frame and (1/Angstrom)\n' 
        elif mode == 2:
            _s += 'Frame Mode 2: (Qx, Qy, Qz) in phi-frame and (1/Angstrom)\n'
        elif mode == 3:
            _s += 'Frame Mode 3: (Qx, Qy, Qz) in cartesian-frame and (1/Angstrom)\n'
        elif mode == 4:
            _s += 'Frame Mode 4: (H, K, L) in hkl-frame and (reciprocal lattice units)\n'

        _s += '\n\nGrid Vectors:\n\n'
        if self.Qmin is not None:
            _s += 'Q_min     : [%.3e %.3e %.3e]\n' % (self.Qmin[0], self.Qmin[1], self.Qmin[2])
        if self.Qmax is not None:
            _s += 'Q_max     : [%.3e %.3e %.3e]\n' % (self.Qmax[0], self.Qmax[1], self.Qmax[2])
        if self.dQN is not None:
            _s += 'Grid Size : [%.3e %.3e %.3e]\n' % (self.dQN[0], self.dQN[1], self.dQN[2])

        if self.gridData is not None:
            _s += '\n\nGrid Statistics:\n\n'
            _s += 'No. of voxels in grid           : \t %.2e\n' % (self.gridData.size)
            _s += 'No. of data points outside grid : \t %.2e\n' % (self.gridOut)
            _s += 'No. of empty voxels             : \t %.2e\n' % ((self.gridOccu == 0).sum())
            
        return _s

    #
    # process part
    #
    
    def processToQ(self):
        """Process selcted images of the full set into (Qx, Qy, Qz, I)

        This function takes the images provided by a FileProcessor instance, and
        converts them to a set of (Qx, Qy, Qz, I) values in the frame mode which
        is set prevously in this class.

        """
        ccdToQkwArgs = {}
        if self.totSet is not None:
            del self.totSet
            gc.collect()

        if not self.fileProcessor:
            raise Exception("No FileProcessor specified.")

        # if images not yet processed, do it
        if getattr(self.fileProcessor, 'images', None) is None:
            self.fileProcessor.process()

        if self.settingAngles is None:
            raise Exception("No setting angles specified.")

        print "\n"
        print "---- Setting angle size :", self.settingAngles.shape
        print "---- CCD Size :", (self.detSizeX, self.detSizeY)
        print "**** Converting to Q"
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
        print "---- Setsize is %d" % self.totSet.shape[0]
        self.totSet[:,3] = np.ravel(self.fileProcessor.getImage())  
       
        if self.detMask is None:
            m = self.fileProcessor.getMask()
            if m is not None:
                if m.sum():
                    self.detMask = m 

        if self.detMask is not None:
            print "---- Masking data"
            totMask = self.detMask.ravel()
            self.totSet = self.totSet[totMask == False,:]
            
            print "---- New setsize is %d " % self.totSet.shape[0]
            
                    
    def process(self):
        """Convinience function to perform all processing"""
        self.processToQ()
        self.processGrid()

    def processGrid(self):
        """Process imageset of (Qx, Qy, Qz, I) into grid 

        This function, process the set of (Qx, Qy, Qz, I) values and grid the data."""
        print "---- Total data is %f MBytes\n" % (self.totSet.nbytes / 1024.0**2)
        
        if self.totSet is None:
            raise Exception("No set of (Qx, Qy, Qz, I). Cannot process grid.")

        # prepare min, max,... from defaults if not set
        if self.Qmin is None:
            self.Qmin = np.array([ self.totSet[:,0].min(), self.totSet[:,1].min(), self.totSet[:,2].min() ])
        if self.Qmax is None:
            self.Qmax = np.array([ self.totSet[:,0].max(), self.totSet[:,1].max(), self.totSet[:,2].max() ])
        if self.dQN is None:
            self.dQN = [100, 100, 100]

        # 3D grid of the data set 
        print "**** Gridding Data."
        t1 = time.time()
        gridData, gridOccu, gridStdErr, gridOut = ctrans.grid3d(self.totSet, self.Qmin, self.Qmax, self.dQN, norm = 1)
        t2 = time.time()
        print "---- DONE (Processed in %f seconds)" % (t2 - t1)
        emptNb = (gridOccu == 0).sum()
        if gridOut != 0:
            print "---- Warning : There are %.2e points outside the grid (%.2e bins in the grid)" % (gridOut, gridData.size)
        if emptNb:
            print "---- Warning : There are %.2e values zero in the grid" % emptNb

       # store intensity, occupation and no. of outside data points of the grid
        self.gridData   = gridData
        self.gridOccu   = gridOccu
        self.gridOut    = gridOut
        self.gridStdErr = gridStdErr

    def makeGridData(self, *args, **kwargs):
        """Convinience to call makeGrid"""
        self.process(*args, **kwargs)


##################################################################################################################
##                                                                                                              ##
##  GridProcessorClass to perform cutting of grid                                                               ##
##                                                                                                              ##
##################################################################################################################

class GridProcessor():
    """Class to process grid data

    This class is used to process grid data, and perform line and sum integration from the 
    gridded data. 

    This class assumes that you have processed the data using an ImageProcessor and requires 
    the use of that when initializing the class."""

    def __init__(ip = None):
        if ip is None:
            raise Exception("This class must be initialized with a valid ImageProcessor")
        else:
            self.ip = ip

    def get1DCut(pos, axis):
        """Make a 1D cut of the grid along a principal axis."""
        ax = self._processAxis()

    def _processAxis(axis):
        if type(axis) == int:
            return axis
        elif type(axis) == str:
            if axis.upper() == "X":
                return 0
            elif axis.upper() == "Y":
                return 1
            elif axis.upper() == "Z":
                return 2
            else:
                raise Exception("Invalid string. Axis must be 'X', 'Y', or 'Z'.")
        else:
            raise Exception("Invalid type in axis. Must be integer or string")
        return Null
    
##################################################################################################################
##                                                                                                              ##
##  CCD Config Parser, to read config file from disk                                                            ##
##                                                                                                              ##
##################################################################################################################        
        
    
class CCDParamsConfigParser(ConfigParser):
    """Class to read config file which defines all CCD parameters"""
    def readAllLocations(self, filename):
        
        if not os.path.isfile(filename):
            locations = []
            if os.name is 'posix':
                if os.environ.has_key('HOME'):
                    locations.append(os.environ['HOME'] + os.path.sep + ".pyspec")
                locations.append('/usr/local/pyspec/etc')
                locations.append('/etc/pyspec')
            
            for l in locations:
                _f = l + os.path.sep + filename
                if os.path.isfile(_f):
                    self.read(_f)
                    return _f
        else:
            self.read(filename)
            return filename

        return None
            
    def getWithDefault(self,section, option, default):
        try:
            r = self.get(section, option)
        except ConfigParser.NoOptionError:
            r = default
            print "**** Using default value for %s." % option
        return r
    
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
        return self._getWithConvert(float, *args, **kwargs)
    def getInt(self, *args, **kwargs):
        return self._getWithConvert(int, *args, **kwargs)
    def getFloat(self, *args, **kwargs):
        return self._getWithConvert(float, *args, **kwargs)

        
####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    from pyspec import spec
    
    sf   = spec.SpecDataFile('/home/tardis/spartzsch/data/ymn2o5_oct10/ymn2o5_oct10_1', 
    			     ccdpath = '/mounts/davros/nasshare/images/oct10')
    scan = sf[360]

    fp = FileProcessor(spec = scan)
    fp.process()

    testData = ImageProcessor(fp)
    #testData.setDetectorPos(detAng = -1.24)
    testData.setBins(4, 4)
    testData.setSpecScan(scan)
    testData.setFrameMode(4)

    testData.process()

    print testData

    

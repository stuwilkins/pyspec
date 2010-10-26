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
from   pyspec import spec
from   pyspec.diffractometer import Diffractometer
from   pyspec.ccd.PrincetonSPE import *
from   pyspec.ccd.plotter import PlotGrid
import gridder
from ConfigParser import ConfigParser

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
        # image treatment
        self.conRoi      = None
        self.frameMode   = 1
	self.setName     = 'Set #'
        self.setNum      = 1
        # gridder options
        self.setGridOptions()
        # plot options
        self.axesLabels  = ['Qx', 'Qy', 'Qz']
        self.plotFlag2D = 7
        self.plotFlag1D = 7
        self.logFlag1D  = 0
        self.logFlag2D  = 0
        self.fit1D      = False
        self.fitType    = 'lor2a'
        self.histBin    = 50


    def readConfigFile(self, filename):
        config = MyConfigParser()
        config.read(filename)
        config.get(self.ccdname, '', vars = {})
        
    #
    # set part
    #

    def setBins(self, binX, binY):
        """Takes binning into acount"""
        self.detPixSizeX *= binX
        self.detPixSizeY *= binY
        self.detSizeX    /= binX
        self.detSizeY    /= binY
        self.detX0       /= binX
        self.detY0       /= binY

    def setBySpec(self, conScan):
        """get settings from the considered pyspec scan object
        wavelength, filenames, angels, kind of set (scan)"""

        self.waveLen       = conScan.wavelength
        self.imFileNames   = conScan.getCCDFilenames()
        self.darkFileNames = conScan.getCCDFilenames(dark = 1)
        self.settingAngles = conScan.getSIXCAngles()
        self.intentNorm    = conScan.Ring
        self.UBmat         = conScan.UB
        self.setName       = 'Scan #'
        self.setNum        = conScan.scan
        self.setSize       = self.settingAngles.shape[0]
   
    def setConRoi(self, conRoi):
        """Sets the considered region of interest [xMin, xStep, yMin, yStep]"""
        self.conRoi = conRoi

    def setFrameMode(self, mode):
        """modes of frames are: theta- (1), phi- (2), cartesian- (3) and hkl-frame (4)"""
        self.frameMode = mode

    def setGridOptions(self, Qmin = None, Qmax = None, dQN = None):
        """Sets the options for the gridding of the dataset
        Size of the cuboid: Qmin, Qmax = [Qx, Qy, Qz]_min, max
        Number of parts   : dQN = [Nqx, Nqy, Nqz]"""

        self.Qmin = Qmin
        self.Qmax = Qmax
        self.dQN  = dQN

    #
    # plot settings
    #

    def setAxesLabels(self, axesLabels):
        """Sets the plotting labels for the axes"""
        self.axesLabels = axesLabels

    def setPlotFlags(self, flag2D = 7, flag1D = 7):
        """Sets the ploting flags for 2D and 1D plots
        binary code, flag & 1: intensity, flag & 2: missing grid parts,
        flag & 4: histogram of occupation of the grid parts"""
        self.plotFlag2D = flag2D
        self.plotFlag1D = flag1D

    def setLogFlags(self, flag1D = 0, flag2D = 0):
        """Sets whether data are plotted on linear (False) or logarithmic (True) scale
        binary code, flag & 1: intensity, flag & 2: missing grid parts"""
        self.logFlag1D = flag1D
        self.logFlag2D = flag2D
        
    def setFit1D(self, fit1D = 0, fitType = 'lor2a'):
        """Sets whether 1D lines get fitted by a Loarenzian squared"""
        self.fit1D   = fit1D
        self.fitType = fitType

    def setHistBin(self, histBin = 50):
        """Sets the number of bins for the histograms"""
        self.histBin

    def setPlotIm(self, plotImSelect = None, plotImHor = 4, plotImVer = 3):
        """Sets the options for ploting the raw images"""
        if plotImSelect == None:
            self.plotImSelect = range(self.setSize)
        else:
            self.plotImSelect = plotImSelect
        self.plotImHor    = plotImHor
        self.plotImVer    = plotImVer
        
    #
    # help function part
    #

    def _readImage(self, imNum):
        """Read in the considered region of interest of the image
        dark image subtraction and normalization by ring current"""

        if imNum % 10 == 0:
            print '%s%d image #%03d: read image' % (self.setName, self.setNum, imNum)

        # get dark image (first CCD image)
        darkMon  = self.intentNorm[0]
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
        conDark  = getAreaSet(darkVal,  self.conRoi) / darkMon
        conIm    = getAreaSet(pointVal, self.conRoi) / pointMon - conDark

        return conIm

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

    def _make1DSum(self):
        """1D Lines of the grid data and occupations by summing in the other directions"""
        gridData1DSum = [self.gridData.sum(1).sum(1),
                         self.gridData.sum(0).sum(1),
                         self.gridData.sum(0).sum(0)]
        gridOccu1DSum = [self.gridOccu.sum(1).sum(1),
                         self.gridOccu.sum(0).sum(1),
                         self.gridOccu.sum(0).sum(0)]

        return gridData1DSum, gridOccu1DSum

    def _make2DSum(self):
        """2D Areas of the grid data and occupations by summing in the other direction"""
        gridData2DSum = [self.gridData.sum(0), self.gridData.sum(1), self.gridData.sum(2)]
        gridOccu2DSum = [self.gridOccu.sum(0), self.gridOccu.sum(1), self.gridOccu.sum(2)]

        return gridData2DSum, gridOccu2DSum
    
    def _make1DCut(self):
        """1D Lines of the grid data and occupations at the position of the maximum intensity"""
        gridData1DCut = [self.gridData[:,self.maxInd[1],self.maxInd[2]],
                         self.gridData[self.maxInd[0],:,self.maxInd[2]],
                         self.gridData[self.maxInd[0],self.maxInd[1],:]]
        gridOccu1DCut = [self.gridOccu[:,self.maxInd[1],self.maxInd[2]],
                         self.gridOccu[self.maxInd[0],:,self.maxInd[2]],
                         self.gridOccu[self.maxInd[0],self.maxInd[1],:]]

        return gridData1DCut, gridOccu1DCut
    
    def _make2DCut(self):
        """2D Areas of the grid data and occupations at the position of the maximum intensity"""
        gridData2DCut = [self.gridData[self.maxInd[0],:,:],
                         self.gridData[:,self.maxInd[1],:],
                         self.gridData[:,:,self.maxInd[2]]]
        gridOccu2DCut = [self.gridOccu[self.maxInd[0],:,:],
                         self.gridOccu[:,self.maxInd[1],:],
                         self.gridOccu[:,:,self.maxInd[2]]]

        return gridData2DCut, gridOccu2DCut
    
    #
    # process part
    #

    def processOneImage(self, imNum):
        """Process one image to (Qx, Qy, Qz, I)
        modes of frames are: theta- (1), phi- (2), cartesian- (3) and hkl-frame (4)"""

        # used mode
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

    def processOneSet(self):
        """Process full set of images (Qx, Qy, Qz, I)
        modes of frames are: theta- (1), phi- (2), cartesian- (3) and hkl-frame (4)"""

        print '\n%s%d: process to (Qx, Qy, Qz, I)' % (self.setName, self.setNum)

        # prepare size of full dataset and get data of first scan
        setSize  = self.setSize
        totFirst = self.processOneImage(0)
        imSize   = totFirst.shape[0]
        npts     = setSize * imSize
        totSet   = np.zeros((npts, 4))

        # go through all images and get data sets
        j = 0
        k = imSize
        totSet[j:k,:] = totFirst
        for i in range(1, setSize):
            j = j + imSize
            k = k + imSize
            totSet[j:k, :] = self.processOneImage(i)

        return totSet

    def makeGridData(self, totData):
        """Grid the data set into a cuboid
        Size of the cuboid: Qmin, Qmax = [Qx, Qy, Qz]_min, max
        Number of parts   : dQN = [Nqx, Nqy, Nqz]"""

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
        gridData, gridOccu, gridOut = gridder.grid3d(totData, Qmin[0], Qmax[0], dQN[0], Qmin[1], Qmax[1], dQN[1], Qmin[2], Qmax[2], dQN[2])
        if gridOut != 0:
            print "Warning : There are %.2e points outside the grid (%.2e points in the grid)" % (gridOut, gridData.size - gridOut)
        emptNb = (gridOccu == 0).sum()
        if emptNb:
            print "Warning : There are %.2e values zero in the grid" % emptNb

        # mask the gridded data set
        gridData = np.ma.array(gridData / gridOccu, mask = (gridOccu == 0))
      
        # store intensity, occupation and no. of outside data points of the grid
        self.gridData = gridData
        self.gridOccu = gridOccu
        self.gridOut  = gridOut

        # calculated the corresponding vectors and maximum intensity position of the grid
        self._calcVecDataSet()
        self._calcMax()
        
    #
    # plot part
    #

    def plotGrid1DSum(self):
        """Plots the 1D Lines of the data grid summed over the other dimensions"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag1D = 7)
        gridPlot.setLogFlags(flag1D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot1DAxes(self.qVal, self.axesLabels)
        gridData1DSum, gridOccu1DSum = self._make1DSum()
        gridPlot.setPlot1DData(gridData1DSum, gridOccu1DSum,
                               plotTitle = '1D Lines, over other directions is summed')
        # plot, get figure and axes back
        fig1, allax1 = gridPlot.plot1DData()

    def plotGrid2DSum(self):
        """Plots the 2D Areas of the data grid summed over the other dimension"""

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
        gridData2DSum, gridOccu2DSum = self._make2DSum()
        gridPlot.setPlot2DData(gridData2DSum, gridOccu2DSum,
                               plotTitle = '2D Areas, over other direction is summed')
        # plot, get figure and axes back
        fig2, allax2 = gridPlot.plot2DData()

    def plotGrid1DCut(self):
        """Plots the 1D Lines of the data grid summed over the other dimensions"""

        gridPlot = PlotGrid()
        # flag options and no. of bins for histogram
        gridPlot.setPlotFlags(flag1D = 7)
        gridPlot.setLogFlags(flag1D = 0)
        gridPlot.setHistBin(20)
        # axes and data configuration
        gridPlot.setPlot1DAxes(self.qVal, self.axesLabels)
        gridData1DCut, gridOccu1DCut = self._make1DCut()
        gridPlot.setPlot1DData(gridData1DCut, gridOccu1DCut,
                               plotTitle = '1D Line Cuts at Maximum Position')
        # plot, get figure and axes back
        fig1, allax1 = gridPlot.plot1DData()

        # fit the 1D intensities
        if self.fit1D == True:
            for i in range(gridData1DCut.shape(0)):
                f = fit.fit(x = self.qVal[i], y = gridData1DCut[i], func = [fitfuncs.linear, getatr(fitfuncs, self.fitType)])
                f.go()
                res = f.result

    def plotGrid2DCut(self):
        """Plots the 2D Areas of the data grid summed over the other dimension"""

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
        gridData2DCut, gridOccu2DCut = self._make2DCut()
        gridPlot.setPlot2DData(gridData2DCut, gridOccu2DCut,
                               plotTitle = '2D Area Cuts at Maximum Position')
        # plot, get figure and axes back
        fig2, allax2 = gridPlot.plot2DData()

    def plotImages(self):
        """Plots the selcted images"""
        
        # prepare plots
        plotImNum = self.plotImHor * self.plotImVer
        j = 0
        allfig = []
        allax  = []
        plotImTitle  = '%s%s' % (self.setName, self.setNum)
        plotImExtent = [self.conRoi[0], self.conRoi[0] + self.conRoi[1],
                        self.conRoi[2], self.conRoi[2] + self.conRoi[3]]

        # go through images numbers which should be plotted
        for i in self.plotImSelect:

            # label for y-axis
            yLabel = 'image # %d' % (i)

            if j%plotImNum == 0:
                # prepare plot window
                fig = plt.figure()
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
    testData.setBins(4, 4)
    testData.setBySpec(scan)
    testData.setConRoi([1, 325, 1, 335])
    testData.setFrameMode(4)
    testData.setGridOptions(Qmin = None, Qmax = None, dQN = [90, 160, 30])
    print testData.processOneImage(40)
    #totSet = testData.processOneSet()
    #testData.makeGridData(totSet)
    # plot options
    testData.setAxesLabels([ur"Qx (\u00c5$^{-1}$)", ur"Qy (\u00c5$^{-1}$)", ur"Qz (\u00c5$^{-1}$)"])
    testData.setPlotFlags(7, 7)
    testData.setLogFlags(0, 3)
    testData.setFit1D(False)
    testData.setHistBin(50)
    #testData.plotGrid1DSum()
    #testData.plotGrid2DSum()
    #testData.plotGrid1DCut()
    #testData.plotGrid2DCut()
    #testData.plotAll()

    #testData.setPlotIm(plotImSelect = None, plotImHor = 4, plotImVer = 3)
    #testData.plotImages()
    #plt.show()


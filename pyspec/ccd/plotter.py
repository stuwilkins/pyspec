#
# plotter.py (c) Stuart B. Wilkins 2010 and (c) Sven Partzsch
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

import numpy as np
import matplotlib.pyplot as plt
from   matplotlib.colors import LogNorm
from   pyspec  import fit, fitfuncs

"""
# Try improting 3d routines
try:
    from enthought.mayavi import mlab
except:
    warnings.warn("*** No Enthought mayavi module, 3D visualization is disabled ***")
    pass

# Try importing 3d routines from matplotlib
try:
    from mpl_toolkits import mplot3d as plt3d
except:
    warnings.warn("*** No matplotlib mplot3d module, 3D figures is disabled ***")
    pass
"""

__version__   = "$Revision$"
__author__    = "Stuart B. Wilkins <stuwilkins@mac.com>"
__date__      = "$LastChangedDate$"
__id__        = "$Id$"

class CCDPlot():
    def __init__(self):
        pass

    

class PlotGrid3D():
    def __init__(self, imProc = None):
        self.imProc = imProc
    def plot3D(self):
        data = self.imProc.gridData
        ma = data.max()
        mb = data.min()
        hm = (np.arange(0.7, 0.95, 0.3) * (ma - mb) + mb).tolist()
        Qmin, Qmax, dQN = self.imProc.getGridOptions()
        extent = []
        for i in range(3):
            extent.append(Qmin[i])
            extent.append(Qmax[i])
        print extent
        mlab.contour3d(data, contours = hm, vmax = ma,
                       transparent = True, vmin = mb,
                       extent = extent)


class PlotImages():
    """Plot CCD-images"""

    def __init__(self, imProc):
        # image Processor to get all the needed data
        self._imProc = imProc
        # image selection
        self._plotSelect = None
        # Figure layout
        self._figSize  = (11, 8.5)
        self._plotHor  = 4
        self._plotVer  = 3
        self._plotOrd  = 'hv'
        # plotting intensities / histograms
        self._plotFlag = 1
        # logarithmic plot flag
        self._logFlag  = 0
        # No. of bins for histograms
        self._histBin  = 50

    #
    # set/get part
    #

    def setPlotSelect(self, plotSelect = None):
        """Set the selcected images

        plotSelect : list with the raw images which will be plotted, all if None"""
        
        if plotSelect == None:
            self._plotSelect = range(self._imProc.setSize)
        else:
            self._plotSelect = plotSelect
    
    def getPlotSelect(self):
        """Get the selcected images

        plotSelect : list with the raw images which will be plotted, all if None"""
        
        return self._plotSelect

    def setPlotLayout(self, figSize = (11, 8.5), plotHor = 4, plotVer = 3, plotOrd = 'hv'):
        """Set the options for the ploting window

        figSize : figure size (width, height) in inches, e.g. (11, 8.5)
        plotHor : no. of horizontal images per window  , e.g. 4
        plotVer : no. of vertical   images per window  , e.g. 3
        plotOrd : order of plotting, horizontal-vertical ('hv') of vertical-horizontal ('vh')"""
        
        self._figSize = figSize
        self._plotHor = plotHor
        self._plotVer = plotVer
        self._plotOrd = plotOrd
        
    def getPlotLayout(self):
        """Get the options for the ploting window

        figSize : figure size (width, height) in inches, e.g. (11, 8.5)
        plotHor : no. of horizontal images per window  , e.g. 4
        plotVer : no. of vertical   images per window  , e.g. 3
        plotOrd : order of plotting, horizontal-vertical ('hv') of vertical-horizontal ('vh')"""
        
        return self._figSize, self._plotHor, self._plotVer, self._plotOrd

    def setPlotFlag(self, flag = 1):
        """Set the ploting flag

        flag : flag to select plots
      
        binary code, flag & 1: intensity, flag & 2: histogram"""
        
        self._plotFlag = flag
    
    def getPlotFlag(self):
        """Get the ploting flag

        flag : flag to select plots
      
        binary code, flag & 1: intensity, flag & 2: histogram"""
        
        return self._plotFlag

    def setLogFlag(self, flag = 0):
        """Set whether data are plotted on linear (0) or logarithmic (1) scale

        flag : flag to select between linear and logarithmic plotting

        binary code, flag & 1: intensity, flag & 2: histogram"""
        
        self._logFlag = flag
        
    def getLogFlag(self):
        """Get whether data are plotted on linear (0) or logarithmic (1) scale

        flag : flag to select between linear and logarithmic plotting

        binary code, flag & 1: intensity, flag & 2: histogram"""
        
        return self._logFlag

    def setHistBin(self, histBin = 50):
        """Set the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        self._histBin = histBin

    def getHistBin(self):
        """Get the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        return self._histBin

    #
    # plot functions
    #

    def _plot2DData(self):
        """Plots the selcted images and histograms

        retrurns
        allfig : list of plt.figure objects of the plotting windows
        allax  : list of plt.axes objects which carry the figures"""
        
        # prepare plots
        plotNum = self._plotHor * self._plotVer
        # plotting order
        if self._plotOrd == 'vh':
            axNum = np.ravel( np.array(range(plotNum)).reshape(self._plotVer, self._plotHor).T ) + 1
        else:
            axNum = np.array(range(plotNum)) + 1
        j = 0
        allfig = []
        allax  = []

        # go through images numbers which should be plotted
        for i in range(self._totImNum):

            # considered image
            image  = self._images[i]
            # image No. as title of each image
            pTitle = self._pTitles[i] 
            
            # prepare plot window
            if j%plotNum == 0:
                fig = plt.figure(figsize = self._figSize)
                fig.suptitle(self._plotTitle, fontsize = 24)
                allfig.append(fig)

            # intensities
            if self._plotFlag & 1:
                # new subplot
                ax = plt.subplot(self._plotVer, self._plotHor, axNum[j%plotNum])
                allax.append(ax)
                # plot the image
                if self._logFlag & 1 :
                    minPosVal = (image > 0.0).min()
                    negPart   = image <= 0
                    negPart   = np.ones(negPart.shape) * minPosVal
                    cax  = ax.imshow(image, norm=LogNorm(), extent = self._plotExtent)
                else:
                    cax = ax.imshow(image, extent = self._plotExtent)
                # add a colorbar
                fig.colorbar(cax, format = '%.1e')
                # make the shape a square
                ax.set_aspect(1./ax.get_data_ratio())
                # show the image with number as title    
                ax.set_title(pTitle, fontsize = 18)
                
                # increment the plot image counter
                j += 1

            # prepare plot window
            if j%plotNum == 0:
                fig = plt.figure(figsize = self._figSize)
                fig.suptitle(self._plotTitle, fontsize = 24)
                allfig.append(fig)

            # histograms
            if self._plotFlag & 2:
                # new subplot
                ax = plt.subplot(self._plotVer, self._plotHor, axNum[j%plotNum])
                allax.append(ax)
                # plot hsitogram of the image
                ax.hist(np.ravel(image), self._histBin, log=self._logFlag & 2, facecolor='green')
                # show the image with number as title
                ax.set_title(pTitle, fontsize = 18)
                # layout of the axes
                ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1e'))
                ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1e'))
                                 
                # increment the plot image counter
                j += 1

        return allfig, allax

    #
    # plot jobs
    #

    def plotImages(self, plotSelect = None, plotType = 'norm'):
        """Select the plotting of the images

        plotSelect : selction of the images for plotting, use object default if None
        plotType   : plot raw images ('raw'), dark image substracted ('dark'),
                     normalized by ring current ('norm'), background substracted ('back')

        retrurns
        allfig : list of plt.figure objects of the plotting windows
        allax  : list of plt.axes objects which carry the figures"""

        # title of the full windows
        self._plotTitle  = '%s%s' % (self._imProc.setName, self._imProc.setNum)

        # read in first image to get conRoi if not set
        if self._imProc.conRoi == None:
            self._imProc._readImage(0)
        self._plotExtent = [self._imProc.conRoi[0], self._imProc.conRoi[0] + self._imProc.conRoi[1],
                            self._imProc.conRoi[2] + self._imProc.conRoi[3], self._imProc.conRoi[2]]

        # selection of images by passed, default or set size
        if plotSelect == None:
            plotSelect = self._plotSelect
        if plotSelect == None:
            plotSelect = range(self._imProc.setSize)
        self._totImNum = len(plotSelect)
            
        # image No. as title of each image
        self._pTitles = ['test']*self._totImNum
        for i in range(self._totImNum):
            self._pTitles[i] = 'image # %d' % (plotSelect[i]) 

        # considered images
        if plotType == 'norm':
            self._images  = self._imProc.fileProcessor.getImage()[plotSelect]

        allfig, allax = self._plot2DData()

        return allfig, allax

class PlotGrid():
    """Plot Grid Class

    This class plots: 
    intensity, occupation of the grid parts (bins), histogram of the occupation
    of a given 2D or 1D grid"""

    def __init__(self, imProc):
        # image Processor to get all the needed data
        self.imProc = imProc
        # set proper labels for axes
        if self.imProc.frameMode != 4:
            self.setAxesLabels([ur"Qx (\u00c5$^{-1}$)", ur"Qy (\u00c5$^{-1}$)", ur"Qz (\u00c5$^{-1}$)"])
        else:
            self.setAxesLabels(['H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)'])
        # plot flags: intensity, missing grid parts, histogram of occupation of the grid parts
        self.plotFlag2D = 7
        self.plotFlag1D = 7
        self.logFlag1D  = 0
        self.logFlag2D  = 0
        self.histBin    = 50
        # Figure Size
        self._defaultFigureSize = (11, 8.5)
        # plot the 1D fits
        self.plot1DFit  = False
        

    #
    # set part
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

    def setHistBin(self, histBin = 50):
        """Set the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        self.histBin = histBin

    def getHistBin(self, histBin = 50):
        """Get the no. of bins for the histograms

        hisBin : no. of bins for the histograms of the occupation numbers"""
        
        return self.histBin

    def setPlot1DAxes(self, valSetAx0, labelsAx0):
        """Sets the axes for the 1D plots"""
        self.valSetAx0p1D = valSetAx0
        self.labelsAx0p1D = labelsAx0

    def setPlot1DData(self, intentLine, occuLine, plotTitle = ''):
        """Sets the intenstiy and occupation of each grid part (bin) for the 1D plots"""
        self.intentLine = intentLine
        self.occuLine   = occuLine
        self.title1D    = plotTitle
        
    def setPlot2DAxes(self, minSetAx1, maxSetAx1, minSetAx0, maxSetAx0, labelsAx1, labelsAx0):
        """Sets the axes for the 2D plots"""
        self.minSetAx1p2D = minSetAx1
        self.maxSetAx1p2D = maxSetAx1
        self.minSetAx0p2D = minSetAx0
        self.maxSetAx0p2D = maxSetAx0
        self.labelsAx1p2D = labelsAx1
        self.labelsAx0p2D = labelsAx0

    def setPlot2DData(self, intentArea, occuArea, plotTitle = ''):
        """Sets the intenstiy and occupation of each grid part (bin) for the 2D plots"""
        self.intentArea = intentArea
        self.occuArea   = occuArea
        self.title2D    = plotTitle
        
    #
    # plot settings
    #
 
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
        
    def setHistBin(self, histBin = 50):
        """Sets the number of bins for the histograms"""
        self.histBin

    def setPlot1DFit(self, plot1DFit):
        """Set plotting the 1D fits"""

        self.plot1DFit = plot1DFit
              
    def getPlot1DFit(self):
        """Get plotting the 1D fits"""

        return self.plot1DFit

    #
    # plot functions
    #

    def plot1DData(self):
        """Plots 1D Lines: intensity, occupation, histogram"""

        # aliases
        valSet     = self.valSetAx0p1D
        labelAx0   = self.labelsAx0p1D
        intentSet  = self.intentLine
        occuSet    = self.occuLine
        
        # figure for 1D plots        
        fig = plt.figure(figsize = self._defaultFigureSize)
        fig.suptitle(self.title1D, fontsize = 24)
        allax = []
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag1D & 1: plotHor += 1
        if self.plotFlag1D & 2: plotHor += 1
        if self.plotFlag1D & 4: plotHor += 1
        plotOffSet = 0
        
        # intensity
        if self.plotFlag1D & 1:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag1D & 1 :
                    ax.semilogy(valSet[i], intentSet[i], '-bo')
                else:
                    ax.plot(valSet[i], intentSet[i], '-bo')
                ax.set_xlabel(labelAx0[i], fontsize = 18)
                ax.set_ylabel('Intensity', fontsize = 18)
                if i == 0:
                    ax.set_title('Intensity of Data', fontsize = 20)
                allax.append(ax)

        # add 1D fits
        if self.plot1DFit == True:
            yFit, fitRes = self.imProc.get1DFit(valSet, intentSet, fitType = None, infoDes = self.title1D)
            for i in range(3):
                allax[i].plot(valSet[i], yFit[i], '-r')
       
        # occupation of the grid parts (bins)
        if self.plotFlag1D & 2:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag1D & 2 :
                    ax.semilogy(valSet[i], occuSet[i], '-bo')
                else:
                    ax.plot(valSet[i], occuSet[i], '-bo')
                ax.set_xlabel(labelAx0[i], fontsize = 18)
                ax.set_ylabel('Occupation No.', fontsize = 18)
                if i == 0:
                    ax.set_title('Occupation of the Bins', fontsize = 20)
                allax.append(ax)

        # histogram for occupation numbers
        if self.plotFlag1D & 4:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag1D & 4 :
                    ax.hist(occuSet[i], self.histBin, log=True, facecolor='green')
                else:
                    ax.hist(occuSet[i], self.histBin, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    ax.set_title('Histogram', fontsize = 20)
                allax.append(ax)

        return fig, allax

    def plot2DData(self):
        """Plots 2D Areas: intensity, occupation, histogram"""

        # aliases
        minSetAx1  = self.minSetAx1p2D
        maxSetAx1  = self.maxSetAx1p2D
        minSetAx0  = self.minSetAx0p2D
        maxSetAx0  = self.maxSetAx0p2D
        labelsAx1  = self.labelsAx1p2D
        labelsAx0  = self.labelsAx0p2D
        intentArea = self.intentArea
        occuArea   = self.occuArea 
        
        # figure for 2D plots        
        fig = plt.figure(figsize = self._defaultFigureSize)
        fig.suptitle(self.title2D, fontsize = 24)
        allax = []
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag2D & 1: plotHor += 1
        if self.plotFlag2D & 2: plotHor += 1
        if self.plotFlag2D & 4: plotHor += 1
        plotOffSet = 0
        
        # intensity
        if self.plotFlag2D & 1:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag2D & 1 :
                    minPosVal = (intentArea[i] > 0.0).min()
                    negPart   = intentArea[i] <= 0
                    negPart   = np.ones(negPart.shape) * minPosVal
                    cax  = ax.imshow(intentArea[i], norm=LogNorm(), extent = [minSetAx1[i], maxSetAx1[i], maxSetAx0[i], minSetAx0[i]])
                else:
                    cax  = ax.imshow(intentArea[i], extent = [minSetAx1[i], maxSetAx1[i], maxSetAx0[i], minSetAx0[i]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(labelsAx1[i], fontsize = 18)
                ax.set_ylabel(labelsAx0[i], fontsize = 18)
                if i == 0:
                    ax.set_title('Intensity of Data', fontsize = 20)
                allax.append(ax)
        
        # occupation of the grid parts (bins)
        if self.plotFlag2D & 2:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag2D & 2:
                    cax  = ax.imshow(occuArea[i], norm=LogNorm(), extent = [minSetAx1[i], maxSetAx1[i], maxSetAx0[i], minSetAx0[i]])
                else:
                    cax  = ax.imshow(occuArea[i], extent = [minSetAx1[i], maxSetAx1[i], maxSetAx0[i], minSetAx0[i]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(labelsAx1[i], fontsize = 18)
                ax.set_ylabel(labelsAx0[i], fontsize = 18)
                if i == 0:
                    ax.set_title('Occupation of the Bins', fontsize = 20)
                allax.append(ax)
                
        # histogram for occupation numbers
        if self.plotFlag2D & 4:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if self.logFlag2D & 4 :
                    ax.hist(np.ravel(occuArea[i]), self.histBin, log=True, facecolor='green')
                else:
                    ax.hist(np.ravel(occuArea[i]), self.histBin, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    ax.set_title('Histogram of Flattend Array', fontsize = 20)
                allax.append(ax)

        return fig, allax

    #
    # plot jobs
    #

    def plotGrid1D(self, calcMode = 'sum', intenFits = None):
        """Select and plots the 1D Lines of the data grid

        calcMode  : select which calculated values are plotted, 'sum', 'cut', 'cutAv'
        intenFits : intensity of the 1D fits, do not consider if None

        retrurns
        fig1   : plt.figure object of the plotting window
        allax1 : list of plt.axes objects which carry the figures
        allRes : all results of the fits, [[a1, b1, cen1, width1, area1],...], [0, 0, 0, 0, 0] if unsuccessful fit"""

        # results for fit of 1D data, None if no fitting
        allRes = np.zeros((3,5))

        # axes and data configuration
        self.setPlot1DAxes(self.imProc.qVal, self.axesLabels)
        if calcMode == 'sum':
            gridData1D, gridOccu1D = self.imProc.get1DSum()
            plotTitle = '1D Lines, over other directions is summed'
        elif calcMode == 'cut':
            gridData1D, gridOccu1D = self.imProc.get1DCut()
            plotTitle = '1D Line cuts at maximum position'
        else:
            gridData1D, gridOccu1D = self.imProc.get1DCutAv()
            plotTitle = '1D average over 9 line cuts around maximum position'
        self.setPlot1DData(gridData1D, gridOccu1D, plotTitle = plotTitle)
        # plot, get figure and axes back
        fig1, allax1 = self.plot1DData()
        # if there are fits, show them
        if intenFits != None:
            for i in range(3):
                if intentFits != None:
                    allax1[i].plot(self.imProc.qVal[i], intenFits[i], '-r')
       
        return fig1, allax1, allRes

    def plotGrid2D(self, calcMode = 'sum'):
        """Select and plots the 2D Areas of the data grid

        calcMode : select which calculated values are plotted, 'sum', 'cut', 'cutAv'

        retrurns
        fig2   : plt.figure object of the plotting window
        allax2 : list of plt.axes objects which carry the figures"""

        if self.imProc.frameMode != 4:
            axLabels =[ur"Qx (\u00c5$^{-1}$)", ur"Qy (\u00c5$^{-1}$)", ur"Qz (\u00c5$^{-1}$)"]
        else:
            axLabels = ['H (r.l.u.)', 'K (r.l.u.)', 'L (r.l.u.)']
        # axes and data configuration
        self.setPlot2DAxes([self.imProc.Qmin[2], self.imProc.Qmin[2], self.imProc.Qmin[1]], 
                           [self.imProc.Qmax[2], self.imProc.Qmax[2], self.imProc.Qmax[1]],
                           [self.imProc.Qmin[1], self.imProc.Qmin[0], self.imProc.Qmin[0]], 
                           [self.imProc.Qmax[1], self.imProc.Qmax[0], self.imProc.Qmax[0]],
                           [self.axesLabels[2], self.axesLabels[2], self.axesLabels[1]],
                           [self.axesLabels[1], self.axesLabels[0], self.axesLabels[0]])
        if calcMode == 'sum':
            gridData2D, gridOccu2D = self.imProc.get2DSum()
            plotTitle = '2D Areas, over other direction is summed'
        elif calcMode == 'cut':
            gridData2D, gridOccu2D = self.imProc.get2DCut()
            plotTitle = '2D Line cuts at maximum position'
        else:
            gridData2D, gridOccu2D = self.imProc.get2DCutAv()
            plotTitle = '2D average over 3 area cuts around maximum position'
        for i in range(3):
            gridData2D[i] = np.ma.array(gridData2D[i], mask = (gridOccu2D[i] == 0))
            gridOccu2D[i] = np.ma.array(gridOccu2D[i], mask = (gridOccu2D[i] == 0))
        self.setPlot2DData(gridData2D, gridOccu2D, plotTitle = plotTitle)
        # plot, get figure and axes back
        fig2, allax2 = self.plot2DData()

        return fig2, allax2

    def plotAll(self):
        """Plots 1D/2D sums and cuts"""
        self.plotGrid1D('sum')
        self.plotGrid2D('sum')
        self.plotGrid1D('cut')
        self.plotGrid2D('cut')
        self.plotGrid1D('cutAv')
        self.plotGrid2D('cutAv')
    

####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    from pyspec import spec
    from pyspec.ccd.transformations import FileProcessor, ImageProcessor

    sf   = spec.SpecDataFile('/home/tardis/spartzsch/data/2010_09_X1A2/ymn2o5_sep10_1', 
			     ccdpath = '/mounts/davros/nasshare/images/sept10')
    scan = sf[244]

    fp = FileProcessor()
    fp.setFromSpec(scan)
    fp.process()

    testData = ImageProcessor()

    testData.setDetectorAngle(-1.24)
    testData.setBins(4, 4)
    testData.setFileProcessor(fp)
    testData.setSpecScan(scan)
    #testData.setConRoi([1, 325, 1, 335])
    testData.setFrameMode(1)

    # grid data
    testData.setGridOptions(Qmin = None, Qmax = None, dQN = [90, 160, 30])
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [200, 400, 100])
    #testData.setGridOptions(Qmin = None, Qmax = None, dQN = [100, 100, 100])
    #testData.makeGridData()

    #testPlotter = PlotGrid3D(testData)
    #testPlotter.plot3D()

    # plotter for grid
    testPlotter = PlotGrid(testData)

    testPlotter.setLogFlags(7, 7)
    testPlotter.setPlot1DFit(True)
    #testPlotter.plotGrid1D('sum')
    #testPlotter.plotGrid1D('cut')
    #testPlotter.plotGrid1D('cutAv')
    #testPlotter.plotGrid2D('sum')
    #testPlotter.plotGrid2D('cut')
    #testPlotter.plotGrid2D('cutAv')
    #testPlotter.plotAll()

    """
    testData = np.array([[range(0,3),range(3,6)],[range(6,9),range(9,12)]])
    testOccu = np.array([[np.zeros(3),[2,5,7]],[[0,4,2],[7,0,8]]])
    testPlot = PlotGrid()
    testPlot.setPlotFlags(7, 7)
    testPlot.setLogFlags(0, 3)
    # axes configuration
    testPlot.setPlot1DAxes([range(2), range(2), range(3)], ['X', 'Y', 'Z'])
    testPlot.setPlot2DAxes([0, 0, 0], [2, 2, 1], [0, 0, 0], [1, 1, 1], ['Y', 'X', 'X'], ['Z', 'Z', 'Y'])
    # data for sums
    testPlot.setPlot1DData([testData.sum(1).sum(1), testData.sum(0).sum(1), testData.sum(0).sum(0)],
                           [testOccu.sum(1).sum(1), testOccu.sum(0).sum(1), testOccu.sum(0).sum(0)],
                           plotTitle = '1D Lines, over other directions is summed')
    testPlot.setPlot2DData([testData.sum(0), testData.sum(1), testData.sum(2)],
                           [testOccu.sum(0), testOccu.sum(1), testOccu.sum(2)],
                           plotTitle = '2D Areas, over other direction is summed')
    # plot data for sums
    fig1, allax1 = testPlot.plot1DData()
    fig2, allax2 = testPlot.plot2DData()
    """

    # plotter for images
    testPlotIm = PlotImages(testData)
    # cofigurations for image plotting
    testPlotIm.setPlotSelect(plotSelect = None)
    testPlotIm.setPlotLayout(figSize = (11, 8.5), plotHor = 3, plotVer = 2, plotOrd = 'vh')
    testPlotIm.setPlotFlag(3)
    testPlotIm.setLogFlag(3)
    testPlotIm.setHistBin(histBin = 100)
    # plot the images
    testPlotIm.plotImages(plotSelect = range(0, 81, 20))
    

    plt.show()

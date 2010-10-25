#
# plotter.py (c) Stuart B. Wilkins 2010 and (c) Sven Partzsch
#
# $Id: transformations.py 40 2010-10-13 15:58:30Z stuwilkins $
# $HeadURL: https://pyspec.svn.sourceforge.net/svnroot/pyspec/trunk/pyspec/ccd/transformations.py $
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

class PlotGrid():
    """Plot Grid Class

    This class plots a given 1D/2D/3D grid
    provides cuts and projections to lower the dimension"""

    def __init__(self):
        self.axesLabels = ['Qx (1/A)', 'Qy (1/A)', 'Qz (1/A)']
        # plot flags: intensity, missing grid parts, histogram of occupation of the grid parts
        self.plotFlag2D = 7
        self.plotFlag1D = 7
        self.normFlag2D = 2
        self.normFlag1D = 2
        self.logFlag2D  = 0
        self.fit1D      = False
        self.histBin    = 50

    #
    # set part
    #

    def setGrid(self, gridData, gridOccu):
        """Set the grid data set
        e.g. 2D grid: [[val00, val10,..., valN0],..., [val0M, val1M,.., valNM]]"""
        self.gridData = np.ma.array(gridData / gridOccu, mask = (gridOccu == 0))
        self.gridOccu = gridOccu
        self.gridSize = gridData.shape
        # use number of grid part (indicies) as temporary vectors
        self.minBox   = np.zeros(3)
        self.maxBox   = np.array(self.gridSize)
        self.dVec     = np.ones(3)
        self._calcVecDataSet()
        self._calcMax()

    def setGridVec(self, minBox, maxBox):
        """Sets the vectors for each point of the grid by the minimum and maximum values"""
        self.minBox   = minBox
        self.maxBox   = maxBox
        self.dVec     = (maxBox - minBox) / self.gridSize
        self._calcVecDataSet()

    def setAxesLabels(self, labels):
        """Set the labels for the axes"""
        self.axesLabels = labels

    def setPlotFlags(self, flag2D, flag1D):
        """Sets the ploting flags for 2D and 1D plots
        binary code, flag & 1: intensity, flag & 2: missing grid parts,
        flag & 4: histogram of occupation of the grid parts"""
        self.plotFlag2D = flag2D
        self.plotFlag1D = flag1D

    def setNormFlags(self, flag2D, flag1D):
        """Sets the normalization flags for 2D and 1D plots
        binary code, flag & 1: intensity, flag & 2: missing grid parts,
        flag & 4: histogram of occupation of the grid parts"""
        self.normFlag2D = flag2D
        self.normFlag1D = flag1D

    def setLogFlags(self, flag2D):
        """Sets whether 2D areas are plotted linear (False) of logarithmic (True)
        binary code, flag & 1: intensity, flag & 2: missing grid parts"""
        self.logFlag2D = flag2D
        
    def setFit1D(self, fit1D):
        """Sets whether 1D lines get fitted by a Loarenzian squared"""
        self.fit1D = fit1D

    def setHistBin(self, histBin):
        """Sets the number of bins for the histograms"""
        self.histBin
        
    #
    # help functions
    #

    def _calcVecDataSet(self):
        """Calculats the vector data set for the grid points"""

        # aliases
        minBox = self.minBox
        maxBox = self.maxBox
        dVec   = self.dVec
        
        # vector data set of the center of each grid part
        self.xVal = np.arange(minBox[0], maxBox[0] - dVec[0]/2, dVec[0]) + dVec[0]/2
        self.yVal = np.arange(minBox[1], maxBox[1] - dVec[1]/2, dVec[1]) + dVec[1]/2
        self.zVal = np.arange(minBox[2], maxBox[2] - dVec[2]/2, dVec[2]) + dVec[2]/2

    def _calcMax(self):
        """Calculates the position of the maximum as indicies"""
        maxN = self.gridData.argmax()
        ind2 = maxN % self.gridSize[2]
        ind1 = maxN / self.gridSize[2] % self.gridSize[1]
        ind0 = maxN / self.gridSize[2] / self.gridSize[1]
        self.maxInd = np.array([ind0, ind1, ind2])

    #
    # plot functions
    #
    
    def plot2DSum(self):
        """Plots 2D Areas, over the third dimension is summed"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        minBox     = self.minBox
        maxBox     = self.maxBox
        
        # figure for 2D plots        
        fig = plt.figure()
        fig.suptitle('2D Areas, summed in the other direction', fontsize = 24)
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag2D & 1: plotHor += 1
        if self.plotFlag2D & 2: plotHor += 1
        if self.plotFlag2D & 4: plotHor += 1
        plotOffSet = 0
        
        axesInd = [[2, 1], [2, 0], [1, 0]]
        sumInd  = [0, 1, 2]
        # intensity
        if self.plotFlag2D & 1:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                Area = intentSet.sum(sumInd[i])
                if self.normFlag2D & 1:
                    Area = Area / Area.max()
                if self.logFlag2D & 1 :
                    cax  = ax.imshow(Area, norm=LogNorm(), extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                else:
                    cax  = ax.imshow(Area, extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(axesLabels[axesInd[i][0]], fontsize = 18)
                ax.set_ylabel(axesLabels[axesInd[i][1]], fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 1:
                        ax.set_title('Rel. Intensity of Data', fontsize = 20)
                    else:
                        ax.set_title('Intensity of Data', fontsize = 20)
        
        # missing grid parts, relative
        if self.plotFlag2D & 2:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                Area = intentSet.mask.sum(sumInd[i])
                if self.normFlag2D & 2:
                    Area = Area / float(sizeBox[sumInd[i]])
                if self.logFlag2D & 2:
                    cax  = ax.imshow(Area, norm=LogNorm(), extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                else:
                    cax  = ax.imshow(Area, extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(axesLabels[axesInd[i][0]], fontsize = 18)
                ax.set_ylabel(axesLabels[axesInd[i][1]], fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 2:
                        ax.set_title('Rel. No. of Missing Data', fontsize = 20)
                    else:
                        ax.set_title('No. of Missing Data', fontsize = 20)
                    
        # histogram for occupation numbers
        if self.plotFlag2D & 4:
            plotOffSet += 1
            if self.normFlag2D & 4:
                normHis = 1
            else:
                normHis = 0
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                ax.hist(self.gridOccu.sum(sumInd[i]), self.histBin, normed=normHis, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 4:
                        ax.set_title('Rel. Histogram of Flattend Array', fontsize = 20)
                    else:
                        ax.set_title('Histogram of Flattend Array', fontsize = 20)
                    

    def plot1DSum(self):
        """Plots 1D Lines over the other two dimensions is summed"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        xVal       = self.xVal
        yVal       = self.yVal
        zVal       = self.zVal
        valSet     = (xVal, yVal, zVal)

        # figure for 1D plots        
        fig = plt.figure()
        fig.suptitle('1D Lines, summed in all other directions', fontsize = 24)
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag1D & 1: plotHor += 1
        if self.plotFlag1D & 2: plotHor += 1
        if self.plotFlag1D & 4: plotHor += 1
        plotOffSet = 0
        
        peakArea  = np.zeros(3)
        peakWidth = np.zeros(3)

        axesInd = [0, 1, 2]
        sumInd  = [[1, 1], [0, 1], [0, 0]]

        # intensity
        if self.plotFlag1D & 1:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                intentVal = intentSet.sum(sumInd[i][0]).sum(sumInd[i][1])
                if self.normFlag1D & 1:
                    intentVal = intentVal / intentVal.max()
                ax.plot(valSet[i], intentVal, '-bo')
                ax.set_xlabel(axesLabels[i], fontsize = 18)
                ax.set_ylabel('Intensity', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 1:
                        ax.set_title('Rel. Intensity of Data', fontsize = 20)
                    else:
                        ax.set_title('Intensity of Data', fontsize = 20)
                if self.fit1D == True:
                    f = fit.fitdata(funcs=[fitfuncs.linear, fitfuncs.lor2a])
                    peakWidth[i] = f.result[3]
                    peakArea[i]  = f.result[4]
       
        # missing grid parts, relative
        if self.plotFlag1D & 2:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                maskVal = intentSet.mask.sum(sumInd[i][0]).sum(sumInd[i][1])
                if self.normFlag1D & 2:
                    maskVal = maskVal / float(sizeBox[sumInd[i][0]])/sizeBox[sumInd[i][1]+1]
                ax.plot(valSet[i], maskVal, '-bo')
                ax.set_xlabel(axesLabels[i], fontsize = 18)
                ax.set_ylabel('Masked', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 2:
                        ax.set_title('Rel. No. of Missing Data', fontsize = 20)
                    else:
                        ax.set_title('No. of Missing Data', fontsize = 20)

        # histogram for occupation numbers
        if self.plotFlag1D & 4:
            plotOffSet += 1
            if self.normFlag1D & 4:
                normHis = 1
            else:
                normHis = 0
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                ax.hist(self.gridOccu.sum(sumInd[i][0]).sum(sumInd[i][1]), self.histBin, normed=normHis, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 4:
                        ax.set_title('Rel. Histogram', fontsize = 20)
                    else:
                        ax.set_title('Histogram', fontsize = 20)

    def plot2DCut(self):
        """Plots 2D Areas, cuts at maximum intensity"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        minBox     = self.minBox
        maxBox     = self.maxBox
        maxInd     = self.maxInd
        
        # figure for 2D plots        
        fig = plt.figure()
        fig.suptitle('2D Areas, cuts at maximum intensity', fontsize = 24)
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag2D & 1: plotHor += 1
        if self.plotFlag2D & 2: plotHor += 1
        if self.plotFlag2D & 4: plotHor += 1
        plotOffSet = 0
        
        axesInd = [[2, 1], [2, 0], [1, 0]]
        # intensity
        if self.plotFlag2D & 1:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    Area = intentSet[maxInd[0],:,:]
                elif i == 1:
                    Area = intentSet[:,maxInd[1],:]
                else:
                    Area = intentSet[:,:,maxInd[2]]
                if self.normFlag2D & 1:
                    Area = Area / Area.max()
                if self.logFlag2D & 1 :
                    cax  = ax.imshow(Area, norm=LogNorm(), extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                else:
                    cax  = ax.imshow(Area, extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(axesLabels[axesInd[i][0]], fontsize = 18)
                ax.set_ylabel(axesLabels[axesInd[i][1]], fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 1:
                        ax.set_title('Rel. Intensity of Data', fontsize = 20)
                    else:
                        ax.set_title('Intensity of Data', fontsize = 20)
        
        # missing grid parts, relative
        if self.plotFlag2D & 2:
            plotOffSet += 1
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    Area = intentSet.mask[maxInd[0],:,:]
                elif i == 1:
                    Area = intentSet.mask[:,maxInd[1],:]
                else:
                    Area = intentSet.mask[:,:,maxInd[2]]
                if self.normFlag2D & 2:
                    Area = Area / Area.max()
                if self.logFlag2D & 2:
                    cax  = ax.imshow(Area, norm=LogNorm(), extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                else:
                    cax  = ax.imshow(Area, extent = [minBox[axesInd[i][0]], maxBox[axesInd[i][0]], minBox[axesInd[i][1]], maxBox[axesInd[i][1]]])
                fig.colorbar(cax)
                ax.set_aspect(1./ax.get_data_ratio())
                ax.set_xlabel(axesLabels[axesInd[i][0]], fontsize = 18)
                ax.set_ylabel(axesLabels[axesInd[i][1]], fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 2:
                        ax.set_title('Rel. No. of Missing Data', fontsize = 20)
                    else:
                        ax.set_title('No. of Missing Data', fontsize = 20)
                    
        # histogram for occupation numbers
        if self.plotFlag2D & 4:
            plotOffSet += 1
            if self.normFlag2D & 4:
                normHis = 1
            else:
                normHis = 0
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    hisData = self.gridOccu[maxInd[0],:,:]
                elif i == 1:
                    hsiData = self.gridOccu[:,maxInd[1],:]
                else:
                    hisData = self.gridOccu[:,:,maxInd[2]]
                ax.hist(hisData, self.histBin, normed=normHis, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    if self.normFlag2D & 4:
                        ax.set_title('Rel. Histogram of Flattend Array', fontsize = 20)
                    else:
                        ax.set_title('Histogram of Flattend Array', fontsize = 20)
                    

    def plot1DCut(self):
        """Plots 1D Lines, cut at position of maximum intensity"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        xVal       = self.xVal
        yVal       = self.yVal
        zVal       = self.zVal
        valSet     = (xVal, yVal, zVal)
        maxInd     = self.maxInd

        # figure for 1D plots        
        fig = plt.figure()
        fig.suptitle('1D Lines, cut at position of maximum intensity', fontsize = 24)
        # how many horizontal plots?
        plotHor = 0
        if self.plotFlag1D & 1: plotHor += 1
        if self.plotFlag1D & 2: plotHor += 1
        if self.plotFlag1D & 4: plotHor += 1
        plotOffSet = 0
        
        peakArea  = np.zeros(3)
        peakWidth = np.zeros(3)

        axesInd = [0, 1, 2]

        # intensity
        if self.plotFlag1D & 1:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    intentVal = intentSet[:,maxInd[1],maxInd[2]]
                elif i == 1:
                    intentVal = intentSet[maxInd[0],:,maxInd[2]]
                else:
                    intentVal = intentSet[maxInd[0],maxInd[1],:]
                if self.normFlag1D & 1:
                    intentVal = intentVal / intentVal.max()
                ax.plot(valSet[i], intentVal, '-bo')
                ax.set_xlabel(axesLabels[i], fontsize = 18)
                ax.set_ylabel('Intensity', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 1:
                        ax.set_title('Rel. Intensity of Data', fontsize = 20)
                    else:
                        ax.set_title('Intensity of Data', fontsize = 20)
                if self.fit1D == True:
                    f = fit.fitdata(funcs=[fitfuncs.linear, fitfuncs.lor2a])
                    peakWidth[i] = f.result[3]
                    peakArea[i]  = f.result[4]
       
        # missing grid parts, relative
        if self.plotFlag1D & 2:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    maskVal = intentSet.mask[:,maxInd[1],maxInd[2]]
                elif i == 1:
                    maskVal = intentSet.mask[maxInd[0],:,maxInd[2]]
                else:
                    maskVal = intentSet.mask[maxInd[0],maxInd[1],:]
                if self.normFlag1D & 2:
                    maskVal = maskVal / maskVal.max()
                ax.plot(valSet[i], maskVal, '-bo')
                ax.set_xlabel(axesLabels[i], fontsize = 18)
                ax.set_ylabel('Masked', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 2:
                        ax.set_title('Rel. No. of Missing Data', fontsize = 20)
                    else:
                        ax.set_title('No. of Missing Data', fontsize = 20)

        # histogram for occupation numbers
        if self.plotFlag1D & 4:
            plotOffSet += 1
            if self.normFlag1D & 4:
                normHis = 1
            else:
                normHis = 0
            for i in range(3):
                ax   = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
                if i == 0:
                    hisData = self.gridOccu[:,maxInd[1],maxInd[2]]
                elif i == 1:
                    hisData = self.gridOccu[maxInd[0],:,maxInd[2]]
                else:
                    hisData = self.gridOccu[maxInd[0],maxInd[1],:]
                ax.hist(hisData, self.histBin, normed=normHis, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    if self.normFlag1D & 4:
                        ax.set_title('Rel. Histogram', fontsize = 20)
                    else:
                        ax.set_title('Histogram', fontsize = 20)

    def plotAll(self):
        """Plots 2D/1D intensities, masked parts, and histograms of sums, and cuts"""
        self.setPlotFlags(7,7)
        self.plot2DSum()
        self.plot1DSum()
        self.plot2DCut()
        self.plot1DCut()

####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

    testData = np.array([[range(0,3),range(3,6)],[range(6,9),range(9,12)]])
    testOccu = np.array([[np.zeros(3),[2,5,7]],[[0,4,2],[7,0,8]]])
    testPlot = PlotGrid()
    testPlot.setGrid(testData, testOccu)
    #testPlot.setPlotFlags(7, 7)
    testPlot.setNormFlags(2, 2)
    testPlot.setLogFlags(0)
    #testPlot.plot2DSum()
    #testPlot.plot1DSum()
    #testPlot.plot2DCut()
    testPlot.plot1DCut()
    #testPlot.plotAll()
    plt.show()

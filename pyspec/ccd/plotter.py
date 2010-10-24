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
from   pyspec  import fit, fitfuncs

class PlotGrid():
    """Plot Grid Class

    This class plots a given 1D/2D/3D grid
    provides cuts and projections to lower the dimension"""

    def __init__(self):
        self.axesLabels = ['Qx', 'Qy', 'Qz']

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

    def setGridVec(self, minBox, maxBox):
        """Sets the vectors for each point of the grid by the minimum and maximum values"""
        self.minBox   = minBox
        self.maxBox   = maxBox
        self.dVec     = (maxBox - minBox) / self.girdSize

    def setAxesLabels(self, labels):
        """Set the labels for the axes"""
        self.axesLabels = labels

    def plot2DPro(self):
        """Plots the 2D projections of the data set
        over the third dimension is summed"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        dVec       = self.dVec
        minBox     = self.minBox
        maxBox     = self.maxBox
        
        # vector data set
        xVal = np.arange(minBox[0], maxBox[0] - dVec[0]/2, dVec[0]) + dVec[0]/2
        yVal = np.arange(minBox[1], maxBox[1] - dVec[1]/2, dVec[1]) + dVec[1]/2
        zVal = np.arange(minBox[2], maxBox[2] - dVec[2]/2, dVec[2]) + dVec[2]/2
        
        fig1 = plt.figure()
        plotHor = 2
        
        # yz Area
        ax = fig1.add_subplot(3, plotHor, 1)
        yzArea = intentSet.sum(0)
        cax = ax.imshow(yzArea)#, extent = [minBox[2], maxBox[2], minBox[1], maxBox[1]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Qy', fontsize = 18)
        ax.set_title('Data', fontsize = 24)
        # xz Area
        ax = fig1.add_subplot(3, plotHor, plotHor+1)
        xzArea = intentSet.sum(1)
        cax = ax.imshow(xzArea)#, extent = [minBox[2], maxBox[2], minBox[0], maxBox[0]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Qx', fontsize = 18)
        # xy Area
        ax = fig1.add_subplot(3, plotHor, 2*plotHor+1)
        xyArea = intentSet.sum(2)
        cax = ax.imshow(xyArea)#, extent = [minBox[1], maxBox[1], minBox[0], maxBox[0]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qy', fontsize = 18)
        ax.set_ylabel('Qx', fontsize = 18)
            
        # yz Area - masked
        ax  = fig1.add_subplot(3, plotHor, 2)
        cax = ax.imshow(intentSet.mask.sum(0)/float(sizeBox[0]))#, extent = [minBox[2], maxBox[2], minBox[1], maxBox[1]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Qy', fontsize = 18)
        ax.set_title('Relativ No. of Missing Data', fontsize = 24)
        # xz Area - masked
        ax  = fig1.add_subplot(3, plotHor, plotHor+2)
        cax = ax.imshow(intentSet.mask.sum(1)/float(sizeBox[1]))#, extent = [minBox[2], maxBox[2], minBox[0], maxBox[0]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Qx', fontsize = 18)
        # xy Area - masked
        ax  = fig1.add_subplot(3, plotHor, 2*plotHor+2)
        cax = ax.imshow(intentSet.mask.sum(2)/float(sizeBox[2]))#, extent = [minBox[1], maxBox[1], minBox[0], maxBox[0]])
        fig1.colorbar(cax)
        ax.set_aspect(1./ax.get_data_ratio())
        ax.set_xlabel('Qy', fontsize = 18)
        ax.set_ylabel('Qx', fontsize = 18)

    def plot1DPro(self):
        """Plots the 1D projections of the data set
        over the other two dimensions is summed"""

        # aliases
        axesLabels = self.axesLabels
        intentSet  = self.gridData
        sizeBox    = self.gridSize
        dVec       = self.dVec
        minBox     = self.minBox
        maxBox     = self.maxBox
        
        # vector data set
        xVal = np.arange(minBox[0], maxBox[0] - dVec[0]/2, dVec[0]) + dVec[0]/2
        yVal = np.arange(minBox[1], maxBox[1] - dVec[1]/2, dVec[1]) + dVec[1]/2
        zVal = np.arange(minBox[2], maxBox[2] - dVec[2]/2, dVec[2]) + dVec[2]/2
        
        fig2 = plt.figure()
        plotHor = 2

        peakArea  = np.zeros(3)
        peakWidth = np.zeros(3)
        # x Line
        ax = fig2.add_subplot(3, plotHor, 1)
        ax.plot(xVal, intentSet.sum(1).sum(1), '-bo')
        ax.set_xlabel('Qx', fontsize = 18)
        ax.set_ylabel('Intensity', fontsize = 18)
        ax.set_title('Data', fontsize = 24)
        #f = fit.fitdata(funcs=[fitfuncs.linear, fitfuncs.lor2a])
        #peakWidth[0] = f.result[3]
        #peakArea[0]  = f.result[4]
        # y Line
        ax = fig2.add_subplot(3, plotHor, plotHor+1)
        ax.plot(yVal, intentSet.sum(0).sum(1), '-bo')
        ax.set_xlabel('Qy', fontsize = 18)
        ax.set_ylabel('Intensity', fontsize = 18)
        #f = fit.fitdata(funcs=[fitfuncs.linear, fitfuncs.lor2a])
        #peakWidth[1] = f.result[3]
        #peakArea[1]  = f.result[4]
        # z Line
        ax = fig2.add_subplot(3, plotHor, 2*plotHor+1)
        ax.plot(zVal, intentSet.sum(0).sum(0), '-bo')
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Intensity', fontsize = 18)
        #f = fit.fitdata(funcs=[fitfuncs.linear, fitfuncs.lor2a])
        #peakWidth[2] = f.result[3]
        #peakArea[2]  = f.result[4]

        # x Line - masked
        ax = fig2.add_subplot(3, plotHor, 2)
        ax.plot(xVal, intentSet.mask.sum(1).sum(1)/float(sizeBox[1])/sizeBox[2], '-bo')
        ax.set_xlabel('Qx', fontsize = 18)
        ax.set_ylabel('Masked', fontsize = 18)
        ax.set_title('Relativ No. of Missing Data', fontsize = 24)
        # y Line - masked
        ax = fig2.add_subplot(3, plotHor, plotHor+2)
        ax.plot(yVal, intentSet.mask.sum(0).sum(1)/float(sizeBox[0])/sizeBox[2], '-bo')
        ax.set_xlabel('Qy', fontsize = 18)
        ax.set_ylabel('Masked', fontsize = 18)
        # z Line - masked
        ax = fig2.add_subplot(3, plotHor, 2*plotHor+2)
        ax.plot(zVal, intentSet.mask.sum(0).sum(0)/float(sizeBox[0])/sizeBox[1], '-bo')
        ax.set_xlabel('Qz', fontsize = 18)
        ax.set_ylabel('Masked', fontsize = 18)

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
    testPlot.plot2DPro()
    testPlot.plot1DPro()
    plt.show()

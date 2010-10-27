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

    This class plots: 
    intensity, occupation of the grid parts (bins), histogram of the occupation
    of a given 2D or 1D grid"""

    def __init__(self):
        # plot flags: intensity, missing grid parts, histogram of occupation of the grid parts
        self.plotFlag2D = 7
        self.plotFlag1D = 7
        self.logFlag1D  = 0
        self.logFlag2D  = 0
        self.histBin    = 50

        self._defaultFigSize = (11, 8.5)

    #
    # set part
    #
    
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
        fig = plt.figure(figsize = self._defaultFigSize)
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
                ax.plot(valSet[i], intentSet[i], '-bo')
                ax.set_xlabel(labelAx0[i], fontsize = 18)
                ax.set_ylabel('Intensity', fontsize = 18)
                if i == 0:
                    ax.set_title('Intensity of Data', fontsize = 20)
                allax.append(ax)
       
        # occupation of the grid parts (bins)
        if self.plotFlag1D & 2:
            plotOffSet += 1
            for i in range(3):
                ax = fig.add_subplot(3, plotHor, i*plotHor+plotOffSet)
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
        fig = plt.figure(figsize = self._defaultFigSize)
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
                    cax  = ax.imshow(intentArea[i], norm=LogNorm(), extent = [minSetAx1[i], maxSetAx1[i], minSetAx0[i], maxSetAx0[i]])
                else:
                    cax  = ax.imshow(intentArea[i], extent = [minSetAx1[i], maxSetAx1[i], minSetAx0[i], maxSetAx0[i]])
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
                    cax  = ax.imshow(occuArea[i], norm=LogNorm(), extent = [minSetAx1[i], maxSetAx1[i], minSetAx0[i], maxSetAx0[i]])
                else:
                    cax  = ax.imshow(occuArea[i], extent = [minSetAx1[i], maxSetAx1[i], minSetAx0[i], maxSetAx0[i]])
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
                ax.hist(np.ravel(occuArea[i]), self.histBin, facecolor='green')
                ax.set_xlabel('No. of Occupations', fontsize = 18)
                ax.set_ylabel('No. of Grid Parts', fontsize = 18)
                if i == 0:
                    ax.set_title('Histogram of Flattend Array', fontsize = 20)
                allax.append(ax)

        return fig, allax
    

####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

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
    plt.show()

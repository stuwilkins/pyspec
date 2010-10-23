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

class plotGrid():
    """Plot Grid Class

    This class plots a given 1D/2D/3D grid
    provides cuts and projections to lower the dimension"""

    def __init__(self):
        self.plotLim    = None
        self.axesLabels = None

    #
    # set part
    #

    def setGrid(self, gridData):
        """Set the grid data set
        e.g. 2D grid: [[val00, val10,..., valN0],..., [val0M, val1M,.., valNM]]"""
        self.gridData = gridData
        self.gridSize = gridData.shape

    def setGridVec(self, minBox, maxBox):
        """Sets the vectors for each point of the grid by the minimum and maximum values"""
        self.minBox   = minBox
        self.maxBox   = maxBox
        self.dVec     = (maxBox - minBox) / self.girdSize

    def setPlotLim(self, lim):
        self.plotLim = lim

    def setAxesLabels(self, labels):
        """Set the labels for the axes"""
        self.axesLabels = labels

    def plot2D(self, areaVal, fig = None, ax = None):
        """Plots the 2D data set"""
        if fig == None:
            fig = figure()
        if ax == None:
            ax  = fig.subplot(111)
        cax = ax.imshow(areaVal)
        fig.colorbar(cax)




####################################
#
# main program, for test
#
####################################


if __name__ == "__main__":

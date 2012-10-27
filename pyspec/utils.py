#
# utils.py (c) Stuart B. Wilkins 2008
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
"""Set of useful utilities for plotting, and life in general!"""
from __future__ import with_statement
import pylab as pl
import numpy as np
import pickle 
import pyspec
import os
from matplotlib.ticker import MultipleLocator, MaxNLocator, FormatStrFormatter

# Some constants
golden_mean = (np.sqrt(5)-1.0)/2.0

class MultiPlot():
    def __init__(self, shape, prefuncs = [], postfuncs = []):
        """Initialize the class

        shape : (tuple) Size of subplots
        """

        self._shape = shape
        self._currentAxis = 1
        self._figures = []
        
        if type(prefuncs) != list:
            self._preFuncs = [prefuncs]
        else:
            self._preFuncs = prefuncs
            
        if type(postfuncs) != list:
            self._postFuncs = [postfuncs]
        else:
            self._postFuncs = postfuncs
        return

    def saveAll(self, fileformat, *args, **kwargs):
        """Save all figures to files

        This routine saves all figures to files. The fileformat
        should be a valid filename with a standard python string
        with place for an integer for the figure number. For example:

             MultiPlot.saveAll('Plot%d.png')

        Any further arguments are passed directly to the savefig function.

        """
        for n,f in enumerate(self._figures):
            figure(f)
            savefig(fileformat % n, *args, **kwargs)
        
    def getAxis(self):
        """Return the next axis"""
        # First check if we have to create a new figure
        if self._figures == []:
            self._figures.append(pl.figure())
            
        if (self._currentAxis - 1) == (self._shape[0] * self._shape[1]):
            self._currentAxis = 1
            self._figures.append(pl.figure())
        ax = self._figures[-1].add_subplot(self.shape[0], self.shape[1], self.currentAxis)
        self._currentAxis = self._currentAxis + 1
        return ax

    def getFigure(self):
        """Return the handle to the current figure"""
        return self._figures[-1]

    def getAllFigures(self):
        """Return all figures as a list"""
        return self._figures

def makeDataSquare(grids, sz = 0, ax = 0, verbose = False):
    """Make the nearest square from size"""
  
    # Note SZ [0] is the inner size
  
    if sz == 0:
        # Guess the split using differential
        x = abs(diff(grids[ax]))
        t = where(x > (0.25 * x.max()))[0]
        sz = int(t[0] + 1)

    l = int(grids[0].size / sz)
  
    new_grids = []
    for g in grids:
        new_grids.append(g[0:sz * l].reshape(l, sz))
    return new_grids
    

def rebin(a, nsize, factor = True, norm = False):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows. If norm is ture then
    normalize by area of rebinning. If factor is true binning rate is used.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    if factor:
        factor = np.asarray(nsize)
        nsize = np.asarray(shape) / factor
    else:
        factor = np.asarray(shape)/np.asarray(nsize)
    
    evList = ['a.reshape('] + \
             ['nsize[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)]
    if norm:
        evList.append(['/factor[%d]'%i for i in range(lenShape)])
    print ''.join(evList)

    return eval(''.join(evList))


def printAll(fname = "Plot", hc = False, small = True, lpcommand = 'lp'):
    """Print all figures"""
    for x in range(1, pl.gcf().number + 1):
        f = pl.figure(x)
        if small:
            f.set_size_inches(4, 3)
        pl.savefig('%s%03d.ps' % (fname, x), papertype = 'letter', orientation = 'portrait')
        if hc:
            os.system('%s %s%03d.ps' % (lpcommand, fname, x))

def multifit(sf, scans, var, *args, **kwargs):
    alldata = np.array([])
    allerrors = np.array([])
    pl.figure(figsize=(8.5, 11))

    firstrun = True
    plotnum = 1

    npp = kwargs.get('npp', (6, 2));
    if kwargs.has_key('npp'):
        del kwargs['bgnd']
    bgnd = kwargs.get('bgnd', None);
    if kwargs.has_key('npp'):
        del kwargs['bgnd']
    guessfollow = kwargs.get('guessfollow', False);
    if kwargs.has_key('guessfollow'):
        del kwargs['guessfollow']

    if type(var) == tuple:
        xvar = var[0]
        pvar = var[1]
        if len(var) > 2:
            mvar = var[2]
        else:
            mvar = None

        if len(var) > 3:
            yvar = var[3]
        else:
            yvar = None

    else:
        xvar = var
        pvar = None
        mvar = None
        yvar = None

    for scan in scans:
        pl.subplot(npp[0], npp[1], plotnum)

        sf[scan].plot(new = False, notitles = True, xcol = pvar, ycol = yvar,mcol = mvar)

        pl.subplots_adjust(hspace=0.4)

        if firstrun == True:
            f = pyspec.fit.fitdata(*args, **kwargs)
            firstrun = False
        else:
            if guessfollow:
                kwargs['guess'] = f.result
            f = pyspec.fit.fitdata(*args, **kwargs)

        exec '_xvar = sf[scan].%s' % var[0] 

        alldata = np.concatenate((alldata, np.array([np.mean(_xvar)]), f.result))
        allerrors = np.concatenate((allerrors, np.array([np.std(_xvar)]), f.stdev))

        pl.title('[%s] %f +- %f' % (scan, np.mean(_xvar), np.std(_xvar)))

        if(plotnum == (npp[0] * npp[1])):
            pl.figure(figsize=(8.5, 11))
            plotnum = 1
        else:
            plotnum += 1
            
    alldata = alldata.reshape(-1, len(f.result)+1)
    allerrors = allerrors.reshape(-1, len(f.result)+1)

    return alldata, allerrors

##
## Pickle routines to easily pickle data
##

def pickleit(filename, *args):
    """Pickle a python object

    filename  : string
       Filename of file to pickle objects to.
    args      : list
       List of python objects to pickle

    This function can be used as a quick convenience function
    for pickling python objects. For example::

       >>> a = 'hello'
       >>> b = array([1, 2, 3, 4])
       >>> c = SpecDataFile('myfile.01')
       >>> pickleit('spam.pckl', a, b, c)

    """

    output = open(filename, 'wb')

    if len(args) == 1:
        if type(args[0]) != list:
            pobject = [args[0]]
        else:
            pobject = args[0]
    else:
        pobject = list(args)

    for o in pobject:
        pickle.dump(o, output)

    print "**** Pickled %d objects to %s" % (len(pobject), filename)
    output.close()

def unpickleit(filename):
    """Unpickle a python object (created with pickleit)
    
    filename : string
       filename of file to unpickle.

    This routine can be used to easily unpickle a file and serves
    as the companion to the :func:'pickle()' function.

    If only one python object is unpickled then that object is
    returned.

    If multiple objects are pickled, then a list of these objects
    is returned."""

    o = []
    f = open(filename, 'rb')
    while(1):
        try:
            o.append(pickle.load(f))
        except EOFError:
            f.close()
            if len(o) == 1:
                # Single value so return first item
                print "**** Unpickled 1 object from %s." % filename 
                return o[0]
            else:
                print "**** Unpickled %d objects from %s." % (len(o), filename)
                return o

    return None

##
## Plot Formatting utilities
##

def makePanelPlot(n = 3, fig = None, 
                  xlmargin = 0.15, ytmargin = 0.10,
                  xrmargin = 0.05, ybmargin = 0.10,
                  ylabels = True):
    """Make a multi panel plot from matplotlib.

    This function, makes a typical panel plot and returns a list
    of the axes objects for plotting. 

    n : int
       Number of panels
    fig : figure object
       Figure object to use (If None creates new figure)
    xmargin : float
       Margin at x-axis
    ymargin : float
       Margin at y-axis

    """
    
    if fig is None:
        fig = pl.figure(figsize = [6, 6 * golden_mean * n])
    
    xsize = (1. - (xlmargin + xrmargin)) 
    ysize = (1. - (ybmargin + ytmargin)) / n

    pos = np.array([xlmargin, ybmargin, xsize, ysize])

    allax = []
    for x in range(n):
        ax = fig.add_axes(pos + np.array([0, ysize * x, 0, 0]))
        if x > 0:
            # Remove ticklabels
            ax.xaxis.set_ticklabels("")
        allax.append(ax)

    return allax

def makeNicePlot(ax, xint = 5, yint = 5, mxint = 4, myint = 4,
                 tickcolor = None, framecolor = None,
                 xformat = None, yformat = None):
    
    """Make nice plot by setting all border widths etc. 
    

    Designed to make a pedantic figure suitable for PRL, PRB etc.
    This function formats the border, tickmarks and sets the 
    axes labels and titles to sensible values for 'production' 
    quality plots.

    xint    : int
       Number of x intervals for minor tics
    yint    : int
       Number of y intervals for minor tics
    xformat : string
       Format string for major labels
    yformat : string
       Format string for minor labels
    """
    
    for i in ax.spines.itervalues():
        i.set_linewidth(2)
    if framecolor:
        for i in ax.spines.itervalues():
            i.set_color(framecolor)

    # Set the x and y axes fontsizes
    ax.xaxis.label.set_fontsize(18)
    ax.yaxis.label.set_fontsize(18)

    # Set the title fontsize and place a little higher
    ax.title.set_fontsize(18)
    ax.title.set_position([0.5, 1.03])

    # Set number formatting
    if xformat:
        ax.xaxis.set_major_formatter(FormatStrFormatter(xformat))

    if yformat:
        ax.yaxis.set_major_formatter(FormatStrFormatter(yformat))

    if xint:
        ax.xaxis.set_major_locator(MaxNLocator(mxint))
        ax.xaxis.set_minor_locator(MaxNLocator(mxint * xint))

    if yint:
        ax.yaxis.set_major_locator(MaxNLocator(myint))
        ax.yaxis.set_minor_locator(MaxNLocator(myint * yint))

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.tick1line.set_markersize(10)
        tick.tick2line.set_markersize(10)
        tick.tick1line.set_markeredgewidth(2)
        tick.tick2line.set_markeredgewidth(2)

        if tickcolor:
            tick.tick1line.set_color(tickcolor)
            tick.tick2line.set_color(tickcolor)

    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(8)
        tick.tick2line.set_markersize(8)
        tick.tick1line.set_markeredgewidth(1)
        tick.tick2line.set_markeredgewidth(1)

        if tickcolor:
            tick.tick1line.set_color(tickcolor)
            tick.tick2line.set_color(tickcolor)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)


def setImageRange(data, limits, bins = 100):
    
    h,b = np.histogram(data, bins = bins)
    b = np.arange(data.min(), data.max(), 
                  (data.max() - data.min()) / bins)
    com = (h * np.arange(h.size)).sum() /  h.sum()
    limits = (np.array(limits) / 100.0) * bins

    if (com - limits[0]) < 0:
        dmin = data.min()
    else:
        dmin = b[int(com - limits[0])]
    
    if (com + limits[1]) >= bins:
        dmax = data.max()
    else:
        dmax = b[int(com + limits[1])]

    return dmin, dmax

def writeString(filename, string, append = True):
    """Write a sring out to a file.

    filename : string
       Filename to write to.
    string : string or object
       String to write to file. If this is a python object then
       :func:'str()' is called to make a string output.
    append : bool
       If true append text to file"""

    if append:
        fmt = 'a'
    else:
        fmt = 'w'
    f = open(filename, fmt)
    if type(string) == type('string'):
        f.write(string)
    else:
        f.write(str(string))

    f.close()


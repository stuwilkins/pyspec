"""
Python module for handling of spec data files (pyspec).

Written by S. Wilkins
(c) Stuart Wilkins 2007


SPEC is a Diffractometer control program from
Certified Scientific Software.
http://www.certif.com/


Example:

	# Load datafile and index
	
	sd = SpecDataFile('fourc.01')		
	
	# Load scan #1000 into scan
	
	scan = sd[1000]				
				
	# Plot scan #1000
	
	scan.plot()		
	
	# Print (Hardcopy) scan #1000
	
	scan.plot().prints()				
	
	# Access to variables
	
	scan.data							
	scan.values['TTH']					 
	scan.TTH							

	# Fit the current plot to a 
	# functions lor2a and linear
	scan.fit([lor2a, linear])			
										
										
"""

import time
import sys
import os
from scipy.optimize import leastsq
from numpy import *
from scipy import *
from pylab import *
import matplotlib.pylab as pylab
from fit import fit

__version__ = "$Revision$"
__author__ = "Stuart B. Wilkins <swilkins@bnl.gov>"
__date__ = "$LastChangedDate$"
__id__ = "$Id$"
__verbose__ = 0x01

NotInFile = "Scan is not in Datafile"
	
def removeIllegals(key):
	illegal = ['/', '-', ' ']
	for j in illegal:
		key = key.replace(j,'')
	if key[0].isalpha() == False:
		key = "X" + key
		
	return key
	
class SpecDataFile:
	"""
	DataFile class for handling spec data files
	"""
	def __init__(self,fn):
		self.filename = fn
		if __verbose__:
			print "**** Opening specfile %s." % self.filename
		self.index()
		self.readHeader()
		print self.getStats()
		
		self.scandata = {}
		
		return
	
	def readHeader(self):
		"""
		Read the spec header from the datafile by
		scanning the datafile.
		
		Currently supported header items:
			'#O'	(Motor positions)
		"""
		
		self.file = open(self.filename, 'r')
		
		if __verbose__:
			print "---- Reading Header."
		
		self.motors = []
		self.file.seek(0,0)
		line = self.file.readline()
		while line[0:2] != "#S":
			if line[0:2] == "#O":
				self.motors = self.motors + line.strip()[4:].split()
			line = self.file.readline()
			
		self.file.close()
		return

	def index(self):
		"""
		Index the datafile by storing the byte-offests for
		all the scans (Lines beginning with '#S')
		
		usage: SpecDataFile.index()
		"""
		self.file = open(self.filename, 'r')
		
		if __verbose__:
			print "---- Indexing scan :       ",
			sys.stdout.flush()
			sys.stderr.flush()
			
		self.file.seek(0,0)
		self.findex = {}
		
		pos = self.file.tell()
		line = self.file.readline()
		while line != "":
			if line[0:2] == "#S":
				a = line.split()
				s = int(a[1])
				if (s % 5) is 0:
					print "\b\b\b\b\b\b\b%5d " % s,
					sys.stdout.flush()
				self.findex[s] = pos
			pos = self.file.tell()
			line = self.file.readline()
		print "\b\b\b\b\b\b\bDONE  "
		
		self.file.close()
		return
	
	def getStats(self, head = "---- "):
		"""
		returns string with statistics on specfile
		
		usage SpecDataFile.showStats()
		"""
		
		string = ""
		string = string + head + "Specfile contains %d scans\n" % len(self.findex)
		string = string + head + "Start scan = %d\n" % min(self.findex.keys())
		string = string + head + "End   scan = %d\n" % max(self.findex.keys())
		
		return string
		
	def _moveto(self, item):
		"""
		Move to a location in the datafile for scan
		"""
		if self.findex.has_key(item):
			self.file.seek(self.findex[item])
		else:
			raise Exception("Scan %s is not in datafile ....." % item)
			
	def __getitem__( self, item):
		""" 
		Convinience routine to use [] to get scan
		"""
		return self.getScan(item, setkeys = True)
		
	def getScan(self, item, setkeys = True):
		"""
		Get a scan from the data file and load it into the
		list of SpecData instances. 
		
		setting 'setkeys = True' will set attributes to the
		specdata item of all the motors and counters
		
		returns the scandata instance
		"""
		if type(item) == int:
			items = [item]
		elif type(item) == list:
			items = item
		else:
			raise Exception("item can only be <int> or <list>")
		
		self.file = open(self.filename, 'r')
		rval = []
		n = 0
		for i in items:
			if self.scandata.has_key(i) is False:
				if __verbose__:
					print "**** Reading scan/item %s" % i
				self._moveto(i)
				self.scandata[i] = SpecScan(self, i, setkeys)
			rval.append(self.scandata[i]) 
		
		if len(rval) > 1:
			for i in range(len(rval)-1):
				rval[0].concatenate(rval[i+1])
			rval = [rval[0]]
		self.file.close()
		
		if len(rval) == 1:
			return rval[0]
		else:
			return rval
		
	def _getLine(self):
		"""
		Read line from datafile
		
		usage   : _getLine()
		returns : string of next line in datafile
		"""
		line = self.file.readline()
		if __verbose__ & 0x10:
				print "xxxx %s" % line.strip()
		return line
		
	def multifit(self, scans, funcs = None, mvar = None):
		
		cnt = 1
		
		results = None
		
		for i in scans:
			print i
			subplot(220 + cnt)
			
			s = self[int(i)]
				
			if funcs != None:	
				onefit = s.plot(notitles = True).fit(funcs)
			else:
				s.plot(notitles = True)
			
			if mvar != None:
				mvarval = mean(s.scandata.values[mvar])
				title("%s = %f" % (mvar, mvarval))
			if cnt == 4:
				cnt = 1
				figure()
			else:
				cnt += 1
				
			if mvar != None:
				onestack = array([mvarval])
			else:
				onestack = array([])
				
			onestack = hstack((onestack, array([i]), onefit[0,:], onefit[1,:]))
			
			if results == None:
				results = onestack
			else:
				results = vstack((results, onestack))
		
		return results
				
class SpecScan:
	def __init__(self, specfile, item, setkeys = True):	
		"""
		Read scan data from SpecData class and set variables
		to all the data contained in the scan
		
		"""
		
		# Keep track of the datafile
		
		self.datafile = specfile
		self.scandata = SpecData()
		self.scanplot = None
		self.setkeys = setkeys
		
		line = specfile._getLine()
		
		if __verbose__:
			print "---- %s" % line.strip()
		
		sline = line.strip().split()
		
		self.scan = int(sline[1])
		self.scan_type = sline[2]
		self.scan_command = ' '.join(sline[2:])
		
		self.header = line
		self.comments = ""
		
		x = 0
		self.values = {}
		self.data = array([])
		
		line = specfile._getLine()
		self.header = self.header + line
		
		#
		# Read the spec header and place the data into this class
		#
		
		while (line[0:2] != "#L") & (line != ""):
			if line[0:2] == "#P":
				# Motor positions
				pos = line.strip().split()
				for i in range(1,len(pos)):
					self.scandata.setValue(removeIllegals(specfile.motors[x]), array([float(pos[i])]))
					x += 1
					
			elif line[0:2] == "#C":
				# Comments 
				self.comments = self.comments + line
			
			elif line[0:2] == "#D":
				try:
					self.scandate = strptime(line[2:].strip())
				except:
					self.scandate = None
					
			line = specfile._getLine()
			self.header = self.header + line	
					
		if line[0:2] == "#L":
			# Comment line just before data
			self.cols = line.strip()[3:].split()
			print "---- %s" % line.strip()

		line = specfile._getLine()
		self.header = self.header + line
		
		print "---- %s" % line.strip()
		
		while (line[0:2] != "#S") & (line != "") & (line[0:4] != "# CM"):
			if line[0] != "#":
				datum = array([])
				d = line.strip().split()
				if len(d) != 0:
					for i in range(len(d)):
						v = array([float(d[i])])
						datum = concatenate((datum, v), 1)
					
					if self.data.size == 0:
						self.data = datum
					else:
						self.data = vstack((self.data, datum))
						
			elif line[0:2] == '#C':
				self.comments = self.comments + line
			else:
				self.header = self.header + line
			
			line = specfile._getLine()
			
		# Now set the motors
		
		for i in range(len(self.cols)):
			self.scandata.setValue(removeIllegals(self.cols[i]), self.data[:,i])

		# Now set the variables into the scan class from the data
		
		if self.setkeys:
			for i in self.scandata.values.keys():
				if __verbose__ & 0x02:
					print "oooo Setting variable %s" % i
				setattr(self,i , self.scandata.values[i])


		return None
	
	def concatenate(self, a):
		# Could put check in here for cols matching ?!?
		
		self.header = self.header + a.header

		self.data = vstack((self.data, a.data))
		for i in range(len(self.cols)):
			self.scandata.setValue(removeIllegals(self.cols[i]), self.data[:,i])
		if self.setkeys:
			for i in self.scandata.values.keys():
				setattr(self,i , self.scandata.values[i])
		return self

	def bin(self, a):
		return

	def plot(self, 	xcol = None, ycol = None, mcol = None, 
					norm = True, doplot = True, errors = True,
					fmt = 'ro', new = False, xint = 200, yint = 200,
					notitles = False, log = False):
		"""
		Plot the SpecScan using matplotlib
		"""
		if self.scanplot == None:
			self.scanplot = SpecPlot(self)
		
		self.scanplot.show( xcol, ycol, mcol,
							norm, doplot, errors, 
							fmt, new, xint, yint,
							notitles, log)
		return self.scanplot
	
	def fit(self, funcs, quiet = False):
		"""
		Overloaded function which calls the current scanplot
		with scanplot.fit(...)
		"""
		if self.scanplot == None:
			raise Exception("You need to plot something before trying to fit it!")
		else:
			return self.scanplot.fit(funcs, quiet)
		
		
	def show(self):
		"""
		Show statistics on scan
		"""
		print "Scan:"
		print "\t%s" % self.scan
		print "Datafile:"
		print "\t%s" % self.datafile.file.name
		self.scandata.show()
		
		
class SpecData:
	def __init__(self):
		self.values = {}
	def __call__(self, key):
		print key
		
	def setValue(self, key, data, setdata = True):
		self.values[key] = data
		if __verbose__ & 0x20:
			print "oooo Setting key %s" % key
	
	def get(self, key):
		if self.values.has_key(key):
			return self.values[key]
		else:
			return None
			
	def show(self, prefix = "", nperline = 6):
		"""
		Show statistics on data (motors, scalars)
		
		example   : SpecData.show()
		"""
		
		j = nperline
		print prefix + "Motors:"
		print prefix
		for i in self.values.keys():
			if self.values[i].shape == (1,):
				print "%10s " % i,
				j -= 1
				if j == 0:
					print prefix
					j = nperline
		
		if j != nperline:
			print prefix
		
		print prefix
		print prefix + "Scan Variables:"
		print prefix
		j = nperline
		for i in self.values.keys():
			if self.values[i].shape != (1,):
				print "%10s " % i,
				j -= 1
				if j == 0:
					print prefix
					j = nperline
		
		return None
		

class SpecPlot:
	def __init__(self, specscan):
		self.scan = specscan
		self.plt = None
		
	def show(self, 	xcol = None, ycol = None, mcol = None, 
					norm = True, doplot = True, errors = True,
					fmt = 'ro', new = True, 
					xint = 200, yint = 200,
					notitles = False, log = False):
		"""
		Plot and display the scan
		
		'xcol', 'ycol' and 'mcol' can be either numbers (negative count from end)
		or string values for the motor names
		
		for 2D plots 'xcol' can be a list or array e.g. [0, 1]
		
		'norm = False' will suppress normalization of the data
		'doplot = False' will suppress the plotting of the figure
		'errors = False' will suppress the plotting of errorbars
		'fmt' is the format string passed onto the plot command
		'new = True' will create a new figure
		'xint = 200' will set the x size of the 2D grid to 200 pts
		'yint = 200' will set the y size of the 2D grid to 200 pts
		'log = True' will set the y axis to logscale
		"""
		
		twod = False
		x2col = None
		
		if ycol == None:
			ycol = -1
		if mcol == None:
			mcol = -2
		
		if xcol is None:
			if self.scan.scan_type.strip() == 'hklmesh':
				# Look at H K and L and see which varies the most
				# This is a cludge, surely there is a better way of 
				# doing this!
				
				v = self.scan.data[:,0:3].var(axis=0)
				v = vstack((arange(0,3), v))
				v = v.transpose()
				v = v[argsort(v[:,1],0)]
				xcol = int(v[2,0])
				x2col = int(v[1,0])
				twod = True
			elif self.scan.scan_type.strip() == 'mesh':
				x2col = 1
				xcol = 0
				twod = True
			else:
				x2col = None
				xcol = 0
		
		if (type(xcol) == list) | (type(xcol) == ndarray):
			x2col = xcol[1]
			xcol = xcol[0]
		else:
			if (type(xcol) != int) & (type(xcol) != str):
				raise Exception("Illegal '%s' : xcol can only be 'int', 'list' or 'ndarray'." % type(xcol))
		
				
		if type(xcol) == str:
			xcol = self.scan.cols.index(xcol)
		if type(x2col) == str:
			x2col = self.scan.cols.index(x2col)
		if type(ycol) == str:
			# Not a number so it is probably a name
			ycol = self.scan.cols.index(ycol)
		if type(mcol) == str:
			mcol = self.scan.cols.index(mcol)
		
		
		if x2col != None:
			twod = True
		
		self.plotx = self.scan.data[:,xcol]
		if x2col != None:
			self.plotx = vstack((self.plotx, self.scan.data[:,x2col]))
			self.plotx = self.plotx.transpose()
			
		if __verbose__:
			print "**** Plotting scan %s (%s)" % (self.scan.scan, self.scan.scan_command)
			print "---- x  = %s" % self.scan.cols[xcol]
			if x2col != None:
				print "---- x2 = %s" % self.scan.cols[x2col]
			
		if norm == True:
			self.ploty = self.scan.data[:,ycol] / self.scan.data[:,mcol]
			self.plote = self.ploty * (sqrt(self.scan.data[:,ycol]) / self.scan.data[:,ycol])
			yl = "%s / %s" % (self.scan.cols[ycol], self.scan.cols[mcol])
			if __verbose__:
				print "---- y  = %s / %s" % (self.scan.cols[ycol], self.scan.cols[mcol])
		else:
			self.ploty = self.scan.data[:,ycol]
			self.plote = sqrt(self.ploty)
			yl = self.scan.cols[ycol]
			if __verbose__:
				print "---- y  = %s" % self.scan.cols[ycol]
		if log == True:
			plotit = semilogy
		else:
			plotit = plot

		if twod:
			try:
				from scipy.delaunay import Triangulation
			except ImportError:
				try:
					from delaunay import Triangulation
					#from scipy.sandbox.delaunay import Triangulation
				except ImportError:
					print """
********************                   ********************
******************** ERROR ERROR ERROR ********************
********************                   ********************

This module uses the delaunay routines to do interpolation
you will have to enable this module and re-install scipy
to be able to use it.
3D plotting will be disabled."""
					return
	
	
		# Plot the data
		if doplot:
			if new:
				pylab.figure()
				
			hold(False)
			
			if twod:
				
				xi, yi = mgrid[min(self.plotx[:,0]):max(self.plotx[:,0]):complex(0, xint)
							  ,min(self.plotx[:,1]):max(self.plotx[:,1]):complex(0, yint)]
				print "---- xint = %d" % xint
				print "---- yint = %d" % yint
				
				tri = Triangulation(self.plotx[:,0], self.plotx[:,1])
				interp = tri.nn_interpolator(self.ploty)
				zi = interp(xi,yi)
				
				pcolor(xi,yi,zi, shading = 'interp')
				
				xlim([min(self.plotx[:,0]), max(self.plotx[:,0])])
				ylim([min(self.plotx[:,1]), max(self.plotx[:,1])])
			
				if not notitles:
					xlabel(self.scan.cols[xcol])
					ylabel(self.scan.cols[x2col])
			else:
				
				if errors:
					hold(True) # Needed for errorbar plots
					self.plt = errorbar(	self.plotx, self.ploty, self.plote, 
										ecolor = 'red', fmt = fmt)
				else:
					self.plt = plot(self.plotx, self.ploty, fmt = fmt)
				
				grid(True)
				
				if not notitles:
					xlabel(self.scan.cols[xcol])
					ylabel(yl)
				
				xlim((min(self.plotx), max(self.plotx)))
			
			if not notitles:
				title("%s [%d]\n%s" % (self.scan.datafile.file.name, self.scan.scan, self.scan.scan_command))
			self.plotted = True
		else:
			self.plotted = False

		
	def prints(self, printer = None):
		"""
		Print the current figure.
		
		lpcommand should be either lp or lpr
		printer will overide the system default printer
		"""
		from tempfile import mkdtemp
		
		if self.plt == None:
			self.show()
			
		tmpdir = mkdtemp()
		
		fig = tmpdir + os.sep + "figure.eps"
		savefig(fig)
		
		texfile = open(tmpdir + os.sep + "figure.tex", "w") 
		
		texfile.write("\documentclass[12pt]{article}\n")
		texfile.write("\usepackage{graphicx}\n")
		texfile.write("\\begin{document}\n")
		texfile.write("\\begin{center}\n")
		texfile.write("\includegraphics[width=0.7\columnwidth]{%s}\n" % fig)
		texfile.write("\end{center}\n")
		texfile.write("\end{document}\n")

		texfile.close()
		
		dvips_opts = ""
		
		if printer is not None:
			dvips_opts = "-P%s" % printer
		
		os.system("cd %s ; latex figure.tex ; dvips %s figure; lp figure.ps" % (tmpdir, dvips_opts))
		os.system("rm -rfv %s" % tmpdir)

	def fit(self, funcs, quiet = False):
		
		f = fit(x = self.plotx, y = self.ploty, funcs = funcs)
		plsq = f.go()
		
		# Generate a new values of x (500) to make plotted
		# functions look good!
		
		step = ( self.plotx.max() - self.plotx.min() ) / 500
		x = arange(self.plotx.min(), self.plotx.max(), step)
			
		self.fitx = x
		self.fity = f.evalfunc(x = x)
		
		if self.plotted:	
			hold(True)
			plot(self.fitx,self.fity, 'b-')
			
		return plsq

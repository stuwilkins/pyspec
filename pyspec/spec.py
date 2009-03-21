# Written by S. Wilkins
# (c) Stuart Wilkins 2007, 2008, 2009
#
# SPEC is a Diffractometer control program from
# Certified Scientific Software.
# http://www.certif.com/
#
"""Python module for handling of spec data files (pyspec).

This module defines classes, objects and routines for reading SPEC datafiles
into python classes. There are three main classes for the data, SpecDataFile, 
SpecScan and SpecData and are arranged as follows:

SpecDataFile, contains a single data file. The datafile is indexed upon
initializing the class to allow for fast searching of large data files.

To retrieve a scan from the datafile, an item is requested from the object. For
example;

    sd = SpecDataFile('data.01')
    scan = sd[100]

The construct sd[100] will return scan 100 as a SpecScan object. By passing a range
of integers, for example;

    scans = sd[[100, 101]]

these data will be concatenated into a single object. This is useful, for example,
when two scans are used to take the data (one wide, one close) and you need to fit
the data.

The SpecScan object contains members for all motors and counters for the scan in
question. These can be accessed in two main ways. For example, to obtain the 
detector counts for the scan as an array;

   Det = sd[1000].Detector

or

   scan = sd[1000]
   Det = scan.Detector

can be used. For situations where you need machine adjustable code, the construct

   Det = scan.values['Detector']

can be used, where the member 'values' is a python dictionary of all scan variables.
In addition to the counters and motors, values defined in the header are also included.
See the SpecScan class documentation for further information. 
							
Finally, the SpecScan class defines some helper routines for plotting, fitting etc. 
which can be used to help in data analysis. See the documentation on the SpecScan class
for more information.			
										
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

def splitSpecString(ips):
	"""Split a spec string which is formated with two spaces"""

	ops = []
	for o in ips.split('  '):
		if o != '':
			ops.append(o.strip())

	return ops
	
class SpecDataFile:
	""" DataFile class for handling spec data files
	
	"""
	def __init__(self,fn):
		"""Initialize SpecDataFile

		Params
		------
		fn : string
		     Filename of spec file to open

		Returns
		-------
		SpecDataFile object

		"""
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
		Read the spec header from the datafile.

		Currently supported header items:
			'#O'	(Motor positions)
		"""
		
		self.file = open(self.filename, 'rb')
		
		if __verbose__:
			print "---- Reading Header."
		
		self.motors = []
		self.file.seek(0,0)
		line = self.file.readline()
		while line[0:2] != "#S":
			if line[0:2] == "#O":
				#self.motors = self.motors + line.strip()[4:].split()
				self.motors = self.motors + splitSpecString(line[4:])
			line = self.file.readline()
			
		self.file.close()
		return

	def index(self):
		"""Index the datafile 

		This routine indexes and sorts the byte-offests for
		all the scans (Lines beginning with '#S')
		
		"""
		self.file = open(self.filename, 'rb')
		
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
		""" Returns string with statistics on specfile."""
		
		string = ""
		string = string + head + "Specfile contains %d scans\n" % len(self.findex)
		string = string + head + "Start scan = %d\n" % min(self.findex.keys())
		string = string + head + "End   scan = %d\n" % max(self.findex.keys())
		
		return string
		
	def _moveto(self, item):
		"""Move to a location in the datafile for scan"""
		
		if self.findex.has_key(item):
			self.file.seek(self.findex[item])
		else:
			raise Exception("Scan %s is not in datafile ....." % item)
			
	def __getitem__( self, item):
		"""Convinience routine to use [] to get scan"""
		return self.getScan(item, setkeys = True)
		
	def getScan(self, item, setkeys = True):
		"""Get a scan from the data file

		This routine gets a scan from the data file and loads it into the
		list of SpecData instances. 
		
		setting 'setkeys = True' will set attributes to the
		specdata item of all the motors and counters
		
		Returns the ScanData object corresponding to the scan requested.

		"""
		if type(item) == int:
			items = [item]
		elif type(item) == list:
			items = item
		else:
			raise Exception("item can only be <int> or <list>")
		
		self.file = open(self.filename, 'rb')
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
		"""Read line from datafile"""
	
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
	"""Class defining a SPEC scan

	This class defines a single spec scan, or a collection of scans either binned or concatenated.
	The class has members for the data and variables contained in the scan. If the optional 'setkeys'
	is defined (see __init__) then the motor positions and variables will be set as class members. For
	example, in a typical scan the motor position TTH can be accessed as a class member by SpecScan.TTH

	This object has some standard members: 
	
	header     : [string]  Header of spec scan
	values     : [dict]    All scan data variables (Motors and Counters)
	data       : [array]   All scan data (cols = variables, rows = datum)
	comments   : [string]  Comments inserted in scan
	scandate   : [time]    Date of start of scan
	
	There are a number of special members. Wherever possible the following members are added to 
	the class:

	Qvec       : [array]   Q Vector
	alphabeta  : [array]   array([alpha, beta])
	wavelength : [float]   incident wavelength
	omega      : [float]   omega ((TTH / 2) - TH)
	azimuth    : [float]   azimuthal angle

	"""
	
	def __init__(self, specfile, item, setkeys = True):	
		"""Read scan data from SpecFile

		Initialize the SpecScan class from a SpecData instance.

		Parameters
		----------
		specfile       : [SpecFile] Instance of a spec data file.
		item           : [item]     Scan item.
		setkeys        : [bool]     If true set items of the class
		                            from all variables.


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
			elif line[0:2] == "#G4":
				try:
					# Geometry data
					pos = line.strip().split()
					self.Qvec = array([float(line[1]), float(line[2]), float(line[3])])
					self.alphabeta = array([float(line[5]), float(line[6])])
					self.wavelength = float(line[4])
					self.omega = float(line[7])
					self.azimuth = float(line[8])
				except:
					print "Unable to read geometry information"

			line = specfile._getLine()
			self.header = self.header + line	
					
		if line[0:2] == "#L":
			# Comment line just before data
			#self.cols = line[3:].strip().split()
			self.cols = splitSpecString(line[3:])
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
		raise Exception("Not yet implemented")
		return

	def plot(self, 	*args, **kwargs):
		"""Plot the SpecScan using matplotlib"""

		if self.scanplot == None:
			self.scanplot = SpecPlot(self)
		
		self.scanplot.show(*args, **kwargs)
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
		
	def __str__(self):
		return self.show()

	def show(self):
		"""Return string of statistics on SpecScan"""
		p = ""
		p = p + "Scan:\n"
		p = p + "\t%s\n" % self.scan
		p = p + "Datafile:\n"
		p = p + "\t%s\n" % self.datafile.file.name
		p = p + self.scandata.show()
		return p
		
class SpecData:
	"""Class defining the data contained in a scan"""
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
		
	def __str__(self):
		return self.show()

	def show(self, prefix = "", nperline = 6):
		"""Return string of statistics on data (motors, scalars)"""
		
		j = nperline
		p = ""
		p = p + prefix + "Motors:\n"
		p = p + prefix
		for i in self.values.keys():
			if self.values[i].shape == (1,):
				p = p + "%10s " % i
				j -= 1
				if j == 0:
					p = p + "\n" + prefix
					j = nperline
		
		if j != nperline:
			p = p + "\n" + prefix 

		p = p + "\n"
		
		p = p + prefix + "\n"
		p = p + prefix + "Scan Variables:\n"
		p = p + prefix + "\n"
		j = nperline
		for i in self.values.keys():
			if self.values[i].shape != (1,):
				p = p + "%10s " % i
				j -= 1
				if j == 0:
					p = p + "\n" + prefix
					j = nperline
		
		p = p + "\n"
		return p
		

class SpecPlot:
	def __init__(self, specscan):
		self.scan = specscan
		self.plt = None
		
	def show(self, 	xcol = None, ycol = None, mcol = None, 
					norm = True, doplot = True, errors = True,
					fmt = 'ro', new = True, 
					xint = 200, yint = 200,
					notitles = False, log = False):
		"""Plot and display the scan
		
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

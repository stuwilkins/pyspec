from scipy import *
from pylab import gca, hold, plot
from scipy.optimize import leastsq
from scipy.odr import Model, Data, ODR
import mpfit

def fitdata(funcs, ax = None, *args, **kwargs):
	"""Force fitting of a function to graphed data

	Parameters
	----------
	funcs    : [list] list of the fit functions
	ax       : [axis] axis to fit (if none the current
                   axis is used.
	
	All the additional *args and **kwargs are passed onto
	the fit class (see fit.__init__ for details)

	Returns
	-------

	fit Object with result.
	"""

	if ax is None:
		l = gca().lines
	else:
		l = ax.lines

	if len(l) == 0:
		print "No data on graph"
		return None
		
	xdata = l[0].get_xdata()
	ydata = l[0].get_ydata()
	
	f = fit(x = xdata, y = ydata, funcs = funcs, *args, **kwargs)	
	f.go()
	
	# Now plot the results
	
	hold(True)
	x, y = f.evalfitfunc(nxpts = 200)
	plot(x, y, 'b-')
	
	return f

class fit:
	"""General class to perform non-linear least squares fitting of data.

	This class serves as a wrapper around various fit methods, requiring a standard 
	function to allow for easy fitting of data.

	Parameters:

	x [array] : x data
	y [array] : y data
	e [array] : y error 
	
	funcs [list] : A list of functions to be fitted. These functions should be formed
	like the standard fitting functions in the pyspec.fitfuncs.

	guess [array] : Array of the same dimensions as the number of parameters containing
	the initial guess. 

	quiet [bool] : If True then don't print to the console any results.

	ifix [array] : An array of the same dimension as the number of parameters containing
	a 1 if the parameter is to be fixed.
	
	optimiser [string] : 'mpfit', 'leastsq', 'ODR'
	
	Set the optimizer used for the fit. 'mpfit' is the MPFIT python module included
	in pyspec. 'leastsq' is the scipy leastsquares optimizer. 'ODR' is the Orthogonal 
	Regresion from ODRPACK.

	"""
	def __init__(self, x = None, y = None, e = None,
		     funcs = None, guess = None, 
		     quiet = False, optimizer = 'mpfit', 
		     ivary = None, ifix = None, 
		     ilimited = None, ilimits = None):
		
		self.optimizer = optimizer
		self.setFuncs(funcs)
		self.quiet = quiet
		self.setData(x, y, e)
		self.guess = guess
		self.fit_result = None
		self.ifix = ifix
		self.ivary = ivary
		self.ilimited = ilimited
		self.ilimits = ilimits
	
	def setData(self, x, y, e):
		self.datax = x
		self.datay = y
		self.datae = e
		
	def setFuncs(self, funcs):
		self.funcs = funcs
			
	def result(self):
		return self.fit_result

	def resultDict(self):
		return None
	
	def _residuals(self, p):
		"""Residuals function for scipy leastsq optimizer"""
		f = self.evalfunc(_toFullParams(p))
		return ravel(self.datay - f)
	
	def _residualsMPFIT(self, p, fjac = None):
		"""Residuals function for MPFIT optimizer"""
		f = self.evalfunc(_toFullParams(p))
		return 0, ravel(self.datay - f)

	def _modelODR(self):
		"""Model function for ODR"""
		return self.evalfunc

	def _toFullParams(p):
		self.guess[self.ifix == 0] = p
		return self.guess
		
	def _toVaryParams(p):
		fp = p[nonzero(p * self.ifix)]
		return fp

	def evalfitfunc(self, nxpts = None, p = None, x = None, mode = 'eval'):
		"""Evaluate the fit functions with the fesult of a fit
		
		Parameters
		----------

		"""
		if x is None:
			x = self.datax
			
			if nxpts is not None:
				step = ( x.max() - x.min() ) / nxpts
				x = arange(x.min(), x.max(), step)

		f = self.evalfunc(x = x, mode = mode, p = p)

		return x, f

	def evalfunc(self, p = None, x = None, mode = 'eval'):
		f = 0.0
		ps = 0
		
		if x == None:
			x = self.datax
		
		if p == None:
			p = self.guess
		
		for i in range(len(self.funcs)):
			pe = ps + len(self.funcs[i](0,0,'params'))
			f += self.funcs[i](x , p[ps:pe]) 
			ps = pe
		
		return f
	
	def fitguess(self, mode = 'guess'):
		out = array([])
		for i in range(len(self.funcs)):
			gp = self.funcs[i](self.datax, self.datay, mode = mode)
			out = concatenate((out, gp), 1)
		return out

	def _run_odr():
		linear = Model(self._modelODR)
		mydata = Data(self.datax, self.datay, 1)
		myodr = ODR(mydata, linear, 
			    beta0 = self.guess,  
			    maxit = 10000)
		
		myoutput = myodr.run()
		
		self._result = myoutput.beta
		self._stdev = myoutput.sd_beta
		
	def _run_mpfit():
		# Make up the PARINFO list
		
		for i in range(len(self.guess)):
			if self.ivary[i]:
				self.ilimited[i] = array([1, 1])
				v = self.ivary[i]/2
				self.ilimits[i] = array([self.guess[i] - v, 
							 self.guess[i] + v])
			else:
				self.ilimited[i] = array([0, 0])
				self.ilimits[i] = array([0., 0.])

		parinfo = []
		for i in range(len(self.guess)):
			pdict = {'value' : self.guess[i],
				 'fixed' : 0, 
				 'limited' : self.ilimited[i],
				 'limits' : self.ilimits[i]}
			parinfo.append(pdict)

		m = mpfit.mpfit(self._residualsMPFIT, self._guess, 
				parinfo = parinfo, quiet = 1)

		self._result = m.params

		if m.perror is not None:
			self._stdev = m.perror
		else:
			self._stdev = zeros(len(self._guess))
			
		if m.covar is not None:
			self.covar = m.covar
			self.pcor = m.covar * 0.
			for i in range(len(self._guess)):
				for j in range(len(self._guess)):
					self.pcor[i,j] = m.covar[i,j]/sqrt(m.covar[i,i]*m.covar[j,j])

	def _run_leastsq():

		plsq = leastsq(self._residuals, self._guess, Dfun = None,  full_output = 1, factor = 0.1)
		self._result = plsq[0]	# Make the stored guess the last value
		self._stdev = sqrt(diag(plsq[1].T))

	def run(self):
		self.go()

	def go(self):
		"""Start the fit"""

		# Get the initial guess by calling the functions if we 
		# have not supplied a guess.

		if self.guess is None:
			self.guess = self.fitguess()
		elif self.guess is 'graph':
			self.guess = self.fitguess('graphguess')
		else:
			raise Exception("Unknown guess method '%s'." % self.guess)

		# Now we have the guess remove the fixed parameters from the fit

		if self.ifix == None:
			self.ifix = zeros(len(self.guess))
		else:
			# Limit fitting parameters to only thoes which are nessesary
			self._params = self.params[nonzero(self.params * self.ifix)]

		if self.ivary == None:
			self.ivary = zeros(len(self.guess))

		if self.ilimited == None:
			self.ilimited = zeros((len(self.guess), 2), dtype = int)
			self.ilimits = zeros((len(self.guess), 2))
	
		if self.quiet == False:
			if self.optimizer == 'ODR':
				print "Fitting with ODRPACK"
				print "--------------------"
			elif self.optimizer == 'mpfit':
				print "Fitting with MPFIT"
				print "------------------"
			elif self.optimizer == 'leastsq':
				print "Fitting with 'scipy' leastsq"
				print "----------------------------"
			print 
			print "Initial Guess:"
			print "-------------------------------"
			print "Parameter       :  Value       "
			print "-------------------------------"
			
			k = 0
			for i in range(len(self.funcs)):
				pnames = self.funcs[i](0,0,'params')
				for j in range(len(pnames)):
					print "%-15s : % 6.4e" % (pnames[j], self.guess[k])
					k += 1
				print "-------------------------------"
				print

		if self.optimizer == 'ODR':
			self._run_odr()
		elif self.optimizer == 'mpfit':
			self._run_mpfit()
		elif self.optimizer == 'leastsq':
			self._run_leastsq()
		else:
			raise Exception("Unknown fitting optimizer '%s'" % self.optimizer)

		# If we have fixed parameters then return the full parameters list
		# and return the full list of errors

		self.guess[self.ifix == 0] = self._result
		self.stdev = zeros(len(self.guess))
		self.stdev[self.ifix == 0] = self._stdev

		if self.quiet == False:
			print "Fit results to function(s)"
			for i in range(len(self.funcs)):
				print "    %s" % self.funcs[i](0,0,'name')
			print "--------------------------------------------"
			print "Parameter       :  Value      :  Error"
			print "--------------------------------------------"
	
			k = 0
			for i in range(len(self.funcs)):
				pnames = self.funcs[i](0,0,'params')
				for j in range(len(pnames)):
					print "%-15s : % 6.4e : % 6.4e" % (pnames[j], self.guess[k], self.stdev[k])
					k += 1
				print "--------------------------------------------"
			if self.optimizer == 'ODR':
				for reason in myoutput.stopreason:
  					print 'ODR Stop = %s' % reason
			if self.optimizer == 'mpfit':
				print 'MPFIT Status = %s' % m.statusNiceText[m.status]
				print 'MPFIT Warning = %s' % m.errmsg
				print 'MPFIT computed in %d iterations' % m.niter
		
		self.fit_result = vstack((self.guess, self.stdev))
		
		return self.fit_result

	def chiSquared(self, norm = True, dist = 'poisson'):
		""" Return the chi^2 value for the fit

		Calculate the chi-squared value for the fit. This is defined as

		\Chi^2 = \Sum_N (x_{i,n} - x_{m,i})

		Where the normalized chi^2 is given by

		\Chi^2_{norm} = \Chi^2 / M

		where M = N - P, where P is the number of parameters.

		If dist is 'poisson' then the data is divided by the model answer.
		i.e.

		\Chi^2_{poisson} = \Sum_N ((x_{i,n} - x_{m,i}) /  x_{m,i})
		
		Parameters
		----------
		norm     : [bool] return normalized chi^2
		dist     : [string] distribution

		Returns
		-------
		chi^2 value

		"""
		
		N = len(self.datax)
		P = len(self.guess)
		self.chi2 = pow(self.evalfunc() - self.datay, 2)
		if dist == 'poisson':
			self.chi2 = self.chi2 / self.evalfunc()
		self.chi2 = self.chi2.sum()
		if norm:
			self.chi2 = self.chi2 / (N - P)
		
		return self.chi2

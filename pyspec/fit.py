#
# fit.py (c) Stuart B. Wilkins 2008
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
"""FIT - Part of the pyspec routines for fitting of data

These routines / classes provide a method for fitting of data using mostly least 
squares methods. There are two main methods here. The "fit" class provides a 
class for fitting of data. The "fitdata" subroutine serves as a wrapper around
the "fit" class. Here the data is taken from the current selected figure.
Basically this means that you can fit any set of data which is on the current figure.

Examples:

f = fit(x = [xdata], y = [ydata], funcs = [func1, func2])
f.run()

result = f.result
errors = f.stdev



"""

from scipy import *
from pylab import gca, hold, plot
from scipy.optimize import leastsq
from scipy.odr import Model, Data, ODR
import mpfit

def fitdata(funcs, ax = None, showguess = False, *args, **kwargs):
   """Force fitting of a function to graphed data
   
   Parameters
   ----------
   funcs    : [list] list of the fit functions
   ax       : [axis] axis to fit (if none the current
   axis is used.
   
   All the additional *args and **kwargs are passed onto
   the fit class (see fit.__init__ for details). 
   
   Returns
   -------

   "fit" Object with result.
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
   
   hold('on')
   if showguess:
      x, y = f.evalfitfunc(nxpts = 200, p = f.guess)
      plot(x, y, 'r-')
   
   x, y = f.evalfitfunc(nxpts = 200)
   plot(x, y, 'g-')

   return f

class fit:
   """General class to perform non-linear least squares fitting of data.
   
   This class serves as a wrapper around various fit methods, requiring a standard 
   function to allow for easy fitting of data.

   Class members
   -------------

   After a sucsesful fit, the class will contain members for the following.
  
   {fit}.result       : result of fit
   {fit}.stdev        : errors on results (sigma)
   {fit}.fit_result   : an array of both results and stdev for convinience
   {fit}.covar        : covarience matrix of result
   {fit}.corr         : correlation matrix of result
   {fit}.r2           : R^2 value for fit.  
   

   """
   def __init__(self, x = None, y = None, e = None,
		funcs = None, guess = None, 
		quiet = False, optimizer = 'mpfit', 
		ivary = None, ifix = None, 
		ilimited = None, ilimits = None,
		xlimits = None, xlimitstype = 'world'):
      """
      Parameters:
      -----------
      x [array]               : x data
      y [array]               : y data
      e [array]               : y error 
      funcs [list]            : A list of functions to be fitted. 
                                These functions should be formed like the standard 
                                fitting functions in the pyspec.fitfuncs.
      guess [array|string]    : Array of the initial guess, or string for guess type:
                                'auto'  : Auto guess of parameters
      quiet [bool]            : If True then don't print to the console any results.
      ifix [array]            : An array containing a '1' for fixed parameters
      ivary [array]           : An array containing either
                                0        = Float parameter
				non zero = Vary parameter only by a max of this value
      xlimits [array]         : An (n x 2) array of the limits in x to fit.
      xlimitstype [string]    : Either 'world' or 'index'
      optimiser [string]      : 'mpfit', 'leastsq', 'ODR'
        
      """
      self.optimizer = optimizer
      self.quiet = quiet
                
      self.guess = self._checkForArray(guess)
      self.result = None
      self.fit_result = None
      
      self.ifix = self._checkForArray(ifix)
      self.ivary = self._checkForArray(ivary)
      self.ilimited = self._checkForArray(ilimited)
      self.ilimits = self._checkForArray(ilimits)

      if xlimits != None:
	 self.xlimits = self._checkForArray(xlimits)
      else:
	 self.xlimits = None
                        
      self.xlimitstype = xlimitstype 
	 
      self.setFuncs(funcs)
      self.setData(x, y, e)

   def _checkForArray(self, a):
      if type(a) == list:
	 return array(a)
      else:
	 return a
        
   def setData(self, x, y, e):

      self.datax = x
      self.datay = y
      if e is None:
	 self.datae = zeros(len(x))
      else:
	 self.datae = e
	    
      if self.xlimits is not None:
                
	 # For fitting we use the self._xdata so apply the limits

	 if len(self.xlimits.shape) == 1:
	    self.xlimits = array([self.xlimits])
	 if len(self.xlimits.shape) != 2:
	    raise Exception("Invalid array of x limits.")

	 if self.xlimitstype == 'world':
	    self._datax = self.datax.copy()
	    self._datay = self.datay.copy()
	    self._datae = self.datae.copy()
	    for limit in self.xlimits:
	       mask = (self._datax > limit[0])
	       mask = mask & (self._datax < limit[1])
	       self._datax = self._datax[mask]
	       self._datay = self._datay[mask]
	       self._datae = self._datae[mask]
                                        
	 elif self.xlimitstype == 'index':
	    self._datax = array([])
	    self._datay = array([])
	    self._datae = array([])
	    for limit in self.xlimits:
	       self._datax = concatenate((self._datax, self.datax[limit[0]:limit[1]]))
	       self._datay = concatenate((self._datay, self.datay[limit[0]:limit[1]]))
	       self._datae = concatenate((self._datae, self.datae[limit[0]:limit[1]]))
	 else:
	    raise Exception('Invalid x limits mode %s.' % self.xlimitstype)
      else:
	 self._datax = self.datax
	 self._datay = self.datay
	 self._datae = self.datae
	 
   def setFuncs(self, funcs):
      self.funcs = funcs
                        
   def result(self):
      return self.fit_result

   def resultDict(self):
      return None

## 
## The functions called by the fitting routines
##
        
   def _residuals(self, p):
      """Residuals function for scipy leastsq optimizer"""
      f = self.evalfunc(self._toFullParams(p), x = self._datax)
      return ravel(self._datay - f)

   def _residualsMPFIT(self, p, fjac = None):
      """Residuals function for MPFIT optimizer"""
      f = self.evalfunc(self._toFullParams(p), x = self._datax)
      return 0, ravel(self._datay - f)

   def _modelODR(self, p = None, x = None):
      """Model function for ODR"""
      return self.evalfunc(self._toFullParams(p), x = x)
   
   def _toFullParams(self, p):
      """Return the full parameter list
      
      Sets the full parameter list, substituting the guess for 
      the fixed parameters
      """
      g = self.guess.copy()
      g[nonzero(self.ifix == 0)] = p
      return g
        
##
## Functions to evaluate the fitting functions
##
        
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
      if x is None:
	 x = self.datax
                
      if p is None:
	 if self.result is not None:
	    p = self.result
	 else:
	    p = self.guess
                
      for i in range(len(self.funcs)):
	 pe = ps + len(self.funcs[i](0,0,'params'))
	 f += self.funcs[i](x , p[ps:pe]) 
	 ps = pe
      return f

   def fitguess(self, mode = 'guess'):
      """Evaluate the fit functions to get the initial guess"""
      out = array([])
      for i in range(len(self.funcs)):
	 gp = self.funcs[i](self.datax, self.datay, mode = mode)
	 out = concatenate((out, gp), 1)
      return out

##
## Routines to run the different optimizers
##

   def _run_odr(self):
      """Run an ODR regression"""
      linear = Model(self._modelODR)
      mydata = Data(self._datax, self._datay, 1)
      myodr = ODR(mydata, linear, 
		  beta0 = self._guess,  
		  maxit = 10000)
      
      myoutput = myodr.run()
      
      self._result = myoutput.beta
      self._stdev = myoutput.sd_beta
      self._covar = myoutput.cov_beta
      self._odr = myoutput
                
   def _run_mpfit(self):
      """Run an MPFIT regression"""
      # Make up the PARINFO list
      
      for i in range(len(self.guess)):
	 if self.ifix[i] == 0:
	    if self.ivary[i]:
	       self.ilimited[i] = array([1, 1])
	       v = abs(self.ivary[i]/2)
	       self.ilimits[i] = array([self.guess[i] - v, 
					self.guess[i] + v])
	    #else:
	    #   self.ilimited[i] = array([0, 0])
	    #   self.ilimits[i] = array([0., 0.])
                                        
      parinfo = []
      for i in range(len(self.guess)):
	 if self.ifix[i] == 0:
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
	 self._covar = m.covar
      else:
	 self._covar = zeros((len(self._guess), len(self._guess)))

      self._mpfit = m

   def _run_leastsq(self):
      """Run a scipy leastsq regression"""
      plsq = leastsq(self._residuals, self._guess, Dfun = None,  full_output = 1, factor = 0.1)
      self._result = plsq[0]  # Make the stored guess the last value
      self._stdev = sqrt(diag(plsq[1].T))
      self._covar = plsq[1].T
      self._leastsq = plsq

##
## Run the optimization
##

   def run(self):
      self.go()

   def go(self):
      """Start the fit"""

      # Get the initial guess by calling the functions if we 
      # have not supplied a guess.

      if self.guess is None:
	 self.guess = self.fitguess()
      elif self.guess == 'graph':
	 self.guess = self.fitguess('graphguess')
      
      # Now we have the guess remove the fixed parameters from the fit
      
      if self.ifix is None:
	 self.ifix = zeros(len(self.guess))
	 self._guess = self.guess
      else:
	 # Limit fitting parameters to only thoes which are nessesary
	 self._guess = self.guess[nonzero(self.ifix == 0)]

      if self.ivary is None:
	 self.ivary = zeros(len(self.guess))
	 
      if self.ilimited is None:
	 self.ilimited = zeros((len(self.guess), 2), dtype = int)
	 self.ilimits = zeros((len(self.guess), 2))
        
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

      self.result = self.guess.copy()
      self.result[nonzero(self.ifix == 0)] = self._result
      self.stdev = zeros(len(self.guess))
      self.stdev[nonzero(self.ifix == 0)] = self._stdev

      # Now deal with the covarience matrix

      if self.ifix.sum():
	 self.covar = zeros((len(self.result), len(self.result)))
	 ii = 0
	 for i in range(len(self.stdev)):
	    jj = 0
	    for j in range(len(self.stdev)):
	       if self.ifix[i] or self.ifix[j]:
		  # We have a fixed paremeter
		  self.covar[i,j] = 0.0
	       else:
		  # Non fixed
		  self.covar[i,j] = self._covar[ii,jj]
		  jj+=1
	    if not self.ifix[i]:
	       ii+=1
      else:
	 self.covar = self._covar

      ##
      ## Calculate some useful
      ##

      # Correlation matrix 
      self.corr = self.covar * 0.0
      for i in range(len(self.corr)):
	 for j in range(len(self.corr)):
	    self.corr[i,j] = self.covar[i,j] / sqrt(self.covar[i,i] * self.covar[j,j])

      # Chi^2
      self.chiSquared()
      
      # R^2
      ssReg = pow(self.evalfunc() - self.datay, 2).sum()
      ymean = self.datay.sum() / len(self.datay)
      ssTot = pow(ymean - self.datay, 2).sum()
      self.r2 = 1.0 - (ssReg / ssTot)
      
      # Make the dictionaryies

      self.resultsDict = self.makeDict(self.result)

      # Print out the result to console

      if self.quiet == False:
	 print self.textResult()
                
      self.fit_result = vstack((self.result, self.stdev))          
      return self.fit_result

   def makeDict(self, values):
      """Make and return a dictionary of the parameters or errors"""

      alld = []
      # First make up the 
      for f in self.funcs:
	 d = {}
	 pnames = f(0,0,'params')
	 for (name, val) in zip(pnames, values):
	    d[name] = val
	 alld.append(d)

      return alld

   def textResult(self):
      """Return a string containing the text results of the fit"""
      p = ""
                
      if self.optimizer == 'ODR':
	 p += "Fitted with ODRPACK\n"
	 p += "--------------------\n"
      elif self.optimizer == 'mpfit':
	 p += "Fitted with MPFIT\n"
	 p +=  "------------------\n"
      elif self.optimizer == 'leastsq':
	 p += "Fitted with 'scipy' leastsq\n"
	 p += "----------------------------\n"
                
      ## Put in here the x limits

      p += "Fit results to function(s)\n"
      for i in range(len(self.funcs)):
	 p += "    %s\n" % self.funcs[i](0,0,'name')
      sep = "--------------------------------------------------------------------------------------------------\n"
      p  += sep
      p  += "Parameter       :  Init. Guess :  Value      :  Error     :  Fixed     :  -ve Limit :  +ve Limit :\n"
      p  += sep
        
      yesno = ['NO', 'YES']
	    
      k = 0
      for i in range(len(self.funcs)):
	 pnames = self.funcs[i](0,0,'params')
	 for j in range(len(pnames)):
	    p += "%-15s : % 6.4e : % 6.4e : % 6.4e :  %-6s    :" % (pnames[j], self.guess[k], self.result[k], self.stdev[k], yesno[int(self.ifix[k])])
	    for n in range(2):
	       if self.ilimited[k, n]:
		  p += " %6.4e :" % self.ilimits[k,n]
	       else:
		  p += "  %-6s    :" % 'None'
	    p += "\n"
	    k += 1
	 p += sep

      # Display the normalized covarience matrix

      allpnames = []
      for f in self.funcs:
	 allpnames += f(0,0,'params')

      p += "\nCorrelation :\n"
      p += sep
      rsidsdev = self.residualSDev()
      for i in range(len(self.covar)):
	 for j in range(i,len(self.covar)):
	    if i != j:
	       if self.covar[i,j] != 0:
		  p += "%-15s <=> %-15s = %6.4e\n" % (allpnames[i], allpnames[j], self.corr[i,j])

      # Display the fit parameters

      p += "\nGoodness of fit:\n"
      p += sep
      p += "%-10s = %6.4e\n" % ("Chi^2", self.chi2)
      p += "%-10s = %f\n" % ("R^2", self.r2)

      p += sep

      # Now display optimizer specific code

      if self.optimizer == 'ODR':
	 for reason in self._odr.stopreason:
	    p += 'ODR Stop = %s\n' % reason
      if self.optimizer == 'mpfit':
	 p += 'MPFIT Status = %s\n' % self._mpfit.statusNiceText[self._mpfit.status]
	 p += 'MPFIT Warning = %s\n' % self._mpfit.errmsg
	 p += 'MPFIT computed in %d iterations\n' % self._mpfit.niter
      return p

   def chiSquared(self, norm = True, dist = 'poisson'):
      """ Return the chi^2 value for the fit
   
      Calculate the chi-squared value for the fit. This is defined as
   
      \Chi^2 = \Sum_N (x_{d,n} - x_{m,n})
   
      Where d is the data and m is the model.
      The normalized chi^2 is given by
   
      \Chi^2_{norm} = \Chi^2 / M
   
      where M = N - P, where P is the number of parameters.
   
      If dist is 'poisson' then the data is divided by the model answer.
      i.e.
      
      \Chi^2_{poisson} = \Sum_N ((x_{d,n} - x_{m,i}) /  x_{m,i})
   
      Parameters
      ----------
      norm     : [bool] return normalized chi^2
      dist     : [string] distribution
      
      Returns
      -------
      chi^2 value
      
      """
      
      N = len(self.datax)
      P = len(self.result)
      self.chi2 = pow(self.evalfunc() - self.datay, 2)
      if dist == 'poisson':
	 self.chi2 = self.chi2 / self.evalfunc()
      
      self.chi2 = self.chi2.sum()
      if norm:
	 self.chi2 = self.chi2 / (N - P)
	 
      return self.chi2

   def residualSDev(self):
      """Calculate the sandard deviation of the residuals"""
      resid = self.evalfunc() - self.datay
      mean = resid.sum() / len(resid)
      stdev = sqrt(pow(resid - mean, 2).sum() / len(resid))
      self.resid_sdev = stdev 
      return stdev

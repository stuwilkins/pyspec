from scipy import *
from pylab import gca, hold, plot
from scipy.optimize import leastsq
from scipy.odr import Model, Data, ODR


def fitdata(funcs = None, guess = None, ifix=None):
	
	l = gca().lines 
	if len(l) == 0:
		print "No data on graph"
		return None
		
	xdata = l[0].get_xdata()
	ydata = l[0].get_ydata()
	
	f = fit(x = xdata, y = ydata,
			funcs = funcs, guess = guess, ifix = ifix)
	
	f.go()
	
	# Now plot the results
	
	hold(True)
	x, y = f.evalfitfunc(nxpts = 200)
	plot(x, y, 'b-')
	
	return f

class fit:

	def __init__(self, 	x = None, y = None, e = None,
						funcs = None, guess = None, 
						quiet = False, optimizer = 'leastsq', ifix = None):
		
		self.optimizer = optimizer
		self.setFuncs(funcs)
		self.quiet = quiet
		self.setData(x, y, e)
		self.guess = guess
		self.fit_result = None
		self.ifix = ifix
	
	def setData(self, x, y, e):
		self.datax = x
		self.datay = y
		self.datae = e
		
	def setFuncs(self, funcs):
		self.funcs = funcs
			
	def result(self):
		return self.fit_result
			
	def _residuals(self, p):
		f = self.evalfunc(p)
		return ravel(self.datay - f)
		
		
	def _dfdp(self, p, y, x, func):
		m=len(x)
		n=len(p)      #dimensions
		ps=p
		prt=zeros([m,n], dtype=float64)
		delta=zeros(n, dtype=float64);        # initialise Jacobian to Zero

		for j in range(n):
			delta[j]=p[j] * 0.01;    #cal delx=fract(dp)*param value(p)
			print delta[j]
			if p[j] == 0:
				delta[j] = 0.01     #if param=0 delx=fraction
		
			p[j] = ps[j] + delta[j]
		
			f1 = fitfunc(func, x, p)
		
			if delta[j] != 0:
				prt[:,j] = (f1 - y) / delta[j]
			p[j]=ps[j];
		
		return prt

	def evalfitfunc(self, nxpts = None, p = None, x = None, mode = 'eval'):
		
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
		print self.funcs
	
		out = array([])
		for i in range(len(self.funcs)):
			gp = self.funcs[i](self.datax, self.datay, mode = mode)
			out = concatenate((out, gp), 1)
		return out

	def go(self):
	
	
		if self.guess is None:
			self.guess = self.fitguess()
		if self.guess is 'graph':
			self.guess = self.fitguess('graphguess')
	
		if self.quiet == False:
			if self.optimizer == 'ODR':
				print "Fitting with ODRPACK"
				print "--------------------"
			else:
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

			#hold(True)
			#x, y = self.evalfitfunc(nxpts = 200)
			#print plot(x, y, 'g-')

		if self.optimizer == 'ODR':
	
			linear = Model(self.evalfunc)
			mydata = Data(self.datax, self.datay, 1)
			myodr = ODR(mydata, linear, 
						beta0=self.guess,  
						maxit = 10000, ifixb=self.ifix)
			myoutput = myodr.run()
			
			self.guess = myoutput.beta
			self.stdev = myoutput.sd_beta
		else:
		
			plsq = leastsq(self._residuals, self.guess, Dfun = None,  full_output = 1, factor = 0.1)

			self.guess = plsq[0]	# Make the stored guess the last value
			self.stdev = sqrt(diag(plsq[1].T))
			#stdev = ones(len(guess))
	
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
			#print plsq[3]
		
		self.fit_result = vstack((self.guess, self.stdev))
		
		return self.fit_result
		
		

.. highlight:: python

============================
Fitting of Experimental Data
============================

************
Introduction
************

These routines / classes provide a method for fitting of data using mostly least
squares methods. There are two main methods here. The ``pyspec.fit.fit`` class provides a
class for fitting of data. The ``pyspec.firt.fitdata`` subroutine serves as a wrapper around the ``pyspec.fit.fit`` class. Here the data is taken from the current selected figure. Basically this means that you can fit any set of data which is on the current figure.

*****************
The ``fit`` Class
*****************

.. autoclass:: pyspec.fit.fit 
   :members:

**********************
Standard Fit Functions
**********************

.. automodule:: pyspec.fitfuncs
   :members:

**********************
User Defined Functions
**********************

User defined functions can easily be written. These should follow the example::

    def linear(x, p, mode='eval'): 
       if mode == 'eval':
          out = (p[0] * x) + p[1]
       elif mode == 'params':
          out = ['grad','offset']
       elif mode == 'name':
          out = "Linear"
       elif mode == 'guess':
          g = peakguess(x, p)
          out = g[4:6]
       else:
          out = []
       return out

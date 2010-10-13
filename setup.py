#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

#LEVMAR_LIBDIR = '/Users/swilkins/Dropbox/Programming/Python/levmar-2.5'
#LEVMAR_INCDIR = '/Users/swilkins/Dropbox/Programming/Python/levmar-2.5'
LEVMAR_LIBDIR = '/home/tardis/swilkins/C-C++/levmar-2.5'
LEVMAR_INCDIR = '/home/tardis/swilkins/C-C++/levmar-2.5'

levmar = Extension('pylevmar', ['src/pylevmar.c'],
                   libraries = ['levmar', 'm', 'blas', 'lapack'],
                   extra_compile_args = ['-g'],
                   library_dirs = ['/usr/lib', LEVMAR_LIBDIR],
                   include_dirs = [LEVMAR_INCDIR, np.get_include()],
                   depends = ['src/pylevmar.h'])

setup(name='pyspec',
      version='0.1',
      description='Python utilities for spec data analysis',
	  author='Stuart Wilkins',
	  author_email='stuwilkins@mac.com',
	  url='http://www.stuwilkins.com',
	  packages=['pyspec', 'pyspec.ccd'],
          ext_modules = [levmar]
      )

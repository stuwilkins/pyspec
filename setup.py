#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

LEVMAR_LIBDIR = '/Users/swilkins/Dropbox/Programming/Python/levmar-2.5'
LEVMAR_INCDIR = '/Users/swilkins/Dropbox/Programming/Python/levmar-2.5'

levmar = Extension('pylevmar', ['pyspec/pylevmar.c'],
                   libraries = ['levmar', 'm'],
		   extra_link_args = ['-framework vecLIB'],
                   extra_compile_args = ['-g'],
                   library_dirs = ['/usr/lib', LEVMAR_LIBDIR],
                   include_dirs = [LEVMAR_INCDIR, np.get_include()],
                   depends = ['pyspec/pylevmar.h'])

setup(name='pyspec',
      version='0.1',
      description='Python utilities for spec data analysis',
	  author='Stuart Wilkins',
	  author_email='stuwilkins@mac.com',
	  url='http://www.stuwilkins.com',
	  packages=['pyspec'],
          ext_modules = [levmar]
      )

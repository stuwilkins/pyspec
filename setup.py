#!/usr/bin/env python

from distutils.core import setup, Extension
from setupext import ext_modules

setup(name='pyspec',
      version='0.2',
      description='Python utilities for spec data analysis',
      author='Stuart Wilkins',
      author_email='stuwilkins@mac.com',
      url='http://www.stuwilkins.com',
      packages=['pyspec', 'pyspec.ccd'],
      ext_modules = ext_modules
)

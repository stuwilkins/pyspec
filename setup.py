#!/usr/bin/env python
""" pyspec : Python routuines and packages for x-ray data analysis

The pyspec package contains x-ray data analysis utilities useful for 
analyzing data from SPEC and CCD cameras installed at synchrotrons.  

Non-linear least squares regression is provided for peak fitting.

CCD transformations are provided to enable the use of a CCD detector
and enable transformations into reciprocal space coordinates.

SPEC is from Certified Scientific Software (http://www.certif.com/)
"""

classifiers = """\
Development Status :: 4 - Beta
Environment :: Console
Intended Audience :: Science/Research
License :: OSI Approved :: GNU General Public License (GPL)
Programming Language :: Python
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Chemistry
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: Microsoft :: Windows
Operating System :: MacOS :: MacOS X
Operating System :: Unix
"""

from distutils.core import setup, Extension
from setupext import ext_modules
import sys

doclines = __doc__.split("\n")

if sys.version_info < (2, 3):
    _setup = setup
    def setup(**kwargs):
        if kwargs.has_key("classifiers"):
            del kwargs["classifiers"]
        _setup(**kwargs)

setup(name='pyspec',
      version='0.2',
      author='Stuart Wilkins',
      author_email='stuwilkins@mac.com',
      maintainer="Stuart B. Wilkins",
      maintainer_email="stuwilkins@mac.com",
      url = 'http://sourceforge.net/projects/pyspec',
      platforms = ["any"],
      description = doclines[0],
      long_description = "\n".join(doclines[2:]),
      packages=['pyspec', 'pyspec.ccd'],
      ext_modules = ext_modules
)

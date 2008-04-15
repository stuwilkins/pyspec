#!/usr/bin/env python

"""Setup script for the SPS module distribution."""

import os, sys, glob

try:
    import numpy
except ImportError:
    text  = "You must have numpy installed.\n"
    text += "See http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103\n"
    raise ImportError, text

import platform
from distutils.core import setup
from distutils.extension import Extension

if platform.system() == 'Linux' :
    extra_compile_args = ['-pthread']
    #extra_compile_args = []
elif platform.system() == 'SunOS' :
    #extra_compile_args = ['-pthreads']
    extra_compile_args = []
elif platform.system() == 'Darwin' :
    extra_compile_args = ['-pthread']
    #extra_compile_args = []
else:
    extra_compile_args = []
    
ext_modules = [Extension(
                        name = 'spslut',
                        sources=['Src/sps_lut.c',
                                 'Src/spslut_py.c'],
                        extra_compile_args = extra_compile_args,
                        include_dirs  = ['Include', numpy.get_include()],
                   )]
if sys.platform == "win32":
    define_macros = [('WIN32',None)]
else:
    define_macros = []
    ext_modules.append( Extension(
                            name = 'sps',
                            sources=['Src/sps.c',
                                     'Src/sps_py.c'],
                            extra_compile_args = extra_compile_args,
                            include_dirs  = ['Include', numpy.get_include()]))
setup (
        name         = "sps",
        version      = "1.0",
        description  = "shared memory and spec",
        author       = "BLISS Group",
        author_email = "rey@esrf.fr",
        url          = "http://www.esrf.fr/computing/bliss/",

        # Description of the modules and packages in the distribution

        extra_path   = 'Pybliss',
        ext_modules  = ext_modules)

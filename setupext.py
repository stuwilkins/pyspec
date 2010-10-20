import os
import ConfigParser
import numpy as np
from distutils.core import Extension

options = {'build_levmar' : False }
levmar  = {'include_dirs' : [],
           'library_dirs' : [],
           'libraries'    : [] }

def parseExtensionSetup(name, config, default):
    
    try: default['include_dirs'] = config.get(name, "include_dirs").split(os.pathsep)
    except: pass
    try: default['library_dirs'] = config.get(name, "library_dirs").split(os.pathsep)
    except: pass
    try: default['libraries'] = config.get(name, "libraries").split(",")
    except: pass

    return default
    
if os.path.exists("setup.cfg"):
    config = ConfigParser.SafeConfigParser()
    config.read("setup.cfg")

    try: options['build_levmar'] = config.getboolean("levmar","build")
    except: pass
    try: options['build_gridder'] = config.getboolean("gridder","build")
    except: pass

    levmar = parseExtensionSetup('levmar', config, levmar)
    levmar['include_dirs'].append(np.get_include())
    levmar['libraries'].append('levmar')

ext_modules = []
if options['build_levmar']:
    ext_modules.append(Extension('pyspec.pylevmar', ['src/pylevmar.c'],
                                 libraries = levmar['libraries'],
                                 extra_compile_args = ['-g'],
                                 library_dirs = levmar['library_dirs'],
                                 include_dirs = levmar['include_dirs'],
                                 depends = ['src/pylevmar.h']))
if options['build_gridder']:
    ext_modules.append(Extension('pyspec.gridder', ['src/gridder.c'],
                                 depends = ['src/gridder.h'],
				 include_dirs = [np.get_include()]))



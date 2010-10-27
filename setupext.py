import os
import ConfigParser
import numpy as np
import copy
from distutils.core import Extension

options = {'build_levmar'    : False ,
           'build_gridder'   : False ,
           'build_princeton' : False }

ext_default  = {'include_dirs' : [np.get_include()],
                'library_dirs' : [],
                'libraries'    : [] }

def parseExtensionSetup(name, config, default):
    default = copy.deepcopy(default)
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
    try: options['build_princeton'] = config.getboolean("princeton","build")
    except: pass

    levmar = parseExtensionSetup('levmar', config, ext_default)
    levmar['libraries'].append('levmar')

    gridder = parseExtensionSetup('gridder', config, ext_default)
    princeton = parseExtensionSetup('princeton', config, ext_default)


ext_modules = []
if options['build_levmar']:
    ext_modules.append(Extension('pylevmar', ['src/pylevmar.c'],
                                 extra_compile_args = ['-g'],
                                 depends = ['src/pylevmar.h'],
                                 **levmar))
if options['build_gridder']:
    ext_modules.append(Extension('gridder', ['src/gridder.c'],
                                 depends = ['src/gridder.h'],
                                 **gridder))

if options['build_princeton']:
    ext_modules.append(Extension('princeton', ['src/princeton.c'],
                                 depends = ['src/princeton.h'],
				 **princeton))



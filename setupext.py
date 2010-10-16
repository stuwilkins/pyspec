
import os
import ConfigParser
import numpy as np

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

    levmar = parseExtensionSetup('levmar', config, levmar)
    levmar['include_dirs'].append(np.get_include())


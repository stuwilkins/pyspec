
import os
import ConfigParser

options = {'build_levmar' : False }

if os.path.exists("setup.cfg"):
    config = ConfigParser.SafeConfigParser()
    config.read("setup.cfg")

    try: options['build_levmar'] = config.getboolean("levmar","build_levmar")
    except: pass

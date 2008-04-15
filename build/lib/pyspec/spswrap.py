#/*##########################################################################
# Copyright (C) 2004-2007 European Synchrotron Radiation Facility
#
# This file is part of the PyMCA X-ray Fluorescence Toolkit developed at
# the ESRF by the Beamline Instrumentation Software Support (BLISS) group.
#
# This toolkit is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option) 
# any later version.
#
# PyMCA is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# PyMCA; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA 02111-1307, USA.
#
# PyMCA follows the dual licensing model of Trolltech's Qt and Riverbank's PyQt
# and cannot be used as a free plugin for a non-free program. 
#
# Please contact the ESRF industrial unit (industry@esrf.fr) if this license 
# is a problem to you.
#############################################################################*/
try:
    import sps
    from sps import CHAR, DOUBLE, FLOAT, LONG, SHORT, UCHAR, ULONG, USHORT
    from sps import TAG_ARRAY, TAG_MCA, TAG_IMAGE, TAG_SCAN, \
                    TAG_INFO, TAG_MASK, TAG_STATUS
    from sps import IS_ARRAY, IS_MCA, IS_IMAGE, STRING
    from sps import error
    from sps import updatecounter
except:
    #windows does not use it
    pass
import threading

spsdefaultoutput ={"axistitles":   '',
                        "xlabel":       '',
                        "ylabel":       '',
                        "title":        '',
                        "nopts":        0,
                        "xbeg":         0,
                        "xend":         0,
                        "plotlist":     [],
                        "datafile":     '/dev/null',
                        "scantype":     16,
                        "aborted":      0}

spsdefaultarraylist={}

spslock = threading.Lock()

def getarrayinfo(spec,shm):
    result = [0,] * 4

    spslock.acquire()
    try:
       result = sps.getarrayinfo (spec,shm)
    except: 
       pass
    spslock.release()
    return result

def getarraylist( spec ):
 
    result = []
    if specrunning(spec):
        spslock.acquire()
        try:
           result = sps.getarraylist( spec )
           spsdefaultarraylist[spec]=result
        except:
           pass 
        spslock.release()
    else:
        if spsdefaultarraylist.has_key(spec):
            return spsdefaultarraylist[spec]
        else:
            spsdefaultarraylist[spec]=[]
            return result
    return result

def isupdated(spec, shmenv):

    result = 0

    spslock.acquire()

    try:
        result = sps.isupdated( spec, shmenv  )
    except:
        pass

    spslock.release()
    return result

def putenv(spec,shmenv,cmd,outp):

    result = None

    spslock.acquire()

    try:
       result = sps.putenv(spec,shmenv,cmd,outp)
    except:
       pass  
   
    spslock.release()

    return result

def getenv(spec,shmenv,key):
    result = ''

    spslock.acquire()

    try:
        #  if key != 'plotlist':
        result = sps.getenv(spec,shmenv,key)
    except sps.error:
        if key in spsdefaultoutput.keys():
            result = spsdefaultoutput[key]
        pass

    spslock.release()

    return result

def updatedone(spec,shmenv):
    result = 0
    
    spslock.acquire()

    try:
       result = sps.updatedone(spec,shmenv)
    except: 
       pass

    spslock.release()
    return result

def getdata(spec,shm):

    result = []

    spslock.acquire()
    try:
        result = sps.getdata(spec,shm)
    except:
        pass
    spslock.release()
    return result

def getdatacol(spec,shm,idx):

    result = []

    spslock.acquire()
    try:
       result = sps.getdatacol(spec,shm,idx)
    except:
       pass
    spslock.release()
    return result

def getdatarow(spec,shm,idx):
    result = []

    spslock.acquire()
    try:
       result = sps.getdatarow(spec,shm,idx)
    except:
       pass
    spslock.release()
    return result

def getspeclist():

    result = []

    spslock.acquire()
    try:
       result = sps.getspeclist()
    except:
       pass
    spslock.release()
    return result

def getkeylist(spec,shmenv):

    result = []

    spslock.acquire()

    try:
        result = sps.getkeylist(spec,shmenv)
    except:
        pass

    spslock.release()

    return result

def specrunning(spec):
    spec_list = getspeclist()
    if spec not in spec_list:
        return 0
    else:
        return 1


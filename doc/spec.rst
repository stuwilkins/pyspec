Module to process SPEC data files
*********************************

Introduction
============

This module defines classes, objects and routines for reading SPEC datafiles
into python classes. There are three main classes for the data, SpecDataFile,
SpecScan and SpecData and are arranged as follows:

SpecDataFile, contains a single data file. The datafile is indexed upon
initializing the class to allow for fast searching of large data files.

To retrieve a scan from the datafile, an item is requested from the object. For
example::

    >>sd = SpecDataFile('data.01')
    >>scan = sd[100]

The construct sd[100] will return scan 100 as a SpecScan object. By passing a range
of integers, for example::

    scans = sd[[100, 101]]

these data will be concatenated into a single object. This is useful, for example,
when two scans are used to take the data (one wide, one close) and you need to fit
the data.

The SpecScan object contains members for all motors and counters for the scan in
question. These can be accessed in two main ways. For example, to obtain the
detector counts for the scan as an array::

   >>>Det = sd[1000].Detector

or::

   >>>scan = sd[1000]
   >>>Det = scan.Detector

can be used. For situations where you need machine adjustable code, the construct::

   >>>Det = scan.values['Detector']

can be used, where the member 'values' is a python dictionary of all scan variables.
In addition to the counters and motors, values defined in the header are also included.
See the SpecScan class documentation for further information.

Finally, the SpecScan class defines some helper routines for plotting, fitting etc.
which can be used to help in data analysis. See the documentation on the SpecScan class
for more information.
These classes form the basis for the SPEC datafile handling. There are three objects, 
the SpecDataFile the SpecScan and the SpecData.

Spec Data File Class
====================

.. autoclass:: pyspec.spec.SpecDataFile
   :members:

Spec Scan Class
===============

.. autoclass:: pyspec.spec.SpecScan
   :members:

Spec Data Class
===============

.. autoclass:: pyspec.spec.SpecData
   :members:




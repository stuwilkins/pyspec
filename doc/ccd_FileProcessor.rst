************************
File Processing Routines
************************

Introduction
------------

The ``FileProcessor`` class converts a list of filenames (or a scan
object containing from which a list of filenames can be found) into a
3D array, or stack of images. These images are corrected by dark
images and are normalized by the beam monitor.

Using the ``FileProcessor``
---------------------------

The ``FileProcessor`` is a class and can be called as::

   >>>fp = FileProcessor()

which creates a blank image processor object. In order to process
images the user must provide a list of filenames for both the images
and dark images. If this list is 1D then each image corresponds to one
image in the stack. If the list is 2D (or in python-speak a list of
lists) then the inner list is taken as a number of sub images which
should be binned together. This is true for dark images also. The
images can be set by::

   >>>

``FileProcessor`` Documentation
------------------------------- 

.. autoclass:: pyspec.ccd.transformations.FileProcessor
   :members:


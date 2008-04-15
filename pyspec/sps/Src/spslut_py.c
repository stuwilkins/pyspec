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
/* spslut_py.c VERSION 4.0 */
/* $Revision: 1.7 $
/* $Log: spslut_py.c,v $
 * Revision 1.7  2005/02/10 23:37:48  sole
 * minor changes
 *
 * Revision 1.6  2005/02/10 17:44:58  sole
 * *** empty log message ***
 *
 * Revision 1.5  2005/02/10 16:17:15  sole
 * Removed some unused variables
 **/
/* CHANGES:

	[05-09-2002] A. Gobbo
	- Included min and max values to 8 bit colormaps 
  
	[11-03-2002] A. Gobbo
	- Included modes BGR and BGRX 

	[12-12-2001] A. Gobbo
	- Dimentions inverted in the returned array
	- Corrected memory leak bug     
*/
#include <stdio.h>      

#include <Python.h> 
#include <numpy/oldnumeric.h>
#include <sps_lut.h>

static PyObject *SPSLUTError;

PyObject *new_pyimage(const char *mode, unsigned xsize, unsigned ysize, 
		      void *data)
{
  return PyString_FromStringAndSize ((const char *)data, 
				     strlen(mode) * xsize * ysize);
}
   
static int natbyteorder()
{
  union {
    struct {
      unsigned char b1;
      unsigned char b2;
      unsigned char b3;
      unsigned char b4;
    } c;
    unsigned long p;
  } val;

  val.p = 1;
  if (val.c.b4 == 1) {
    return SPS_MSB;
  } else {
    return SPS_LSB;
  }
}

static PyObject *spslut_transform(self, args)
     PyObject *self, *args;
     
{
  char *key, *value;
  void *data;
  int type, cols, rows, reduc, fastreduc, meth, autoscale, mapmin=0, mapmax=255;
  int palette_code;
  double gamma, min, max;
  XServer_Info Xservinfo;
  void *palette;
  int prows, pcols, pal_entries;
  int byteorder, pixelsize;
  unsigned long redmask, greenmask, bluemask;
  void *r/*, *res*/;
  char *mode;
  PyArrayObject *src;
  PyObject *in_src;
  PyObject *res,*aux;
  int array_output=0;
  unsigned char *as_pointer, *as_r;
  int as_dim[3];
  PyArrayObject *as_aux;
  
  if (!PyArg_ParseTuple(args, "O(ii)(id)sii(dd)|(ii)i", &in_src, &reduc, 
			&fastreduc, &meth, &gamma, &mode, &palette_code,
			&autoscale, &min, &max,&mapmin, &mapmax, &array_output)) 
		return NULL;
  
  if (strcmp(mode, "RGB") == 0) {
    Xservinfo.red_mask = 0x0000ff;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0xff0000;
    Xservinfo.pixel_size = 3;
    Xservinfo.byte_order = natbyteorder();
  } else if (strcmp(mode, "RGBX") == 0) {
    Xservinfo.red_mask = 0x0000ff;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0xff0000;
    Xservinfo.pixel_size = 4;
    Xservinfo.byte_order = natbyteorder();
  } 
  
  //###CHANGED - ALEXANDRE 11/03/2002 - Qt uses different order than Tkinter
  else if (strcmp(mode, "BGR") == 0) {
    Xservinfo.red_mask = 0xff0000;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0x0000ff;
    Xservinfo.pixel_size = 3;
    Xservinfo.byte_order = natbyteorder();
  } else if (strcmp(mode, "BGRX") == 0) {
    Xservinfo.red_mask = 0xff0000;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0x0000ff;
    Xservinfo.pixel_size = 4;
    Xservinfo.byte_order = natbyteorder();
  } 

    else if (strcmp(mode, "L") == 0 || strcmp(mode, "P") == 0  ) {
    Xservinfo.pixel_size = 1;
    Xservinfo.byte_order = natbyteorder();
    //mapmin = 0;
    //mapmax = 255;
  } else {
    PyErr_SetString(SPSLUTError, "Mode must be RGB, RGBX, BGR, BGRX, L or P");
    return NULL;
  }

  if (!(src = (PyArrayObject*) PyArray_ContiguousFromObject(in_src, 
   			    PyArray_NOTYPE, 2, 2))) {
    PyErr_SetString(SPSLUTError, "Input Array could is not a 2x2 array");
    return NULL;
  }

  switch (src->descr->type_num) {
  case PyArray_UINT:
    type = SPS_ULONG; break;
  case PyArray_USHORT:
    type = SPS_USHORT; break;
  case PyArray_LONG: 
  case PyArray_INT: 
    type = SPS_LONG; break;
  case PyArray_SHORT: 
    type = SPS_SHORT; break;
  case PyArray_UBYTE: 
    type = SPS_UCHAR; break;
  case PyArray_SBYTE: 
    type = SPS_CHAR; break;
  case PyArray_FLOAT: 
    type = SPS_FLOAT; break;
  case PyArray_DOUBLE: 
    type = SPS_DOUBLE; break;
  default:
    PyErr_SetString(SPSLUTError, "Input Array type not supported");
    return NULL;
  }
  
  data = src->data;
  cols = src->dimensions[1];  /*FIX THIS cols and rows are turned around */
  rows = src->dimensions[0];  /*###CHANGED - ALEXANDRE 24/07/2001*/


  
  r = SPS_PaletteArray (data, type, cols, rows, reduc, fastreduc, meth, gamma, 
			autoscale, mapmin, mapmax, Xservinfo, palette_code,
			&min, &max, &pcols, &prows, &palette, &pal_entries);
  if (r == 0) {
    PyErr_SetString(SPSLUTError, "Error while trying to calculate the image");
    return NULL;
  }
  if (!array_output){
   /*###CHANGED - ALEXANDRE 24/07/2001*/
  aux=new_pyimage(mode, (unsigned) pcols, (unsigned) prows, r);
  res = Py_BuildValue("(O(i,i)(d,d))",aux,pcols, prows, min, max);
  free(r);
  Py_DECREF(aux);

  /*###CHANGED - ALEXANDRE 28/06/2002*/  
  Py_DECREF(src);


  return res;
  }
  as_dim[0] = strlen(mode);
  as_dim[1] = prows * pcols;
  as_aux = (PyArrayObject*) PyArray_FromDims(2,as_dim, PyArray_UBYTE);
  if (as_aux == NULL){
      free(r);
      Py_DECREF(src);
      return NULL;  
  }
  as_pointer = (char *) as_aux -> data; 
  as_r = (char *) r;  
  memcpy(as_pointer, as_r, as_dim[0] * as_dim[1]);
  free(r);
  res = Py_BuildValue("(O(i,i)(d,d))",as_aux,pcols, prows, min, max);
  Py_DECREF(src);
  Py_DECREF(as_aux);
  return res;
}


static PyObject *spslut_transformarray(self, args)
     PyObject *self, *args;
     
{
  void *data;
  int type, cols, rows, reduc, fastreduc, meth, autoscale, mapmin=0, mapmax=255;
  int palette_code;
  double gamma, min, max;
  XServer_Info Xservinfo;
  void *palette;
  int prows, pcols, pal_entries;
  void *r/*, *res*/;
  char *mode;
  unsigned char *as_pointer, *as_r;
  int as_dim[3];
  PyArrayObject *src;
  PyObject *in_src;
  PyArrayObject *aux;

  if (!PyArg_ParseTuple(args, "O(ii)(id)sii(dd)|(ii)", &in_src, &reduc, 
			&fastreduc, &meth, &gamma, &mode, &palette_code,
			&autoscale, &min, &max,&mapmin, &mapmax)) 
		return NULL;
  
  if (strcmp(mode, "RGB") == 0) {
    Xservinfo.red_mask = 0x0000ff;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0xff0000;
    Xservinfo.pixel_size = 3;
    Xservinfo.byte_order = natbyteorder();
  } else if (strcmp(mode, "RGBX") == 0) {
    Xservinfo.red_mask = 0x0000ff;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0xff0000;
    Xservinfo.pixel_size = 4;
    Xservinfo.byte_order = natbyteorder();
  } 
  
  //###CHANGED - ALEXANDRE 11/03/2002 - Qt uses different order than Tkinter
  else if (strcmp(mode, "BGR") == 0) {
    Xservinfo.red_mask = 0xff0000;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0x0000ff;
    Xservinfo.pixel_size = 3;
    Xservinfo.byte_order = natbyteorder();
  } else if (strcmp(mode, "BGRX") == 0) {
    Xservinfo.red_mask = 0xff0000;
    Xservinfo.green_mask = 0x00ff00;
    Xservinfo.blue_mask = 0x0000ff;
    Xservinfo.pixel_size = 4;
    Xservinfo.byte_order = natbyteorder();
  } 

    else if (strcmp(mode, "L") == 0 || strcmp(mode, "P") == 0  ) {
    Xservinfo.pixel_size = 1;
    Xservinfo.byte_order = natbyteorder();
    //mapmin = 0;
    //mapmax = 255;
  } else {
    PyErr_SetString(SPSLUTError, "Mode must be RGB, RGBX, BGR, BGRX, L or P");
    return NULL;
  }

  if (!(src = (PyArrayObject*) PyArray_ContiguousFromObject(in_src, 
   			    PyArray_NOTYPE, 2, 2))) {
    PyErr_SetString(SPSLUTError, "Input Array could is not a 2x2 array");
    return NULL;
  }

  switch (src->descr->type_num) {
  case PyArray_UINT:
    type = SPS_ULONG; break;
  case PyArray_USHORT:
    type = SPS_USHORT; break;
  case PyArray_LONG: 
  case PyArray_INT: 
    type = SPS_LONG; break;
  case PyArray_SHORT: 
    type = SPS_SHORT; break;
  case PyArray_UBYTE: 
    type = SPS_UCHAR; break;
  case PyArray_SBYTE: 
    type = SPS_CHAR; break;
  case PyArray_FLOAT: 
    type = SPS_FLOAT; break;
  case PyArray_DOUBLE: 
    type = SPS_DOUBLE; break;
  default:
    PyErr_SetString(SPSLUTError, "Input Array type not supported");
    return NULL;
  }
  
  data = src->data;
  cols = src->dimensions[1];  /*FIX THIS cols and rows are turned around */
  rows = src->dimensions[0];  /*###CHANGED - ALEXANDRE 24/07/2001*/


  
  r = SPS_PaletteArray (data, type, cols, rows, reduc, fastreduc, meth, gamma, 
			autoscale, mapmin, mapmax, Xservinfo, palette_code,
			&min, &max, &pcols, &prows, &palette, &pal_entries);
  if (r == 0) {
    PyErr_SetString(SPSLUTError, "Error while trying to calculate the image");
    return NULL;
  }
   /*###CHANGED - ALEXANDRE 24/07/2001*/
  as_dim[0] = strlen(mode);
  as_dim[1] = prows * pcols;
/*  printf("dim[0] = %d dim[1] = %d \"%s\" \n",as_dim[0],as_dim[1],mode); */
  aux = (PyArrayObject*) PyArray_FromDims(2,as_dim, PyArray_CHAR);
  if (aux == NULL){
      free(r);
      Py_DECREF(src);
      return NULL;  
  }
  as_pointer = (char *) aux -> data; 
  as_r = (char *) r;  
  memcpy(as_pointer, as_r, as_dim[0] * as_dim[1]);
  free(r);
  Py_DECREF(src);
  PyArray_Return(aux);
}



/* The simple palette always returns 4 bytes per entry */
static PyObject *spslut_palette(self, args)
     PyObject *self, *args;
     
{
  int entries, palette_code;
  XServer_Info Xservinfo;
  void *r;
  char *mode;
  PyObject *res;

  if (!PyArg_ParseTuple(args, "ii", &entries, &palette_code)) 
    return NULL;
  
  mode = "RGBX";
  Xservinfo.red_mask = 0x0000ff;
  Xservinfo.green_mask = 0x00ff00;
  Xservinfo.blue_mask = 0xff0000;
  Xservinfo.pixel_size = 4;
  Xservinfo.byte_order = natbyteorder();

  r = SPS_SimplePalette ( 0, entries - 1, Xservinfo, palette_code); 
  
  if (r == 0) {
    PyErr_SetString(SPSLUTError, "Error calculating the palette");
    return NULL;
  }
  
  return new_pyimage(mode, 1, entries, r);
}


                                                                               
static PyMethodDef SPSLUTMethods[] = {
  { "transform", spslut_transform, METH_VARARGS},
  { "palette", spslut_palette, METH_VARARGS},
  { "transformArray", spslut_transformarray, METH_VARARGS},
  { NULL, NULL}
};
  
void initspslut()
{
  PyObject *d, *m;
  /* Create the module and add the functions */
  m = Py_InitModule ("spslut", SPSLUTMethods);
  
  /* Add some symbolic constants to the module */
  d = PyModule_GetDict(m);

  PyDict_SetItemString(d, "LINEAR", PyInt_FromLong(SPS_LINEAR));
  PyDict_SetItemString(d, "LOG", PyInt_FromLong(SPS_LOG));
  PyDict_SetItemString(d, "GAMMA", PyInt_FromLong(SPS_GAMMA));

  PyDict_SetItemString(d, "GREYSCALE", PyInt_FromLong(SPS_GREYSCALE));
  PyDict_SetItemString(d, "TEMP", PyInt_FromLong(SPS_TEMP));
  PyDict_SetItemString(d, "RED", PyInt_FromLong(SPS_RED));
  PyDict_SetItemString(d, "GREEN", PyInt_FromLong(SPS_GREEN));
  PyDict_SetItemString(d, "BLUE", PyInt_FromLong(SPS_BLUE));
  PyDict_SetItemString(d, "REVERSEGREY", PyInt_FromLong(SPS_REVERSEGREY));
  PyDict_SetItemString(d, "MANY", PyInt_FromLong(SPS_MANY));

  SPSLUTError = PyErr_NewException("spslut.error", NULL, NULL);
  PyDict_SetItemString(d, "error", SPSLUTError);
  import_array();
}

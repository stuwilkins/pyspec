/* 
 
 pyspec.ccd.princeton  
 (c) 2010 Stuart Wilkins <stuwilkins@mac.com>
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 
 $Id$
 
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include "princeton.h"

static PyObject* load_multiprinceton(PyObject *self, PyObject *args, PyObject *kwargs){
  static char *kwlist[] = {"filenames", "dark"};
  PyObject *array = NULL;
  PyObject *filenames = NULL;
  const char* darkimage;
  char *filename;
  int nimages, i;

  if(!PyArg_ParseTuppleAndKeywords(args, kwargs, "O|s", kwlist
				   &filenames, &darkimage)){
    return NULL;
  }

  if(!PyList_Check(filenames)){
    return NULL;
  }

  nimages = PyList_Size(filenames);
  for(i=0;i<nimages;i++){
    filename = PyString_AsString(PyList_GetItem(filenames));
    if(!filename){
      return NULL;
    }

    fprintf(stderr, "Filename %d = %s\n", i, filename);
  }

  return NULL;
}

static PyObject* load_princeton(PyObject *self, PyObject *args, PyObject *kwargs){
	PyObject *array				= NULL;
	//PyObject *result			= NULL;
	const char *filename		= NULL;
	const int normData			= 0;
	static char *kwlist[]		= {"filename", "norm", NULL};
	FILE *fp					= NULL;
	
	uint16_t dims[2]			= {0, 0};
	int32_t nImages				= 0;
	
	uint16_t *data				= NULL;
	
	npy_intp arraySize[3]		= {0, 0, 0};
	
	//double *arrayPtr			= NULL;
	uint16_t *arrayPtr			= NULL;
	//uint16_t *dataPtr			= NULL;
	//int i;
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "s|i", kwlist,
									 &filename, &normData)){
        return NULL;
	}
	
	if((fp = fopen(filename, "rb")) == NULL){
		PyErr_SetString(PyExc_IOError, "Unable to open file.");
		return NULL;
	}
	
	/* 
	 Read the size of image and the number of images
	 */
	
	fseek(fp, XDIM_OFFSET, SEEK_SET);
	fread(&dims[0], sizeof(uint16_t), 1, fp);
	fseek(fp, YDIM_OFFSET, SEEK_SET);
	fread(&dims[1], sizeof(uint16_t), 1, fp);
	fseek(fp, NFRAMES_OFFSET, SEEK_SET);
	fread(&nImages, sizeof(int32_t), 1, fp);
	
	/*
	
	data = (uint16_t*)malloc(dims[0] * dims[1] * nImages * sizeof(uint16_t));
	if(!data){
		PyErr_SetString(PyExc_MemoryError, "Unable to allocate memory.");
		return NULL;
	}
	
	fseek(fp, DATA_OFFSET, SEEK_SET);
	if(!fread(data, sizeof(uint16_t), dims[0] * dims[1] * nImages, fp)){
		PyErr_SetString(PyExc_IOError, "Unable to read data.");
		free(data);
		return NULL;
	}
	*/
	
	arraySize[0] = nImages;
	arraySize[1] = dims[1];
	arraySize[2] = dims[0];
	array = PyArray_SimpleNew(3, arraySize, NPY_USHORT);
	
	arrayPtr = (uint16_t *)PyArray_DATA(array);
	
	fseek(fp, DATA_OFFSET, SEEK_SET);
	if(!fread(arrayPtr, sizeof(uint16_t), dims[0] * dims[1] * nImages, fp)){
		PyErr_SetString(PyExc_IOError, "Unable to read data.");
		free(data);
		return NULL;
	}
	
	/*
	dataPtr = data;
	
	for(i=0;i<(dims[0] * dims[1] * nImages);i++){
		*(arrayPtr++) = (double)*(dataPtr++);
	}
	
	if(!data){
		free(data);
	}
	*/
	return array;
}

PyMODINIT_FUNC initprinceton(void)  {
	(void) Py_InitModule3("princeton", _princetonMethods, _princetonDoc);
	import_array();  // Must be present for NumPy.  Called first after above line.
}

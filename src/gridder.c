/* 
 
 pyspec.gridder  
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
#include "gridder.h"

static PyObject* gridder_3D(PyObject *self, PyObject *args, PyObject *kwargs){
	PyObject *gridout = NULL, *Nout = NULL;
	PyObject *gridI = NULL;
	PyObject *_I;
	
	static char *kwlist[] = { "data", "xrange", "yrange", "zrange", "norm", NULL };

	npy_intp data_size;
	npy_intp dims[3];
	
	double grid_start[3];
	double grid_stop[3];
	int grid_nsteps[3];
	int norm_data = 0;

	unsigned long n_outside;
	
	if(!PyArg_ParseTupleAndKeywords(args, kwargs, "O(ddd)(ddd)(iii)|i", kwlist,
					&_I, 
					&grid_start[0], &grid_start[1], &grid_start[2],
					&grid_stop[0], &grid_stop[1], &grid_stop[2],
					&grid_nsteps[0], &grid_nsteps[1], &grid_nsteps[2],
					&norm_data)){
	  return NULL;
	}	
	
	gridI = PyArray_FROMANY(_I, NPY_DOUBLE, 0, 0, NPY_IN_ARRAY);
	if(!gridI){
		return NULL;
	}
	
	data_size = PyArray_DIM(gridI, 0);
	
	dims[0] = grid_nsteps[0];
	dims[1] = grid_nsteps[1];
	dims[2] = grid_nsteps[2];
	gridout = PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);
	if(!gridout){
		return NULL;
	}
	Nout = PyArray_ZEROS(3, dims, NPY_ULONG, 0);
	if(!Nout){
		return NULL;
	}
	
	n_outside = c_grid3d(PyArray_DATA(gridout), PyArray_DATA(Nout), 
			     PyArray_DATA(gridI),
			     grid_start, grid_stop, data_size, grid_nsteps, norm_data);
	
	return Py_BuildValue("OOl", gridout, Nout, n_outside); 
}

unsigned long c_grid3d(double *dout, unsigned long *nout, double *data, 
		       double *grid_start, double *grid_stop, int max_data, 
		       int *n_grid, int norm_data){
	int i;
	unsigned long *nout_ptr;
	double *dout_ptr;
	double *data_ptr;
	
	double pos_double[3];
	double grid_len[3], grid_step[3];
	int grid_pos[3];
	int pos;
	unsigned long n_outside = 0;
	
	fprintf(stderr, "Gridding in 3D : grid pts = %i x %i x %i, data pts = %i\n", 
		n_grid[0], n_grid[1], n_grid[2], max_data);
	
	dout_ptr = dout;
	nout_ptr = nout;
	
	data_ptr = data;
	for(i = 0;i < 3; i++){
		grid_len[i] = grid_stop[i] - grid_start[i];
		grid_step[i] = grid_len[i] / (n_grid[i]);
	}
	
	for(i = 0; i < max_data ; i++){
		pos_double[0] = (*data_ptr - grid_start[0]) / grid_len[0];
		data_ptr++;
		pos_double[1] = (*data_ptr - grid_start[1]) / grid_len[1];
		data_ptr++;
		pos_double[2] = (*data_ptr - grid_start[2]) / grid_len[2];
		if((pos_double[0] >= 0) && (pos_double[0] < 1) && 
		   (pos_double[1] >= 0) && (pos_double[1] < 1) &&
		   (pos_double[2] >= 0) && (pos_double[2] < 1)){
			data_ptr++;	
			grid_pos[0] = (int)(pos_double[0] * n_grid[0]);
			grid_pos[1] = (int)(pos_double[1] * n_grid[1]);
			grid_pos[2] = (int)(pos_double[2] * n_grid[2]);
			pos =  grid_pos[0] * (n_grid[1] * n_grid[2]);
			pos += grid_pos[1] * n_grid[2];
			pos += grid_pos[2];
			dout[pos] = dout[pos] + *data_ptr;
			nout[pos] = nout[pos] + 1;
			data_ptr++;
		} else {
		  //fprintf(stderr, "Data point (%lf,%lf,%lf) is out of range\n", pos_double[0], pos_double[1], pos_double[2]);
		  n_outside++;
		  data_ptr+=2;
		}
	}

	if(norm_data){
	  for(i = 0; i < (n_grid[0] * n_grid[1] * n_grid[2]); i++){
	    if(nout[i] > 0){
	      dout[i] = dout[i] / nout[i];
	    } else {
	      dout[i] = 0.0;
	    }
	  }
	}
	
	return n_outside;
}

PyMODINIT_FUNC initgridder(void)  {
	(void) Py_InitModule3("gridder", _gridderMethods, _gridderDoc);
	import_array();  // Must be present for NumPy.  Called first after above line.
}


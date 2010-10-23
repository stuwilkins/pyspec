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

static PyObject* gridder_3D(PyObject *self, PyObject *args){
	PyObject *gridout = NULL, *Nout = NULL;
	PyObject *gridI = NULL;
	PyObject *_I;
	
	npy_intp data_size;
	npy_intp dims[3];
	
	double grid_start[3];
	double grid_stop[3];
	int grid_nsteps[3];

	unsigned long n_outside;
	
	if(!PyArg_ParseTuple(args, "Oddiddiddi", 
			     &_I, 
			     &grid_start[0], &grid_stop[0], &grid_nsteps[0],
			     &grid_start[1], &grid_stop[1], &grid_nsteps[1],
			     &grid_start[2], &grid_stop[2], &grid_nsteps[2])){
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
			     grid_start, grid_stop, data_size, grid_nsteps);
	
	return Py_BuildValue("OOl", gridout, Nout, n_outside); 
}

static PyObject* gridder_2D(PyObject *self, PyObject *args){
	PyObject *gridout = NULL, *Nout = NULL;
	PyObject *gridI = NULL;
	PyObject *_I;
	
	npy_intp data_size;
	npy_intp dims[2];
	
	double grid_start[2];
	double grid_stop[2];
	int grid_nsteps[2];
	
	if(!PyArg_ParseTuple(args, "Oddiddi", 
						 &_I, 
						 &grid_start[0], &grid_stop[0], &grid_nsteps[0],
						 &grid_start[1], &grid_stop[1], &grid_nsteps[1])){
		return NULL;
	}	
	
	gridI = PyArray_FROMANY(_I, NPY_DOUBLE, 0, 0, NPY_IN_ARRAY);
	if(!gridI){
		fprintf(stderr, "gridI\n");
		return NULL;
	}
	
	data_size = PyArray_DIM(gridI, 0);
	
	dims[0] = grid_nsteps[0];
	dims[1] = grid_nsteps[1];
	gridout = PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
	if(!gridout){
		return NULL;
	}
	Nout = PyArray_ZEROS(2, dims, NPY_ULONG, 0);
	if(!Nout){
		return NULL;
	}
	
	c_grid2d(PyArray_DATA(gridout), PyArray_DATA(Nout), PyArray_DATA(gridI),
			 grid_start, grid_stop, data_size, grid_nsteps);
	
	return Py_BuildValue("OO", gridout, Nout); 
}

static PyObject* gridder_1D(PyObject *self, PyObject *args){
	PyObject *gridout = NULL, *Nout = NULL;
	PyObject *gridI = NULL;
	PyObject *_I;
	
	npy_intp data_size;
	npy_intp dims[1];
	
	double grid_start;
	double grid_stop;
	int grid_nsteps;
	
	if(!PyArg_ParseTuple(args, "Oddi", &_I, &grid_start, &grid_stop, &grid_nsteps)){
		return NULL;
	}	
	
	gridI = PyArray_FROMANY(_I, NPY_DOUBLE, 0, 0, NPY_IN_ARRAY);
	if(!gridI){
		return NULL;
	}
	
	data_size = PyArray_DIM(gridI, 0);
	
	dims[0] = grid_nsteps;
	gridout = PyArray_ZEROS(1, dims, NPY_DOUBLE, 0);
	if(!gridout){
		return NULL;
	}
	Nout = PyArray_ZEROS(1, dims, NPY_ULONG, 0);
	if(!Nout){
		return NULL;
	}
	
	c_grid1d(PyArray_DATA(gridout), PyArray_DATA(Nout), PyArray_DATA(gridI),
			 grid_start, grid_stop, data_size, grid_nsteps);
	
	return Py_BuildValue("OO", gridout, Nout); 
}

int c_grid2d(double *dout, unsigned long *nout, double *data, double *grid_start, double *grid_stop, int max_data, int *n_grid){
	int i;
	unsigned long *nout_ptr;
	double *dout_ptr;
	double *data_ptr;
	
	double pos_double[2];
	double grid_len[2], grid_step[2];
	int grid_pos[2];
	int pos;
	
	fprintf(stderr, "Gridding in 2D : grid pts = %i x %i, data pts = %i\n", n_grid[0], n_grid[1], max_data);
	
	dout_ptr = dout;
	nout_ptr = nout;
	
	data_ptr = data;
	for(i = 0;i < 2; i++){
		grid_len[i] = grid_stop[i] - grid_start[i];
		grid_step[i] = grid_len[i] / (n_grid[i]);
	}
	
	for(i = 0; i < max_data ; i++){
		pos_double[0] = (*data_ptr - grid_start[0]) / grid_len[0];
		data_ptr++;
		pos_double[1] = (*data_ptr - grid_start[1]) / grid_len[1];
		if((pos_double[0] >= 0) && (pos_double[0] < 1) && 
		   (pos_double[1] >= 0) && (pos_double[1] < 1)){
			data_ptr += 2;	
			grid_pos[0] = (int)(pos_double[0] * n_grid[0]);
			grid_pos[1] = (int)(pos_double[1] * n_grid[1]);
			pos = grid_pos[0] * n_grid[1] + grid_pos[1];
			dout[pos] = dout[pos] + *data_ptr;
			nout[pos] = nout[pos] + 1;
			data_ptr++;
		} else {
		  //fprintf(stderr, "Data point x = %lf (%lf,%lf) is out of range\n", *data_ptr, pos_double[0], pos_double[1]);
		  data_ptr+=3;
		}
	}
	
	return 0;
}

unsigned long c_grid3d(double *dout, unsigned long *nout, double *data, double *grid_start, double *grid_stop, int max_data, int *n_grid){
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
	
	return n_outside;
}

int c_grid1d(double *dout, unsigned long *nout, double *data, double grid_start, double grid_stop, int max_data, int n_grid){
	int i;
	unsigned long *nout_ptr;
	double *dout_ptr;
	double *data_ptr;
	
	double posx_double;
	double grid_len, grid_step;
	int grid_pos;
	
	fprintf(stderr, "Gridding in 1D : grid pts = %i, data pts = %i\n", n_grid, max_data);
	
	dout_ptr = dout;
	nout_ptr = nout;
		
	data_ptr = data;
	grid_len = grid_stop - grid_start;
	grid_step = grid_len / (n_grid);
	
	for(i = 0; i < max_data ; i++){
		posx_double = (*data_ptr - grid_start) / grid_len;
		if((posx_double >= 0) && (posx_double < 1)){
			data_ptr += 3;
			grid_pos = (int)(posx_double * n_grid);
			//fprintf(stderr, "%i\n", grid_pos);
			dout[grid_pos] = dout[grid_pos] + *data_ptr;
			nout[grid_pos] = nout[grid_pos] + 1;
			data_ptr++;
		} else {
		  //fprintf(stderr, "Data point x = %lf (%lf) is out of range\n", *data_ptr, posx_double);
		  data_ptr+=4;
		}
	}
	
	return 0;
}

PyMODINIT_FUNC initgridder(void)  {
	(void) Py_InitModule3("gridder", _gridderMethods, _gridderDoc);
	import_array();  // Must be present for NumPy.  Called first after above line.
}


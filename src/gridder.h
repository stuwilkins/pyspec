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

#ifndef __GRIDDER_H
#define __GRIDDER_H

int c_grid1d(double *dout, unsigned long *nout, double *data, double grid_start, double grid_stop, int max_data, int n_grid);
int c_grid2d(double *dout, unsigned long *nout, double *data, double *grid_start, double *grid_stop, int max_data, int *n_grid);
unsigned long c_grid3d(double *dout, unsigned long *nout, double *data, double *grid_start, double *grid_stop, int max_data, int *n_grid);

static PyObject* gridder_3D(PyObject *self, PyObject *args);
static PyObject* gridder_2D(PyObject *self, PyObject *args);
static PyObject* gridder_1D(PyObject *self, PyObject *args);

static char *_gridderDoc = \
"Python functions to perform gridding (binning) of experimental data.\n\n" 
"Members:\n\n"
"  grid1d = 4 element tuple representing the default\n"
"                        opts into dder an ddif.\n"
"  levmar.INIT_MU      = Initial value for mu\n"
"  levmar.STOP_THRESH  = Stopping threshold\n"
"  levmar.DIFF_DELTA   = Differential delta\n";

static PyMethodDef _gridderMethods[] = {
	{"grid1d", gridder_1D, METH_VARARGS, 
		""},
	{"grid2d", gridder_2D, METH_VARARGS, 
		""},
	{"grid3d", gridder_3D, METH_VARARGS, 
		""},
	{NULL, NULL, 0, NULL}     /* Sentinel - marks the end of this structure */
};

#endif


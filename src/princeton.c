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
 
 $Id: gridder.c 52 2010-10-20 21:59:42Z stuwilkins $
 
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#inlucde "princeton.h"

static PyObject* load_princeton(PyObject *self, PyObject *args, PyObject *kwargs){

}

int read_image(FILE *fp, unsigned short *xdim, unsigned short *ydim, 


PyMODINIT_FUNC initprinceton(void)  {
	(void) Py_InitModule3("princeton", _princetonMethods, _princetonDoc);
	import_array();  // Must be present for NumPy.  Called first after above line.
}
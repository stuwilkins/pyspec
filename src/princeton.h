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

#include <stdint.h>

#ifndef __PRINCETON_H
#define __PRINCETON_H

#define HDRNAMEMAX		120 
#define USERINFOMAX		1000 
#define COMMENTMAX		80
#define LABELMAX		16
#define FILEVERMAX		16
#define DATEMAX			10 
#define ROIMAX			10
#define TIMEMAX			7 

#define DATA_OFFSET		4100
#define XDIM_OFFSET		42
#define YDIM_OFFSET		656
#define TYPE_OFFSET		108
#define NFRAMES_OFFSET	1446

static PyObject* load_princeton(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject* load_multiprinceton(PyObject *self, PyObject *args, PyObject *kwargs);

static char *_princetonDoc = \
"Python functions to read Princeton Instruments SPE files.\n\n";

static PyMethodDef _princetonMethods[] = {
	{"readSPE", (PyCFunction)load_princeton, METH_VARARGS | METH_KEYWORDS, 
	 ""},
	{"readMultiSPE", (PyCFunction)load_multiprinceton, METH_VARARGS | METH_KEYWORDS,
	 ""},
	{NULL, NULL, 0, NULL}     /* Sentinel - marks the end of this structure */
};

#endif

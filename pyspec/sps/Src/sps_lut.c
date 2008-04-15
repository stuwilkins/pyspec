#/*##########################################################################
# Copyright (C) 2004-2006 European Synchrotron Radiation Facility
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
#include <stdio.h>
/* AS + AM */
//#ifdef __linux__
#if 0
	#warning "Assuming LINUX and using the functions lrint() and lrintf ()."
	#define	_ISOC9X_SOURCE	1
	#define _ISOC99_SOURCE	1

	#define	__USE_ISOC9X	1
	#define	__USE_ISOC99	1

	#include	<math.h>
	
/*#elif (defined (WIN32) || defined (_WIN32))
	#include	<math.h>*/
#elif 0

	/*	Win32 doesn't seem to have these functions. 
	**	Therefore implement inline versions of these functions here.
	*/
	
	__inline long int 
	lrint (double flt) ;
	{	int intgr;

		_asm
		{	fld flt
			fistp intgr
			} ;
			
		return intgr ;
	} 
	
	__inline long int 
	lrintf (float flt)
	{	int intgr;

		_asm
		{	fld flt
			fistp intgr
			} ;
			
		return intgr ;
	}

#else
	#include	<math.h>

	#define	lrint(dbl)		((int)(dbl))
	#define	lrintf(flt)		((int)(flt))

#endif

#include <limits.h>
/* #include <malloc.h> */
#include <sps.h>
#include <sps_lut.h>
#include <blissmalloc.h>

#ifndef HPUX
#define log10f log10
#define powf pow
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e+308
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.40282347e+38
#endif


/* machine dependent */
typedef union {
  struct {
    unsigned char dummy;
    unsigned char R;
    unsigned char G;
    unsigned char B;
  } c;
  unsigned int p;
} RGB24bits;

typedef union {
  struct {
    unsigned char b1;
    unsigned char b2;
    unsigned char b3;
    unsigned char b4;
  } c;
  unsigned int p;
} swaptype;

int SPS_Size_VLUT (int t)
{
  switch (t) {
  case SPS_USHORT: return(sizeof(unsigned short));
  case SPS_ULONG:  return(sizeof(unsigned int));
  case SPS_SHORT:  return(sizeof(short));
  case SPS_LONG:   return(sizeof(int));
  case SPS_UCHAR:  return(sizeof(unsigned char));
  case SPS_CHAR:   return(sizeof(char));
  case SPS_STRING: return(sizeof(char));
  case SPS_DOUBLE: return(sizeof(double));
  case SPS_FLOAT:  return(sizeof(float));
  default:        return(0);
  }
}

void *CreatePalette( int type, int meth, double min, double max, double gamma,
		     int mapmin, int mapmax, 
                     XServer_Info Xservinfo, int palette_code);

unsigned char *SPS_SimplePalette ( int min, int max, 
                                   XServer_Info Xservinfo, int palette_code) 
{
  int type = SPS_USHORT, meth = SPS_LINEAR;
  int mapmin = 0, mapmax = 0; /* Is not used in this case USHORT LINEAR*/
  double dmin = min, dmax = max, gamma = 0.0;
  if (Xservinfo.pixel_size == 1) 
    Xservinfo.pixel_size = 3;
  return CreatePalette(type, meth, dmin, dmax, gamma, mapmin, mapmax, 
		       Xservinfo, palette_code);
  

}

void *SPS_PaletteArray (void *data, int type, int cols, int rows, 
          int reduc, int fastreduc, int meth, double gamma, int autoscale, 
	  int mapmin, int mapmax, 
          XServer_Info Xservinfo, int palette_code,
	  double *min, double *max, int *pcols, int *prows,
	  void **pal_return, int *pal_entries)
{
  int calcminmax;
  void *ndata, *Xdata;
  void *palette = NULL;
  double use_min, use_max;
  double minplus = 0;

  *pal_entries = 0;
  *pal_return = NULL;

  if (Xservinfo.pixel_size != 1) {
    mapmin = 0;
    mapmax = 0xffff; /* do we really need more then 256 colors - maybe 0xff */
  }

  /* Calculate the min max minplus of the data only if necessary       */
  /* Calc minplus is very fast if min alread > 0 - always calc minplus */
  calcminmax = (autoscale ? 1 : 0 ) | ((meth != SPS_LINEAR) ? 2 : 0);
 
  if (calcminmax)
    SPS_FindMinMax(data, type, cols, rows, min, max, &minplus, calcminmax); 

  /* Reduce the data with reduction factor - nothing done for reduc == 1 */
  ndata = SPS_ReduceData(data, type, cols, rows, reduc, pcols, prows,
			 fastreduc);
  if (ndata == NULL) 
    return NULL;

  if (meth == SPS_LINEAR) {
    use_min = *min;
    use_max = *max;
  } else if (type == SPS_USHORT || type == SPS_SHORT || type == SPS_CHAR || 
	     type == SPS_UCHAR) {
    use_min = *min; /* Check -  we treat signed types as unsigned ???     */
    use_max = *max; /* Does the palette look like the user expects ???    */
  } else {
    if (minplus == 0) {
      use_min = use_max = 1; /* No value above 0 */
    } else {
      use_min = minplus ; /* Same as min if min > 0 */
      use_max = ( *max > minplus ) ? *max : minplus;
    }
  }

  /* Create the palette if we do not have a hardware palette */
  palette = CreatePalette( type, meth, use_min, use_max, gamma,
                           mapmin, mapmax, Xservinfo, palette_code);
   
  /* Produce an array with data between mapmin and mapmax for reference into
     the palette */
  Xdata = SPS_MapData(ndata, type, meth, *pcols, *prows, use_min, use_max, 
		      gamma, mapmin, mapmax, Xservinfo.pixel_size, palette);
  if (Xdata == NULL) 
    return NULL;
  
  if (ndata != data)
    free (ndata);

  if (Xservinfo.pixel_size != 1) {
    if (type == SPS_USHORT || type == SPS_SHORT || type == SPS_CHAR || 
	type == SPS_UCHAR) {
      *pal_return = (void *) (((unsigned char *)palette) +
			      (int )(Xservinfo.pixel_size * *min));
      *pal_entries = *max - *min + 1;
    } else {
      *pal_return = (void *) (((unsigned char *)palette) + 
			      (int )(Xservinfo.pixel_size * mapmin));
      *pal_entries = mapmax - mapmin + 1;
    }
  }
  
  if (meth != SPS_LINEAR) {
    *min = minplus;
  }

  return Xdata;
}

#define FINDMINMAX(ty, maxval) \
{\
 register ty *c1=(ty *)data;\
 register int i;\
 register ty mmin, mmax, mminplus;\
 int size = cols*rows; \
 mmax = mmin = *c1;\
 mminplus = maxval;\
 if (dominmax && dominplus) { \
   for (i=size;i;i--,c1++) {\
     if (*c1 < mmin)\
       mmin = *c1;\
     if (*c1 > mmax)\
       mmax = *c1;\
     if ((*c1 < mminplus) && (*c1 > 0))\
       mminplus = *c1;\
   }\
 } else if (dominmax) {\
   for (i=size;i;i--,c1++) {\
     if (*c1 < mmin)\
       mmin = *c1;\
     if (*c1 > mmax)\
       mmax = *c1;\
   }\
 } else if (dominplus) {\
   if (mmin > 0) \
     mminplus = mmin;\
   else \
     for (i=size;i;i--,c1++) {\
       if ((*c1 < mminplus) && (*c1 > 0))\
         mminplus = *c1;\
   }\
 }\
 dmin = (double)mmin;\
 dmax = (double)mmax;\
 dminplus = (double)mminplus;\
}

void SPS_FindMinMax(void *data, int type, int cols, int rows,
		  double *min, double *max, double *minplus, int flag)
{
 int dominmax = flag & 1;
 int dominplus = flag & 2;
 double dmin,dmax,dminplus;
 
 switch (type) {
   case SPS_DOUBLE :
     FINDMINMAX(double, DBL_MAX);
     break;
   case SPS_FLOAT :
     FINDMINMAX(float, FLT_MAX);
     break;
   case SPS_LONG :
     FINDMINMAX(int, LONG_MAX);
     break;
   case SPS_ULONG :
     FINDMINMAX(unsigned int, ULONG_MAX);
     break;
   case SPS_SHORT :
     FINDMINMAX(short, SHRT_MAX);
     break;
   case SPS_USHORT :
     FINDMINMAX(unsigned short, USHRT_MAX);
     break;
   case SPS_CHAR :
     FINDMINMAX(char, SCHAR_MAX);
     break;
   case SPS_UCHAR :
     FINDMINMAX(unsigned char, UCHAR_MAX);
     break;
 }

 if (dominmax) {
   *min = dmin;
   *max = dmax;
 }

 if (dominplus)
   *minplus = dminplus;
 
}

#define CALCDATA(ty, mty, mapty, logfct, powfct)\
{\
 ty vmin, vmax;\
 register ty val;\
 register mapty *Xptr;\
 register ty *ptr;\
 register mty Au = A;\
 register mty Bu = B;\
 register mapty cmin = (mapty) mapmin, cmax = (mapty) mapmax;\
 int size;\
 vmin = (ty)Xmin;\
 vmax = (ty)Xmax;\
 if (Xmin > Xmax) {\
   vmin = (ty)Xmax;\
   vmax = (ty)Xmin;\
 }\
 size = cols*rows;\
 ptr = (ty *)data;\
 if (meth == SPS_LINEAR) {\
   if (mapbytes == 1) {\
     register mapty *Xend = (mapty *)Xdata+size; \
     for (Xptr=(mapty *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax)\
         *Xptr = cmax ;\
       else if  (val > vmin) {\
         *Xptr = (mapty)(Au * (mty) val + Bu) ;\
       }else\
         *Xptr = cmin ;\
     }\
   } else if (mapbytes == 3) {\
     register unsigned char * Xend = (unsigned char *)Xdata + (3 * size); \
     register unsigned int *palette = (unsigned int *)pal;\
     register unsigned char *Xptr;\
     register RGB24bits pval;\
     for (Xptr=(unsigned char *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax) {\
         pval.p = *(palette + mapmax) ;\
       } else if  (val > vmin) {\
         pval.p = *(palette + lrint(Au * (mty) val + Bu)) ;\
       } else {\
         pval.p = *(palette + mapmin) ;\
       } \
       *Xptr = pval.c.R;  Xptr++;\
       *Xptr = pval.c.G;  Xptr++;\
       *Xptr = pval.c.B;\
     } \
   } else {\
     register mapty *Xend = (mapty *)Xdata+size; \
     register mapty *palette = (mapty *)pal;\
     for (Xptr=(mapty *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax)\
         *Xptr = *(palette + mapmax) ;\
       else if  (val > vmin)\
         *Xptr = *(palette + lrint(Au * (mty) val + Bu)) ;\
       else\
         *Xptr = *(palette + mapmin) ;\
     }\
   }\
 }\
 else if (meth == SPS_LOG) {\
   if (mapbytes == 1) {\
     register mapty *Xend = (mapty *)Xdata+size;\
     for (Xptr=(mapty *)Xdata;Xptr != Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if  (val >= vmax)\
         *Xptr = cmax;\
       else if (val > vmin)\
         *Xptr = (mapty)(A * logfct((mty)val) + B);\
       else\
         *Xptr = cmin;\
     }\
   } else if (mapbytes == 3) {\
     register unsigned char * Xend = (unsigned char *)Xdata + (3 * size); \
     register unsigned int *palette = (unsigned int *)pal;\
     register unsigned char *Xptr;\
     register RGB24bits pval;\
     for (Xptr=(unsigned char *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax) {\
         pval.p = *(palette + mapmax) ;\
       } else if  (val > vmin) {\
         pval.p = *(palette + lrint(A * logfct((mty) val) + B)) ;\
       } else {\
         pval.p = *(palette + mapmin) ;\
       } \
       *Xptr = pval.c.R;  Xptr++;\
       *Xptr = pval.c.G;  Xptr++;\
       *Xptr = pval.c.B;\
     } \
   } else {\
     register mapty *Xend = (mapty *)Xdata+size; \
     register mapty *palette = (mapty *)pal;\
     for (Xptr=(mapty *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax)\
         *Xptr = *(palette + mapmax) ;\
       else if  (val > vmin) {\
         *Xptr = *(palette + lrint(A * logfct((mty) val) + B)) ;\
       } else\
         *Xptr = *(palette + mapmin) ;\
     }\
   }\
 }\
 else if (meth == SPS_GAMMA) {\
   if (mapbytes == 1) {\
     register mapty *Xend = (mapty *)Xdata+size;\
     for (Xptr=(mapty *)Xdata;Xptr != Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax)\
         *Xptr = cmax;\
       else if  (val > vmin)\
         *Xptr = (mapty)(A *powfct((mty)val,(mty)gamma)+B);\
       else\
         *Xptr = cmin;\
     }\
   } else if (mapbytes == 3) {\
     register unsigned char * Xend = (unsigned char *)Xdata + (3 * size); \
     register unsigned int *palette = (unsigned int *)pal;\
     register unsigned char *Xptr;\
     register RGB24bits pval;\
     for (Xptr=(unsigned char *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax) {\
         pval.p = *(palette + mapmax) ;\
       } else if  (val > vmin) {\
         pval.p = *(palette + lrint(A * powfct((mty) val, (mty) gamma) + B)) ;\
       } else {\
         pval.p = *(palette + mapmin) ;\
       } \
       *Xptr = pval.c.R;  Xptr++;\
       *Xptr = pval.c.G;  Xptr++;\
       *Xptr = pval.c.B;\
     } \
   } else {\
     register mapty *Xend = (mapty *)Xdata+size;\
     register mapty *palette = (mapty *)pal;\
     for (Xptr=(mapty *)Xdata;Xptr != Xend;Xptr++,ptr++) {\
       val=*ptr;\
       if (val >= vmax)\
         *Xptr = *(palette + mapmax) ;\
       else if  (val > vmin)\
         *Xptr = *(palette + lrint (A * powfct((mty)val,(mty)gamma)+B)) ;\
       else\
         *Xptr = *(palette + mapmin) ;\
     }\
   }\
 }\
}

#define CALCDATA_NOMAP(ty, mapty)\
{\
  register ty val;\
  register ty *ptr = (ty *)data;\
  int size = cols*rows;\
  if (mapbytes == 3) {\
    register unsigned char * Xend = (unsigned char *)Xdata + (3 * size); \
    register unsigned int *palette = (unsigned int *)pal;\
    register unsigned char *Xptr;\
    register RGB24bits pval;\
    for (Xptr=(unsigned char *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
      val=*ptr;\
      pval.p = *(palette + val);\
      *Xptr = pval.c.R;  Xptr++;\
      *Xptr = pval.c.G;  Xptr++;\
      *Xptr = pval.c.B;\
    } \
  } else {\
    register mapty *Xend = (mapty *)Xdata+size; \
    register mapty *palette = (mapty *)pal;\
    register mapty *Xptr;\
    for (Xptr=(mapty *)Xdata;Xptr!=Xend;Xptr++,ptr++) {\
      val=*ptr;\
      *Xptr = *(palette + val);\
    }\
  }\
}

#define CALCPREDATA(ty, mty, mapty, logfct, powfct, premin , premax)\
{\
 ty vmin, vmax;\
 ty val;\
 mapty *Xptr;\
 mapty *fb;\
 int lval;\
 register ty *ptr;\
 register mty Au = A;\
 register mty Bu = B;\
 mapty cmin=mapmin, cmax=mapmax;\
 register mapty *Xend; \
 vmin = (ty)Xmin;\
 vmax = (ty)Xmax;\
 if (Xmin > Xmax) {\
   vmin = (ty)Xmax;\
   vmax = (ty)Xmin;\
 }\
 fb=Xptr=(mapty*)malloc (sizeof(mapty)*(premax-premin+1));\
 Xend = Xptr + (premax-premin+1);\
 for (lval=premin;lval<=(int)vmin;lval++) \
   *Xptr++=cmin;\
   if (meth == SPS_LINEAR) {\
   for (;lval<(int)vmax;lval++) \
     *Xptr++ = (mapty)(Au * (mty) lval + Bu);\
 } else if (meth == SPS_LOG) {\
   for (;lval<(int)vmax;lval++) \
     *Xptr++ = (mapty)(Au * logfct((mty)lval) + Bu);\
 } else if (meth == SPS_GAMMA) {\
   for (;lval<(int)vmax;lval++) \
     *Xptr++ = (mapty)(Au * powfct((mty)lval,(mty)gamma) + Bu);\
 }\
 for (;Xptr < Xend;) \
   *Xptr++ = cmax;\
\
 ptr = (ty *)data;\
 Xend = (mapty *)Xdata+cols*rows; \
 for (Xptr=(mapty *)Xdata;Xptr!=Xend;Xptr++,ptr++) \
   *Xptr=fb[*ptr];\
 free (fb);\
}

#define NOSTEPS 150
#define BINSTEPS 256

#define FASTLOG(ty, mapty)\
{ \
  int i;\
  ty xi,x[256],*ptr;\
  double ydbl;\
  ty *xptr,*xmiddle=x+128;\
  double ymin=mapmin,ymax=mapmax;\
  unsigned char idx, *Xend = (mapty *)Xdata+cols*rows,*Xptr;\
  int nosteps=mapmax - mapmin + 1;\
  int v64=64,v32=32,v16=16,v8=8,v4=4,v2=2,v1=1,v01=1,v02=0;\
  \
  for (i=0;i<nosteps;i++) {\
    ydbl = ymin + (double) i ;\
    x[i] = (ty)pow((double)10,(ydbl-B)/A);\
  }\
\
  for (i=nosteps;i<BINSTEPS;i++) \
    x[i] = x[NOSTEPS-1];\
\
  if (x[0] > x[nosteps-1]) {\
    v64=-64;v32=-32;v16=-16;v8=-8;v4=-4;v2=-2;v1=-1;v01=0;v02=1;\
  }\
\
  ptr = (ty *)data;\
  for (Xptr=(mapty *)Xdata;Xptr != Xend;Xptr++,ptr++) {\
    xptr = xmiddle;\
    xi = *ptr;\
    if (xi > *xptr)\
      xptr += v64;\
    else\
      xptr -= v64;\
	\
    if (xi > *xptr)\
      xptr += v32;\
    else\
      xptr -= v32;\
    \
    if (xi > *xptr)\
      xptr += v16;\
    else\
      xptr -= v16;\
    \
    if (xi > *xptr)\
      xptr += v8;\
    else\
      xptr -= v8;\
    \
    if (xi > *xptr)\
      xptr += v4;\
    else\
      xptr -= v4;\
    \
    if (xi > *xptr)\
      xptr += v2;\
    else\
      xptr -= v2;\
    \
    if (xi > *xptr)\
      xptr += v1;\
    else\
      xptr -= v1;\
    \
    if (xi < *xptr)\
      xptr -= v01;\
    else     \
      xptr -= v02;\
\
    idx = (unsigned char)(xptr - x);\
    if (idx > nosteps - 1)\
      *Xptr = (mapty)mapmax;\
    else\
      *Xptr = idx + mapmin;\
  }\
}


unsigned char *SPS_MapData(void *data, int type, int meth, int cols, int rows,
			   double Xmin, double Xmax, double gamma,
  			   int mapmin, int mapmax, int mapbytes, void *pal)
{
 
 double A, B, lmin, lmax;
 void *Xdata;
 int databytes ;

 /* mapbyte == 1 means that there is no palette - we output char */
 databytes = mapbytes ? mapbytes : 1;

 Xdata = (void *)malloc(databytes * cols * rows);
 if (Xdata == NULL) {
   fprintf(stderr, "Malloc Error in CalcData (size = %d), Exit\n",
	   cols*rows);
   return((unsigned char *)NULL);
 }
 
 if ((Xmax-Xmin) != 0) {
   if (meth == SPS_LINEAR) {
     lmin = Xmin;
     lmax = Xmax;
   }
   if (meth == SPS_LOG) {
     lmin = log10(Xmin);
     lmax = log10(Xmax);
   }
   if (meth == SPS_GAMMA) {
     lmin = pow(Xmin, gamma);
     lmax = pow(Xmax, gamma);
   }
   A = (mapmax - mapmin) / (lmax - lmin);
   B = mapmin - ((mapmax - mapmin) * lmin)/(lmax-lmin);
 }
 else {
   A = 1.0;
   B = 0.0;
 }

 switch (type) {
   case SPS_DOUBLE :
     if (mapbytes == 1) {
		/*###CHANGED - ALEXANDRE 11/09/2002*/
       //if (meth == SPS_LOG) {
       //  FASTLOG(double, unsigned char);
       //} else {
         CALCDATA(double, double, unsigned char, log10, pow);
       //}
     } else if (mapbytes == 2) {
       CALCDATA(double, double, unsigned short, log10, pow);
     } else if (mapbytes == 4 || mapbytes == 3) {
       CALCDATA(double, double, unsigned int, log10, pow);
     }
     break;
   case SPS_FLOAT :
     if (mapbytes == 1) {
		/*###CHANGED - ALEXANDRE 11/09/2002*/
       //if (meth == SPS_LOG) {
       //  FASTLOG(float, unsigned char);
       //} else {
         CALCDATA(float, float, unsigned char, log10f, powf);
       //}
     } else if (mapbytes == 2) {
       CALCDATA(float, float, unsigned short, log10f, powf);
     } else if (mapbytes == 4 || mapbytes == 3) {
       CALCDATA(float, float, unsigned int, log10f, powf);
     }
     break;
   case SPS_LONG :
     if (mapbytes == 1) {
       if (meth == SPS_LOG) {
         FASTLOG(int, unsigned char);
       } else {
         CALCDATA(int, float, unsigned char, log10f, powf);
       }
     } else if (mapbytes == 2){
       CALCDATA(int, float, unsigned short, log10f, powf);
     } else if (mapbytes == 4 || mapbytes == 3){
       CALCDATA(int, float, unsigned int, log10f, powf);
     } 
     break;
   case SPS_ULONG :
     if (mapbytes == 1) {
       if (meth == SPS_LOG) {
         FASTLOG(unsigned int, unsigned char);
       } else {
         CALCDATA(unsigned int, float, unsigned char, log10f,powf);
       }
     } else if (mapbytes == 2){
       CALCDATA(unsigned int, float, unsigned short, log10f, powf);
     } else if (mapbytes == 4 || mapbytes == 3){
       CALCDATA(unsigned int, float, unsigned int, log10f, powf);
     } 
     break;
   case SPS_SHORT :
     if (mapbytes == 1) {
       if (cols*rows > 100000) { 
	 CALCPREDATA(short, float, unsigned char, log10f,powf,(-32768),32767);
       } else {
	 CALCDATA(short, float, unsigned char, log10f, powf);
       }
     } else if (mapbytes == 2) {
	CALCDATA_NOMAP(unsigned short, unsigned short);
     } else if (mapbytes == 4 || mapbytes == 3) {
	CALCDATA_NOMAP(unsigned short, unsigned int);
     }
     break;
   case SPS_USHORT :
     if (mapbytes == 1) {
       if (cols*rows > 100000) { 
	 CALCPREDATA(unsigned short,float,unsigned char,log10f, powf,0,65535);
       } else {
	 CALCDATA(unsigned short, float, unsigned char, log10f,powf);
       }
     } else if (mapbytes == 2) {
	CALCDATA_NOMAP(unsigned short, unsigned short);
     } else if (mapbytes == 4 || mapbytes == 3) {
	CALCDATA_NOMAP(unsigned short, unsigned int);
     }
     break;
   case SPS_CHAR :
     if (mapbytes == 1) {
       CALCPREDATA(char, float, unsigned char, log10f, powf, (-128), 127);
     } else if (mapbytes == 2){
       CALCDATA_NOMAP(unsigned char, unsigned short);
     } else if (mapbytes == 4 || mapbytes == 3){
       CALCDATA_NOMAP(unsigned char, unsigned int);
     }
     break;
   case SPS_UCHAR :
     if (mapbytes == 1) {
       CALCPREDATA(unsigned char, float, unsigned char, log10f, powf, 0, 255);
     } else if (mapbytes == 2){
       CALCDATA_NOMAP(unsigned char, unsigned short);
     } else if (mapbytes == 4 || mapbytes == 3){
       CALCDATA_NOMAP(unsigned char, unsigned int);
     }
     break;
 }

 return(Xdata);
}

#define CALCREDUCFAST(datat) \
{\
 int l;\
 register int i, red=reduc;\
 register datat *gtr=data;\
 register datat *ptr= (datat *) ndata;\
\
 for (l=ph;l;l--) {\
   for (i=pw;i;i--,ptr++,gtr+=red) \
     *ptr = *gtr;\
   gtr+=fastjump;\
 }\
}

#define CALCREDUC(datat,calct) \
{\
 calct *line, r2;\
 int linesize, k, l;\
 datat *rline;\
 register int i, j;\
 register datat *gtr=data;\
 register calct *ptr;\
\
 r2 = reduc*reduc;\
 line = (calct*) malloc (linesize=(sizeof(calct) * pw)); \
 for (rline=(datat *)ndata,l=ph;l;l--,rline+=pw) {\
   memset(line, 0,linesize);\
   if (reduc == 2) {\
     for (k=reduc;k;k--) {\
       for (i=pw,ptr=line;i;i--,ptr++) {\
         *ptr += (calct)*gtr++;\
         *ptr += (calct)*gtr++;\
       }\
       gtr+=jump;\
     }\
   }\
   else {\
     for (k=reduc;k;k--) {\
       for (i=pw,ptr=line;i;i--,ptr++) {\
         for (j=reduc;j;j--,gtr++)\
           *ptr += (calct)*gtr;\
       }\
       gtr+=jump;\
     }\
   }\
   {\
     register datat *str;\
     register calct *ltr;\
     register int k;\
     for (str=rline,ltr=line,k=pw;k;k--,str++, ltr++)\
       *str = (datat) (*ltr / r2);\
   }\
 }\
 free(line);\
}

void *SPS_ReduceData (void *data, int type, 
		      int cols, int rows, int reduc, 
		      int *pcols, int *prows, int fastreduction)
{
 int pw, ph, jump, fastjump;
 void *ndata;
 int length = SPS_Size_VLUT(type);
 
 if (reduc == 1) {
   *pcols = cols;
   *prows = rows;
   return(data);
 }

 pw=*pcols = cols / reduc;
 if (pw == 0) {
   pw=*pcols = 1;
 }   
 ph=*prows = rows / reduc;
 if (ph == 0) {
   ph=*prows = 1;
 }   
 jump = cols%reduc;
 fastjump = jump + cols*(reduc-1);

 ndata = (void *)malloc(length * pw * ph);
 if (ndata == (void *)NULL) {
   fprintf(stderr, "Malloc Error in CalcReduction (size = %d), Exit\n",
           length * pw * ph);
   return NULL;
 }

 if (fastreduction) {
   switch (type) {
   case SPS_DOUBLE :
     CALCREDUCFAST(double);
     break;
   case SPS_FLOAT :
     CALCREDUCFAST(float);
     break;
   case SPS_LONG :
     CALCREDUCFAST(int);
     break;
   case SPS_ULONG :
     CALCREDUCFAST(unsigned int);
     break;
   case SPS_SHORT :
     CALCREDUCFAST(short);
     break;
   case SPS_USHORT :
     CALCREDUCFAST(unsigned short);
     break;
   case SPS_CHAR :
     CALCREDUCFAST(char);
     break;
   case SPS_UCHAR :
     CALCREDUCFAST(unsigned char);
     break;
   } 
 } else {
   switch (type) {
   case SPS_DOUBLE :
     CALCREDUC(double,double);
     break;
   case SPS_FLOAT :
     CALCREDUC(float,double);
     break;
   case SPS_LONG :
     CALCREDUC(int,int);
     break;
   case SPS_ULONG :
     CALCREDUC(unsigned int,unsigned int);
     break;
   case SPS_SHORT :
     CALCREDUC(short,int);
     break;
   case SPS_USHORT :
     CALCREDUC(unsigned short,unsigned int);
     break;
   case SPS_CHAR :
     CALCREDUC(char,short);
     break;
   case SPS_UCHAR :
     CALCREDUC(unsigned char,unsigned short);
     break;
   }
 }
 return(ndata);
}

void FillSegment(int pcbyteorder, XServer_Info Xservinfo,
                 unsigned int *val, int from, int to,
                 double R1,double G1,double B1,double R2,double G2,double B2, 
                 int rbit,int gbit,int bbit,int rshift,int gshift,int bshift)
{
 unsigned int *ptr;
 unsigned int R, G, B;
 double Rcol, Gcol, Bcol, Rcst, Gcst, Bcst;
 double coef, width, rwidth, gwidth, bwidth; 
 swaptype value;

/* R = R1 + (R2 - R1) * (i-from) / (to - from)
  palette_col = (int)(R * (2**rbit-1) + 0.5) << rshift |
                (int)(G * (2**gbit-1) + 0.5) << gshift |
                (int)(B * (2**bbit-1) + 0.5) << bshift
*/

 Rcol = (1<<rbit) - 1;
 Rcst = Rcol * R1 + 0.5;
 Gcol = (1<<gbit) - 1;
 Gcst = Gcol * G1 + 0.5;
 Bcol = (1<<bbit) - 1;
 Bcst = Bcol * B1 + 0.5;
 width = (double)(to - from);
 rwidth = Rcol * (R2 - R1) / width;
 gwidth = Gcol * (G2 - G1) / width;
 bwidth = Bcol * (B2 - B1) / width;

 if (pcbyteorder == SPS_LSB) {
   if (Xservinfo.byte_order == SPS_LSB) {
     if (Xservinfo.pixel_size == 3) {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         value.p = (R << rshift) | (G << gshift) | (B << bshift);
         *ptr++ = value.c.b1 << 8 | value.c.b2 << 16 | value.c.b3 << 24;
       }
     } else {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         *ptr++ = (R << rshift) | (G << gshift) | (B << bshift);
       }
     }
   } else {
     if (Xservinfo.pixel_size == 2) {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         value.p = (R << rshift) | (G << gshift) | (B << bshift);
         *ptr++ = value.c.b1 << 8 | value.c.b2;
       }
     } else {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         value.p = (R << rshift) | (G << gshift) | (B << bshift);
         *ptr++ = value.c.b1 << 24 | value.c.b2 << 16 | value.c.b3 << 8;
       }
     }
   }
 } else {
   if (Xservinfo.byte_order == SPS_LSB) {
     if (Xservinfo.pixel_size == 2) {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         value.p = (R << rshift) | (G << gshift) | (B << bshift);
         *ptr++ = value.c.b4 << 8 | value.c.b3;
       }
     } else {
       for (ptr=val+from,coef=0;coef<to-from;coef++) {
         R = Rcst + rwidth * coef;
         G = Gcst + gwidth * coef;
         B = Bcst + bwidth * coef;
         value.p = (R << rshift) | (G << gshift) | (B << bshift);
         *ptr++ = value.c.b4 << 16 | value.c.b3 << 8 | value.c.b2;
       }
     }
   } else {
     for (ptr=val+from,coef=0;coef<to-from;coef++) {
       R = Rcst + rwidth * coef;
       G = Gcst + gwidth * coef;
       B = Bcst + bwidth * coef;
       *ptr++ = (R << rshift) | (G << gshift) | (B << bshift);
     }
   }
 }
}

unsigned int *CalcPalette (XServer_Info Xservinfo, int palette_type) 
{
  static unsigned int *full_palette = NULL;
  static old_type = -1;
  static old_mapbytes = -1;
  unsigned int col;
  int rbit, gbit, bbit, rshift, gshift, bshift, pcbyteorder;
  swaptype val;
  
  if (full_palette && 
      (old_type != palette_type || old_mapbytes != Xservinfo.pixel_size)){
    free(full_palette);
    full_palette = NULL;
  }
  
  if (full_palette == NULL) {
    full_palette = (void*) malloc (0x10000 * sizeof (unsigned int));
    if (full_palette == NULL) {
      fprintf(stderr, "Error - can not malloc memory in FillPalette\n");
      return NULL;
    }
    old_type = palette_type;
    old_mapbytes = Xservinfo.pixel_size;

    val.p = 1;
    if (val.c.b4 == 1) {
      pcbyteorder = SPS_MSB;
    } else {
      pcbyteorder = SPS_LSB;
    }

    col = Xservinfo.red_mask;
    rshift = 0;
    while ((col & 1) == 0) {
      col = col >> 1;
      rshift++;  
    }
    rbit=0;
    while ((col & 1) == 1) {
      col = col >> 1;
      rbit++;
    }

    col = Xservinfo.green_mask;
    gshift = 0;
    while ((col & 1) == 0) {
      col = col >> 1;
      gshift++;  
    }
    gbit=0;
    while ((col & 1) == 1) {
      col = col >> 1;
      gbit++;
    }

    col = Xservinfo.blue_mask;
    bshift = 0;
    while ((col & 1) == 0) {
      col = col >> 1;
      bshift++;  
    }
    bbit=0;
    while ((col & 1) == 1) {
      col = col >> 1;
      bbit++;
    }

    if (palette_type == SPS_GREYSCALE) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x10000, 0, 0, 0, 1, 1, 1, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_TEMP) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x4000, 0, 0, 1, 0, 1, 1, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0x4000, 0x8000, 0, 1, 1, 0, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0x8000, 0xc000, 0, 1, 0, 1, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0xc000, 0x10000, 1, 1, 0, 1, 0, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_RED) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x10000, 0, 0, 0, 1, 0, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_GREEN) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x10000, 0, 0, 0, 0, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_BLUE) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x10000, 0, 0, 0, 0, 0, 1, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_REVERSEGREY) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x10000, 1, 1, 1, 0, 0, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    } else if (palette_type == SPS_MANY) {
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0, 0x2aaa, 0, 0, 1, 0, 1, 1, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0x2aaa, 0x5555, 0, 1, 1, 0, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0x5555, 0x8000, 0, 1, 0, 1, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0x8000, 0xaaaa, 1, 1, 0, 1, 0, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0xaaaa, 0xd555, 1, 0, 0, 1, 1, 0, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
      FillSegment(pcbyteorder, Xservinfo,
                  full_palette, 0xd555, 0x10000, 1, 1, 0, 1, 1, 1, 
                  rbit, gbit, bbit, rshift, gshift, bshift);
    }
  }
  return full_palette;
}

FillPalette (XServer_Info Xservinfo, 
             void *palette, int fmin, int fmax, 
	     int palette_type, int meth, double gamma)
{
  double A, B, round_min;
  double lmin, lmax;
  unsigned int *full_palette;
  int idx;
  
  /*
   SPS_LINEAR:   mapdata = A * data + B
   SPS_LOG   :   mapdata = (A * log(data)) + B
   SPS_GAMMA :   mapdata = A * pow(data, gamma) + B
  */
  if (fmin == 0 && meth != SPS_LINEAR) 
    fmin = 1;
  
  if ((fmax - fmin) != 0) {
    if (meth == SPS_LINEAR) {
      lmin = fmin;
      lmax = fmax;
    }
    if (meth == SPS_LOG) {
      lmin = log10(fmin);
      lmax = log10(fmax);
    }
    if (meth == SPS_GAMMA) {
      lmin = pow(fmin, gamma);
      lmax = pow(fmax, gamma);
    }
    A = 0xffff / (lmax - lmin);
    B = - (0xffff * lmin) / (lmax - lmin);

    if (meth == SPS_LINEAR) {
      round_min = A * fmin + B;
    }
    if (meth == SPS_LOG) {
      round_min = (A * log10(fmin)) + B;
    }
    if (meth == SPS_GAMMA) {
      round_min = (A * pow(fmin,gamma)) + B;
    }
    
    if (round_min < 0.0 && round_min > -1E-5 )
      B += round_min;

  }
  else {
    A = 1.0;
    B = 0.0;
  }

  /* The full palette has always 0x10000 entries of longs; */
  full_palette = CalcPalette (Xservinfo, palette_type);
  
  /* Squeeze the palette into the data range */
  if (Xservinfo.pixel_size == 2) {
    register unsigned short *pal = palette;
    register unsigned short *palend = palette;
    register int j = 0;

    pal += fmin ; palend += fmax;
    if (meth == SPS_LINEAR) {
      j = 0;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint (A * j++)); 
      }
    } else if (meth == SPS_LOG) {
      j = fmin;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint (A * log10 (j++) + B)); 
      }
    } else if (meth == SPS_GAMMA) {
      j = fmin;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint (A * pow(j++, gamma) + B)); 
      }
    }
  } else if (Xservinfo.pixel_size == 4 || Xservinfo.pixel_size == 3) {
    register unsigned int *pal = palette;
    register unsigned int *palend = palette;
    register int j = 0;

    pal += fmin ; palend += fmax;
    if (meth == SPS_LINEAR) {
      j = 0;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint(A * j++)); 
      }
    } else if (meth == SPS_LOG) {
      j = fmin;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint(A * log10 (j++) + B)); 
      }
    } else if (meth == SPS_GAMMA) {
      j = fmin;
      while (pal <= palend) {
	*pal++ = *(full_palette + lrint(A * pow(j++, gamma) + B)); 
      }
    }
  } 
}

void *CreatePalette( int type, int meth, double min, double max, double gamma,
		     int mapmin, int mapmax, 
                     XServer_Info Xservinfo, int palette_type)
{
  int pmin, pmax; /* palette min max values */
  int fmin, fmax; /* fill from these min max values */
  int newsize;
  static void *palette = NULL; 
  static int palette_size = 0;
  int memcorr = 2;
  void *old_palette, *palend;
  int palbytes;

  if (Xservinfo.pixel_size == 1) 
    return NULL;   /* Hardware Palette */

  /* The palette of 3 byte results is 4 byte long */
  palbytes = (Xservinfo.pixel_size == 3) ? 4 : Xservinfo.pixel_size;
  
  if ( type == SPS_FLOAT || type == SPS_DOUBLE || type == SPS_LONG ||
       type == SPS_ULONG ) {
    /* In this case we map first to mapmin and mapmax and use these as an 
       index in the palette */
    fmin = pmin = 0 ; fmax = pmax = mapmax - mapmin;
    meth = SPS_LINEAR; /* We will always map to linear palettes as the mapping
			is not linear - this gives us a higher dynamic range*/
  } else if (type == SPS_USHORT)  {
    /* In all these cases we use the image values directly as an index in the 
       palette */
    pmin = 0 ; pmax = 0xffff;
    fmin = (int) min ; if (fmin < 0) fmin = 0;
    fmax = (int) max ; if (fmax > 0xffff) fmax = 0xffff;
  } else if (type == SPS_UCHAR)  {  
    pmin = 0 ; pmax = 0xff;
    fmin = (int) min ; if (fmin < 0) fmin = 0;
    fmax = (int) max ; if (fmax > 0xff) fmax = 0xff;
  } else if (type == SPS_SHORT )  {
    pmin = 0 ; pmax = 0xffff;
    fmin = (int) min + 0x8000; if (fmin < 0) fmin = 0;
    fmax = (int) max + 0x8000; if (fmax > 0xffff) fmax = 0xffff;
    memcorr = 3;
  } else if (type == SPS_CHAR )  {  
    pmin = 0 ; pmax = 0xff;
    fmin = (int) min + 0x80; if (fmin < 0) fmin = 0;
    fmax = (int) max + 0x80; if (fmax > 0xff) fmax = 0xff;
    memcorr = 3;
  }
    
  /* Size of the alloc is the size of the memory group * 1.5 if we have
     unsigned values. This is done to be able to do the swap with one simple
     memcopy. memcorr is either 2 or 3.
     For 3 mapbytes == 3 the palette is still 4 bytes long;
  */
  newsize = (memcorr * palbytes ) / 2 * (pmax - pmin + 1);
  
  if (palette && newsize > palette_size) {
    free (palette);
    palette = NULL;
  }

  if (palette == NULL) {
    palette = (void *) malloc(newsize);
    if (palette == NULL) {
      fprintf(stderr, "Malloc Error in CreatePalette (size = %d)\n", newsize);
      return NULL;
    }
    palette_size = newsize;
  }
  
  /* Prepare the swap by putting everything 1/2 size higher up */
  if (memcorr == 3) {
    old_palette = palette;
    palette = (void *) ((char *) palette + newsize / 3);
  }
  
  /* Now let's fill the palette */
  FillPalette (Xservinfo, palette, fmin, fmax, palette_type, meth, gamma);
  
  /* Now pad the low and high values */
  if (pmin < fmin) {
    if (Xservinfo.pixel_size == 2) {
      register unsigned short *dest = ((unsigned short *) palette) + pmin; 
      register unsigned short src = *(((unsigned short *) palette) + fmin); 
      register unsigned short *end  = ((unsigned short *) palette) + fmin;
      while (dest < end)
	*dest++ = src;
    } else if (Xservinfo.pixel_size == 4 || Xservinfo.pixel_size == 3) {
      register unsigned int *dest = ((unsigned int *) palette) + pmin; 
      register unsigned int src = *(((unsigned int *) palette) + fmin); 
      register unsigned int *end  = ((unsigned int *) palette) + fmin;
      while (dest < end)
	*dest++ = src;
    }
  }

  if (pmax > fmax) {
    if (Xservinfo.pixel_size == 2) {
      register unsigned short *dest = ((unsigned short *) palette) + fmax +1; 
      register unsigned short src = *(((unsigned short *) palette) + fmax); 
      register unsigned short *end  = ((unsigned short *) palette) + pmax;
      while (dest <= end)
	*dest++ = src;
    } else if (Xservinfo.pixel_size == 4 || Xservinfo.pixel_size == 3) {
      register unsigned int *dest = ((unsigned int *) palette) + fmax + 1; 
      register unsigned int src = *(((unsigned int *) palette) + fmax); 
      register unsigned int *end  = ((unsigned int *) palette) + pmax;
      while (dest <= end)
	*dest++ = src;
    }
  }

  /* Now the palette has to be swaped over when we have signed image values */
  if (memcorr == 3) {
    palette = old_palette;
    palend = (void *) ((char *) palette + newsize / 3 * 2);
    memcpy (palette, palend, newsize / 3);
  }

  return palette;
}

double SPS_GetZdata(void *data, int type, int cols, int rows, int x, int y)
{
 int ind;

 ind = y*cols + x;
 if (ind >= (cols*rows))
   ind = cols*rows-1;
 switch (type) {
   case SPS_DOUBLE :
     return(*((double *)data + ind));
     break;
   case SPS_FLOAT :
     return((double)(*((float *)data + ind)));
     break;
   case SPS_LONG :
     return((double)(*((int *)data + ind)));
     break;
   case SPS_ULONG :
     return((double)(*((unsigned int *)data + ind)));
     break;
   case SPS_SHORT :
     return((double)(*((short *)data + ind)));
     break;
   case SPS_USHORT :
     return((double)(*((unsigned short *)data + ind)));
     break;
   case SPS_CHAR :
     return((double)(*((char *)data + ind)));
     break;
   case SPS_UCHAR :
     return((double)(*((unsigned char *)data + ind)));
     break;
 }
}

void SPS_PutZdata(void *data, int type, int cols, int rows, int x, int y, 
		  double z)
{
  int ind;
  
  ind = y*cols + x;
  if (ind >= (cols*rows))
    ind = cols*rows-1;
  switch (type) {
  case SPS_DOUBLE :
    *((double *)data + ind) = z;
    break;
  case SPS_FLOAT :
    *((float *)data + ind) = z;
    break;
  case SPS_LONG :
    *((int *)data + ind) = z;
    break;
  case SPS_ULONG :
    *((unsigned int *)data + ind) = z;
    break;
  case SPS_SHORT :
    *((short *)data + ind) = z;
    break;
  case SPS_USHORT :
    *((unsigned short *)data + ind) = z;
    break;
  case SPS_CHAR :
    *((char *)data + ind) = z;
    break;
  case SPS_UCHAR :
    *((unsigned char *)data + ind) = z;
    break;
  }
  return;
}

#define CALCSTAT(ty, calct)\
{\
 register calct integr=0;\
 register ty *ptr;\
\
 for (ptr=(ty *)data,i=n;i;ptr++,i--)\
   integr += (calct)(*ptr);\
 aver = (double)integr / (double)n;\
 for (ptr=(ty *)data,i=n;i;ptr++,i--) {\
   val = (double)(*ptr) - aver;\
   std += (val*val);\
 }\
 integ = (double)integr;\
}

void SPS_CalcStat(void *data, int type, int cols, int rows,
              double *integral, double *average, double *stddev)
{
 double integ;
 int n;
 register double std=0.0, val, aver;
 register int i;
 int length = SPS_Size_VLUT(type);

 n = cols*rows;

 switch (type) {
   case SPS_DOUBLE :
     CALCSTAT(double,double);
     break;
   case SPS_FLOAT :
     CALCSTAT(float,double);
     break;
   case SPS_LONG :
     CALCSTAT(int,double);
     break;
   case SPS_ULONG :
     CALCSTAT(unsigned int,double);
     break;
   case SPS_SHORT :
     CALCSTAT(short,double);
     break;
   case SPS_USHORT :
     CALCSTAT(unsigned short,double);
     break;
   case SPS_CHAR :
     CALCSTAT(char,int);
     break;
   case SPS_UCHAR :
     CALCSTAT(unsigned char, unsigned int);
     break;
 }

 std = std / (double)(n-1);
 std = sqrt(std);
 *integral = integ;
 *average = aver;
 *stddev = std;
}

#define DATADIST(ty)\
{\
 register ty *ptr = (ty *)data;\
 register double *dtr=*ydata;\
 register int i, ind;\
 for (i=n;i;i--, ptr++) {\
   ind = (int)(((double)*ptr-min) / step);\
   dtr[ind]++;\
 }\
}

void SPS_GetDataDist(void *data, int type, int cols, int rows,
			 double min, double max, 
			 int nbar, double **xdata, double **ydata)
{
 double step, val, start, *ptr;
 int n=cols*rows;


 step = (max - min) / (double)nbar;

 if (step == 0.0) {
   *xdata = (double *)malloc(sizeof(double));
   if (*xdata == (double *)NULL) {
     fprintf(stderr, "Malloc Error in GetDataDistribution 1 (size=%d), Exit\n",
                     sizeof(double));
     exit(2);
   }

   *ydata = (double *)malloc(sizeof(double)*2); /* on a honte */
   if (*ydata == (double *)NULL) {
     fprintf(stderr, "Malloc Error in GetDataDistribution 2 (size=%d), Exit\n",
                     sizeof(double)*2);
     exit(2);
   }

   (*ydata)[0] = (*ydata)[1] = (double)(cols*rows);
   (*xdata)[0] = (double)max;

   return;
 }

 *xdata = (double *)malloc(sizeof(double)*nbar);
 if (*xdata == (double *)NULL) {
   fprintf(stderr, "Malloc Error in GetDataDistribution 3 (size=%d), Exit\n",
                   sizeof(double)*nbar);
   exit(2);
 }

 *ydata = (double *)malloc(sizeof(double)*(nbar+1)); /* on a honte */
 if (*ydata == (double *)NULL) {
   fprintf(stderr, "Malloc Error in GetDataDistribution 4 (size=%d), Exit\n",
                   sizeof(double)*(nbar+1));
   exit(2);
 }

 start = min + 0.5 * step;

 memset(*ydata, 0,(nbar+1)*sizeof(double)); /*** NOT SURE ***/
 for (ptr=*xdata,val=start;val<max;val+=step,ptr++)
   *ptr = val;

 switch (type) {
   case SPS_DOUBLE :
     DATADIST(double);
     break;
   case SPS_FLOAT :
     DATADIST(float);
     break;
   case SPS_LONG :
     DATADIST(int);
     break;
   case SPS_ULONG :
     DATADIST(unsigned int);
     break;
   case SPS_SHORT :
     DATADIST(short);
     break;
   case SPS_USHORT :
     DATADIST(unsigned short);
     break;
   case SPS_CHAR :
     DATADIST(char);
     break;
   case SPS_UCHAR :
     DATADIST(unsigned char);
     break;
 }

 (*ydata)[nbar-1] += (*ydata)[nbar]; /* on a honte  pour *ptr = max */\
}












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
/* FIXTHIS - double the type defines to make sps_lut independent */
#ifndef SPS_DOUBLE 
#define SPS_DOUBLE      0
#define SPS_FLOAT       1
#define SPS_LONG        2
#define SPS_ULONG       3
#define SPS_SHORT       4
#define SPS_USHORT      5
#define SPS_CHAR        6
#define SPS_UCHAR       7
#define SPS_STRING      8
#endif

/* contian description of how the colors are represented for the Xserver */
typedef struct XServer_Info {
  int byte_order;
  int pixel_size;
  unsigned int red_mask;
  unsigned int green_mask;
  unsigned int blue_mask;
} XServer_Info;

/* GILLES DOC is as usual absent ... but at least these are the different
   possibilities in the code:

   PCByteorder ServerByteorder  pixel_size  remapping
   ----------------------------------------------------
   LSB         LSB              3           b3 b2 b1 00
   LSB         LSB              not 3       none
   LSB         MSB              2           00 00 b1 b2
   LSB         MSB              not 2       b1 b2 b3 00
   LSB         LSB              2           00 00 b4 b3
   LSB         LSB              not 2       00 b4 b3 b2
   LSB         MSB              always      none

with the union {b1 b2 b3 b4} is a long and MSB if from long = 1 follows that
b1 = 1  (in our case LINUX is LSB and Solaris and HP are MSB)

Normally you will have to get this information from the XServer
*/

/*
  High level function to transform data of a certain type (SPS_FLOAT,
  SPS_DOUBLE, SPS_LONG, SPS_ULONG, SPS_SHORT, SPS_USHORT, SPS_CHAR, SPS_UCHAR)
  into a representation suitable for display on a screen. 

  The size of the input array is given in cols and rows. 

  The data array can be reduced by a factor reduc. The number of output
  pixel in each direction will be divided by this factor. The values
  prows and pcols tell the size of this new array in pixel. The flag fastreduc
  can be set to do a faster reduction by not doing an average of the
  pixel to be reduced but skipping the pixels.

  The meth can be set to SPS_LINEAR, SPS_LOG or SPS_GAMMA to define in with
  way the data is mapped to the colors. The following formulas are used:

   SPS_LINEAR:   mapdata = A * data + B
   SPS_LOG   :   mapdata = (A * log(data)) + B
   SPS_GAMMA :   mapdata = A * pow(data, gamma) + B

  The autoscale flag has to be set to know if the minimum and maximum
  data values should be found automatically or if they are given in 
  min, max. If autoscale is used the resulting minimum and maximum values
  are returned in min and max. If the method is not SPS_LINEAR, the routine
  will not return the smallest data value but the smallest positive data
  values. 
    ! Check for SPS_SHORT SPS_USHORT SPS_CHAR SPS_UCHAR behaviour !

  The values of mapmin and mapmax matter mainly when a hardware palette 
  is used. Both values have to be in [0,255]. You can decide the minimum 
  and maximum output value in this case. If you specify for example mapmin=100
  and mapmax=250 you can leave colors for the window manager. If you would
  like to use your private colormap you could use mapmin=0 mapmax=250. In this
  example we leave 5 colors for other elements in our application starting
  from index 251 to 255. 

  The mapbytes value tell the routine the depth of the display you will use.
  The distribution of colors in this bytes is assumed to be fixed. You
  can have the following cases:
      0 : Hardware palette. The colors are given by the hardware. Only a
          value between mapmin and mapmax is produced.
      2 : The RGB parts of the 16 bits are supposed to be 5 bits per color.
          The remaining bit is set to 0.
      3 : The RGB parts of the 24 bits are supposed to be 8 bits per color.
      4 : The RGB parts of the 32 bits are supposed to be 8 bits per color.
          The remaining 8 bits are set to 0.

   The additional parameter maptype is given to distinguish different 
   distributions of the colorbits:
      SPS_NORMAL: Standard distribution (i.e. for 16 bits (5/5/5) for (R/G/B)
      SPS_16BIT : Special case for 16 bit depth (normally 15 bits)
                  distribution is (5/5/6) for (R/G/B)
   The additional parameter maporder represent the byte order of the X server
   used. it can be either :
      SPS_LSB : low byte first
      SPS_MSB : high byte first

   The palette code describes the colors of the palette and is only used if
   we do not use a hardware palette.
      SPS_GREYSCALE : Gradient from white to black. There are 32 different
            colors possible on a 2 byte display and 256 on an 4 byte display.
      SPS_TEMP : A color range from red to blue. There are 128 different 
            colors possible on a 2 byte display and 1024 on a 4 byte display.

   The palette is returned in palette. The palette is either 2 or 4 bytes
   deep and has pal_entries number of entries. This can be used to present
   a color bar to the user. To be able to asign data numbers to every entry
   in the color table you have to understand the different palette modes.
   They are to some extend explained in the documentation for SPS_MapData 
   The palette data must not be freed!!!!
   
   The return values points to an array of 1, 2, or 4 byte pixels of the
   size prows pcols. The data should be freed when not longer in use.
   A return value of NULL indicates an error.
 */

void *SPS_PaletteArray (void *data, int type, int cols, int rows, 
          int reduc, int fastreduc, int meth, double gamma, int autoscale, 
	  int mapmin, int mapmax, 
          XServer_Info Xservinfo, int palette_code,
	  double *min, double *max, int *pcols, int *prows,
	  void **palette, int *pal_entries);

/* 
   Function to create a linear palette with 3 bytes used in mapbyte == 0
   mode where we have a hardware palette which is managed outside the
   library. DO NOT FREE THE RETURN POINTER - it is used again in the library.
   The returned colormap contains 3 bytes per entry and max - min + 1 entries.
   The colormap is put between the values mapmin and mapmax. Outside this 
   region the colors are set to the border values.
*/

unsigned char *SPS_SimplePalette ( int min, int max, 
                                   XServer_Info Xservinfo, int palette_code); 


/*
  Get the data element from an array <data> of type and size cols * rows
  the position of the element is x cols and y rows
*/

double SPS_GetZdata(void *data, int type, int cols, int rows, int x, int y);

/*
  Puts the data element from an array <data> of type and size cols * rows
  the position of the element is x cols and y rows to z
*/
void SPS_PutZdata(void *data, int type, int cols, int rows, int x, int y, 
		  double z);

/* 
   Calculate some statics on the array data of type.
   The size of the array is given with cols and rows.
   Calculated are: integral, average and std deviation 
*/
void SPS_CalcStat(void *data, int type, int cols, int rows,
              double *integral, double *average, double *stddev);



/* 
   Calculate a histogram of the image. array in data type type with size
   rows * columns. From minimum to maximum in nbar steps. The result is
   put into two arrays xdata and ydata. With xdata being the xvalue and
   yvalue beiing the y-value of every histogram point.
*/
void SPS_GetDataDist(void *data, int type, int cols, int rows,
			 double min, double max, 
			 int nbar, double **xdata, double **ydata);

/* --------------------------- Lower Level ------------------------------ */
/* The following functions are called from the above functions            */

/* 
   Maps the data array of type type and size rows * cols either with
   a linear, logarithmic or gamma corrected scale. meth in this case SPS_LOG
   SPS_LINEAR or SPS_GAMMA. mapbytes and the data type will influence the 
   outputdata in the following way:
   
   SPS_LINEAR:   mapdata = A * data + B
   SPS_LOG   :   mapdata = (A * log(data)) + B
   SPS_GAMMA :   mapdata = A * pow(data, gamma) + B

   mapbytes indicates the depth of the palette or is zero if a hardware
   palette is used. It can take the following values:

       0 : Xmin will be maped to mapmin and Xmax to mapmax.
           mapmin and mapmax have to bin in [0,255]. This mode is thought to
	   be used with hardware plalettes and the parameter pal is ignored.
       2 : The color table is 2 bytes deep. In this case, for every color 
           there are 5 bits and the last bit is ignored. The mapping is 
	   different depending on the type of the input data.
	   
           SPS_LONG, SPS_ULONG, SPS_DOUBLE, SPS_FLOAT : The mapping is done
	   according to the above formula. The color palette is supposed
	   to contain the coresponding RGB values for each of the values 
	   between mapmin and mapmax.

	   SPS_USHORT SPS_SHORT SPS_CHAR SPS_UCHAR : No mapping is done !!!
           The actual data values are used direcly as index in the color
	   table.

   3 or 4: The color table is 4 byte deep. There are 3 times 8 bit for
           RGB and 8 bit are ignored. The mapping is described above.
   
   Maybe a less confusing way how these parameters are interpreted is to
   think of the mapping and indexing process as sequential processes.
   In general the source data is first mapped to a range of values which
   are used as indexes in a table of colors.
   When there is a hardware palette (mode mapbytes = 0) then this routine
   has to produce values between 0 and 255. There values can be further
   restricted with mapmin and mapmax to allow other colors to stay in same
   color table. This can be used to give the window manager enough colors
   to operate without creating a private colormap. A standard setting is
   mapmin = 100 mapmax = 250.

   When there is no hardware palette, a palette has to be provided to the
   routine. Every entry of this palette is either 2bytes or 4bytes deep
   depending on the screen color depth. 
   Because of performance reasons we do not want to map data first to 
   other values and use these calculated values as an index into the
   palette if we can use the data directly as an index in the color 
   table. Therefore there is a distinction between types which contain
   data with many different values (SPS_LONG SPS_ULONG SPS_FLOAT SPS_DOUBLE)
   and data with a restricted number of different values (SPS_SHORT SPS_USHORT
   SPS_CHAR SPS_UCHAR). This distinction complicates of course the calculation
   of the palette which has to be provided. 
   
   The palette in the case of input types with restricted values must take 
   into account scaling and the different ways of mapping (linear, log, gamma)
   as there is no other mapping taken place. The routine SPS_PaletteArray does
   take all this into acount and should be used in almost all cases.

*/ 

unsigned char *SPS_MapData(void *data, int type, int meth, int cols, int rows,
			   double Xmin, double Xmax, double gamma,
			   int mapmin, int mapmax, int mapbytes, void *pal);


/* 
   Produce a new array reduced by a factor of reduc. The reduction can be
   done either fast (by setting fastreduction to 1) which will skip reduc - 1
   pixel or accurate by taking the average of all the pixels. The old array
   has row rows and col columns. The new array will have prows and pcols 
*/
void *SPS_ReduceData (void *data, int type, int cols, int rows, int reduc, 
		      int *pcols, int *prows, int fastreduction);
     
/*
  Search through an array startig at <data>  of type <type>. The size
   of the array being rows * cols. The flag tells the function if 
   it should calculate the min/max (1), the positive minimum (2) or both (3).
   Results are returned in min max minplus
*/
void SPS_FindMinMax(void *data, int type, int cols, int rows,
		    double *min, double *max, double *minplus, int flag);


#define SPS_LINEAR 0
#define SPS_LOG 1
#define SPS_GAMMA 2


#define SPS_GREYSCALE 1
#define SPS_TEMP 2
#define SPS_RED 3
#define SPS_GREEN 4
#define SPS_BLUE 5
#define SPS_REVERSEGREY 6
#define SPS_MANY 7 

#define SPS_NORMAL 0
#define SPS_16BIT 1

#define SPS_LSB 0
#define SPS_MSB 1




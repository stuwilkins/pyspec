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
/* Small documentation for the library:
   Your program might get a list of available arrays with the two functions
   ------------------------------------------------------------------------ 
     SPS_GetNextSpec (flag) 
     SPS_GetNextArray (specversion, flag)
     
     Just call the function SPS_GetNextSpec until a NULL pointer is returned.
     The first time you call the function you have to set the flag to 0 to
     indicate it should start from the beginning. For every spec version call
     the function SPS_GetNextArray to get all the shared memory arrays from
     this spec version. 

     Example: 
         for (i=0; spec_version = SPS_GetNextSpec (i) ; i++)
            for (j=0; array = SPS_GetNextArray (spec_version, j) ; j++)
              printf("%s %s\n",spec_version, array);


   Your program can know in which state a certain spec version is in
   -----------------------------------------------------------------
     SPS_GetSpecState (version)
     
     The exact meaning of this flag will be documented

   Your program can get the data from the shared memory array with
   ---------------------------------------------------------------
      SPS_GetDataCopy (spec version, array, data format, &rows, &cols)
      
   You can get a copy of the array data in the format you want. Possible 
   choices for the format are: SPS_DOUBLE SPS_FLOAT SPS_LONG SPS_ULONG
   SPS_SHORT SPS_USHORT SPS_CHAR SPS_UCHAR where DOUBLE and FLOAT are 
   foating point numbers and the other integer values. The U is a short
   cut for unsigned. LONG is stored in 4 byte, SHORT in 2 byte and CHAR in
   on byte. 
   Only one private buffer per shared memory array is possible with this
   function. This buffer is reused at the next call to SPS_GetDataCopy.
   The buffer can be given back with the function

      SPS_FreeDataCopy (spec version, array name)

   but as there is only one per spec version this is rarely necessary. 

      SPS_GetDataRow (version, array, data format, row, columns, &actual cols)
      SPS_GetDataCol (version, array, data format, col, rows, &actual rows)

   Will only give you one row or one column. This functions will use the
   same buffer as SPS_GetDataCopy. The buffer can be freed with 
   SPS_FreeDataCopy but it is reused anyhow. 
   
   To know if the values in your buffer still correspond to the data in
   the shared memory array, there is the function

      SPS_IsUpdated (spec version, array name)

   This function will return the integer value 1 if it is necessary to 
   call SPS_GetDataCopy again. The function SPS_IsUpdated will tell you
   if the shared memory array has been updated since the last time you
   called this function. A return value of -1 indicated that this shared 
   memory array is currently not accessible. To get this information,
   the function has to call the external program ipcs. It is therefore 
   rather slow. Please lower the rate at which you call SPS_IsUpdated
   in this case.


   Example:
      for (;;) {
        if (SPS_IsUpdated ("spec","a_test") == 1) {
          printf("Changed:\n");
          dbuf = SPS_GetDataCopy ("spec","a_test", SPS_DOUBLE, &rows, &cols);
          for (k_row = 0; k_row < rows; k_row ++) {
            for (k_col = 0; k_col < cols < 7; k_col ++)
              printf ("%10.10g ",*(double_buf + k_row * cols + k_col));
            printf("\n");
          }
        }
        sleep(1);
      }

      SPS_IsUpdated (spec version, array name)

   This function is an alternative form of the above. It returns the 
   current update counter for this array or -1 if not longer updated. It
   is used in some cases where the simpler function SPS_IsUpdated is not 
   sufficient (For example if you have multiple subroutines in your program 
   which access the shared memory independently. The SPS routines can not 
   distinguish between calls from the same subroutine twice or single calls 
   from different subroutines)
   You should keep the old counter value and compare the returned value to
   the old. If it is changed you have to update your array.


   Functions to read shared memory string arrays for resources
   -----------------------------------------------------------

   The shared memory string array must be organized in the format:
      identifier=value 
   in a single row of a two dimensional array for these functions to work.

      SPS_GetEnvStr (spec_version, array, identifier)

   returns the value of this identifier or NULL if not found. The function

      SPS_PutEnvStr (spec_version, array, identifier, value)

   puts a value in a row of a 2 dimensional shared memory string array. If
   the identifier was already present before, its value will be updated. 

   Example:
       46.SPEC> shared string array text[10][100] #10 strings of max length 100
       47.SPEC> p text
          xaxis=hello you
          yaxis=Counts
          new=23.45

       SPS_GetEnvStr("spec", "text", "xaxis") returns "hello you"

   More functions to read data from the shared memory
   -------------------------------------------------  

      SPS_GetArrayInfo (spec version, array, &rows, &cols, &type, &flag)

   informs you about the properties of the shared memory array.

      SPS_GetDataPointer (spec version, array name, write_flag)
 
   returns a pointer to the data in the shared memory array. Do not keep
   this pointer around. The memory might have been deleted and recreated 
   and this pointer will not point anymore to the right data. We recommand
   to call the SPS_GetPointer function before each operation on the
   shared memory. The write_flag tells the function if you would like
   to update the shared memory contents or just read it.

   If you want to make a copy of the data in the shared memory array to
   your own buffer or the copy your buffer to the shared memory you can use 
   the functions: 

      SPS_CopyFromShared (spec version, array name, 
                     &buffer, data format, items_in_buffer);
      SPS_CopyToShared (spec version, array name, 
                     &buffer, data format, items_in_buffer);
      SPS_CopyRowFromShared (spec version, array name, &buffer, data_format, 
               row to copy, no of columns to copy, &actual columns copied)
      SPS_CopyColFromShared (spec version, array name, &buffer, data_format, 
               column to copy, no of rows to copy, &actual rows copied)
      SPS_CopyRowToShared (spec version, array name, &buffer, data_format, 
               row to copy, no of columns to copy, &actual columns copied)
      SPS_CopyColToShared (spec version, array name, &buffer, data_format, 
               column to copy, no of rows to copy, &actual rows copied)

   These functions are a little bit different from
   the SPS_GetDataCopy function in that you will have to manage your buffer
   yourself. (Allocate the memory, free the memory, ...)

      SPS_AttachToArray (spec version, array_name, write_flag)
     
   can be called anytime. It will try to attach to the shared memory. If it
   fails to attach to the memory it will still remember the fact that you want
   to stay attached to the shared memory and will attach you the next time 
   you call one of the SPS_ functions for this spec version and array. The
   write_flag has to be set to 1 if you want to write to this array in
   the future. 

      SPS_DetachFromArray (spec version, array name)

   Detaches from the array or informs the library that we don't want to 
   stay attached to the array. 
   If you write to the shared memory you have to tell all the processes
   which are reading he shared memory that the contents has changed. To do
   that you can call
 
     SPS_UpdateDone (spec version, array name)


   If you need to know the number of bytes for a certain type call:

     SPS_Size (type)

   Creating your own shared memory arrays with
   -------------------------------------------

      SPS_CreateArray (spec version, array name,rows,cols, data format, flags)
   
   This function creates shared memory arrays like SPEC would do that. For
   every spec version there is a shared memory array which holds information
   about this spec version. This function will create this information
   array together with the actual data array. If you want to create 
   multiple arrays it is much more efficient to use always the same 
   name for the spec version. You can not create a shared memory array for
   a spec version which already exists and has not been created by your 
   process. 

     SPS_CleanUpAll()
   
   is used to delete all the shared memory you created and free all the 
   memory the libray allocated. Do not use any buffer pointer you got
   from the library from then onwards. You can continue to run and use
   the SPS_ functions again, but in most cases this incstruction will be 
   called just before quiting the program.

*/

#define SPS_IS_ARRAY    0x0002
#define SPS_IS_MCA      (0x0004|SPS_IS_ARRAY)
#define SPS_IS_IMAGE    (0x0008|SPS_IS_ARRAY)

/* structure flags */
#define SPS_TAG_STATUS   0x0001
#define SPS_TAG_ARRAY    0x0002
#define SPS_TAG_MASK     0x000F  /* User can't change these bits */
#define SPS_TAG_MCA      0x0010
#define SPS_TAG_IMAGE    0x0020
#define SPS_TAG_SCAN     0x0040
#define SPS_TAG_INFO     0x0080

/* array data types */
#define SPS_DOUBLE      0
#define SPS_FLOAT       1
#define SPS_LONG        2
#define SPS_ULONG       3
#define SPS_SHORT       4
#define SPS_USHORT      5
#define SPS_CHAR        6
#define SPS_UCHAR       7
#define SPS_STRING      8



/*
  Input: Type code
  Returns: Size of this type i.e. SPS_LONG returns 4
*/
int SPS_Size(int type);

/*
  char * SPS_GetNextSpec (int flag)
  Input: Flag to know if this is the first call to SPS_GetNextSpec
         0: Get first in list, 1: Get next 
  Returns: Name of the SPEC version or NULL if no more in the list
*/

char * SPS_GetNextSpec (int flag) ;

/*
  char * SPS_GetNextArray (char * version, int flag)
  Input: version : name of SPEC version. 
         Flag to know if this is the first call to SPS_GetNextArray
         0: Get first in list, 1: Get next 
  Returns: Name of the SPEC array name or NULL if no more in the list
*/

char * SPS_GetNextArray (char * fullname, int flag);
/* 
   Read the state the particular SPEC version is in for the moment
   Input: version : specversion with PID if necessary (spec(1234) or fourc)
   Returns: State
*/

long int SPS_GetSpecState (char *version);

/* 
   Attaches to a SPEC array. 
   Input: fullname specversion with PID if necessary (spec(1234) or fourc)
          array_name: The name of the SPEC array (i.e. MCA_DATA)
	  write_flag: One if you intend to modify the shared memory
   Output: 0 if OK
           1 if error
*/

int SPS_AttachToArray (char *fullname, char *array, int write_flag);


int SPS_DetachFromArray (char *fullname, char* array);

/* 
   Gives you the pointer to the data area of the SPEC array. If the process
   is not currently attached it will be attached after the call.
   Input: fullname : Spec version
          array : Name of the array
	  write_flag : Tells the routine if we would like to write to the array
   Returns: NULL error
            void * to the data area. Do not remember this pointer as most 
	    of the SPS_ functions can change itr (In case the other 
	    party quit and recreated the shared memory) 
	    
*/

void * SPS_GetDataPointer (char *fullname, char* array, int write_flag) ;

/* 
   You should return the pointer to the library if you do not use it anymore.
   If you returned the pointer as many times as you got it you will be 
   detached from the shared memory and the pointer is not valid anymore.
   Input: fullname : Spec version
          array : Name of the array
   Returns: 0 success
            1 error
*/

int SPS_ReturnDataPointer (void *pointer);

/*
  Used to read a shared string array in where every line is in the
  format identifier=value
    Input: spec_version : The spec version
           array_name : the shared memory string array which holds the
                        lines id=value
	   identifier : the identifier
    Result: NULL if we could not connect to the array or if the identifier
                 did not exist
            value : Do not free this pointer. Do not just keep this pointer
                    around but make a copy of the contents as the storage
                    space is reused at every call to this function.
*/

char * 
SPS_GetEnvStr (char *spec_version, char *array_name, char *identifier);


/*
  Used to write a shared string array in where every line is in the
  format identifier=value
    Input: spec_version : The spec version
           array_name : the shared memory string array which holds the
                        lines id=value
	   identifier : the identifier
	   set_value : The value you would like to put in the identifier 
    Result: 1 error
            0 success

*/

int
SPS_PutEnvStr (char *spec_version, char *array_name, 
	       char *identifier, char *set_value);

/*
  Get all the keys in an environment array

    Input: spec_version : The spec version
           array_name : the shared memory string array which holds the
                        lines id=value
	   flag       : 0  start from first key
                        >0 get next key 		
    Result: NULL if no more keys in environment
*/

char *
SPS_GetNextEnvKey (char *spec_version, char *array_name, int flag); 

/* 
   Copies the data from the shared memory SPEC array to the user's buffer.
   The type of the data in the shared array and the type the user wants
   in his buffer is taken into account. The routine should work efficient
   even if both types are equal. The user has to provide the number of 
   items in the buffer. This number is used to check for buffer overflow.
   If a possible buffer overflow is detected only the items_in_buffer
   are copied. 
   Input: fullname : name of the specversion
          array: name of the array in SPEC
          buffer: pointer to our buffer
	  my_type: In which format do we want the results
	  items_in_buffer: buffersize in number of items of type my_type 
	                   in buffer
   Returns: 0 success
            1 overflow , but copy done
           -1 error Nothing done
*/

int SPS_CopyFromShared (char *fullname, char *array, void *buffer, int my_type,
		      int items_in_buffer);

/* 
   Copies the data to the shared memory SPEC array from the user's buffer.
   The type of the data in the shared array and the type the user wants
   in his buffer is taken into account. The routine should work efficient
   even if both types are equal. The user has to provide the number of 
   items in the buffer. This number is used to check for buffer overflow.
   If a possible buffer overflow is detected only the items_in_buffer
   are copied. 
   Input: fullname : name of the specversion
          array: name of the array in SPEC
          buffer: pointer to our buffer
	  my_type: In which format do we want the results
	  items_in_buffer: buffersize in number of items of type my_type 
	                   in buffer
   Returns: 0 success
            1 overflow , but copy done
           -1 error Nothing done
*/

int SPS_CopyToShared (char *fullname, char *array, void *buffer, int my_type,
		      int items_in_buffer);

/*
  Copy the data in the shared memory array to a private buffer. Only one
  private buffer per shared memory array is allowed. You can call 
  SPS_FreeDataCopy to free the memory by the buffer but it is not necessary.
  If the buffer is not freed it is reused the next time you call the
  function SPS_GetDataCopy for fullname-array. 
  Input:
         fullname : The name of the specversion
	 array : The name of the array in this SPEC version
	 my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 rows : pointer to integer. Will be filled with number of rows 
	        if not NULL
	 cols : pointer to integer. Will be filled with number of cols 
	        if not NULL
  Returns: Pointer to the data array (DO NOT USE free() to free this buffer)
*/

void * SPS_GetDataCopy (char *fullname, char *array, int my_type, 
			int *rows_ptr, int *cols_ptr);

int SPS_FreeDataCopy (char *fullname, char *array);

/*
  Copy a row of data from the shared memory array to a buffer. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 buf : A pointer to the buffer.
         my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 row : integer. Which row do you want.
	 col : integer. How many columns you want to copy. If this is 0 then 
	       all the columns are copied. 
	 act_cols : pointer to int. If not NULL then this is filled with the
               number of columns actually copied.
  Returns: 0 success
           1 error
*/

int SPS_CopyRowFromShared (char *name, char *array, void *buf, 
     int my_type, int row, int col, int *act_cols);

/*
  Copy a column of data from the shared memory array to a buffer. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 buf : A pointer to the buffer.
         my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 col : integer. Which column do you want.
	 row : integer. How many rows you want to copy. If this is 0 then 
	       all the columns are copied. 
	 act_rows : pointer to int. If not NULL then this is filled with the
               number of rows actually copied.
  Returns: 0 success
           1 error
*/

int SPS_CopyColFromShared (char *name, char *array, void *buf, 
     int my_type, int col, int row, int *act_rows);

/*
  Copy a row of data to the shared memory array from a buffer. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 buf : A pointer to the buffer.
         my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 row : integer. Which row do you want.
	 col : integer. How many columns you want to copy. If this is 0 then 
	       all the columns are copied. 
	 act_cols : pointer to int. If not NULL then this is filled with the
               number of columns actually copied.
  Returns: 0 success
           1 error
*/

int SPS_CopyRowToShared (char *name, char *array, void *buf, 
     int my_type, int row, int col, int *act_cols);

/*
  Copy a column of data to the shared memory array from a buffer. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 buf : A pointer to the buffer.
         my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 col : integer. Which column do you want.
	 row : integer. How many rows you want to copy. If this is 0 then 
	       all the columns are copied.  
	 act_rows : pointer to int. If not NULL then this is filled with the
               number of rows actually copied.
  Returns: 0 success
           1 error
*/

int SPS_CopyColToShared (char *name, char *array, void *buf, 
     int my_type, int col, int row, int *act_rows);

/*
  Copy a row of data from the shared memory array to a private buffer. Only 
  one private buffer per shared memory array is allowed. You can call 
  SPS_FreeDataCopy to free the memory by the buffer but it is not necessary.
  If the buffer is not freed it is reused the next time you call the
  function SPS_GetDataRow for fullname-array. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 row : integer. Which row do you want.
	 col : integer. How many columns to you want. If this is 0 then all
	       the columns are copied. 
	 act_cols : pointer to int. If not NULL then this is filled with the
               number of columns actually copied.
  Returns: Pointer to the data array (DO NOT USE free() to free this buffer)
*/

void * SPS_GetDataRow (char *name, char *array, int my_type, 
			int row, int col, int *act_cols);

/*
  Copy a column of data from the shared memory array to a private buffer. 
  Only one private buffer per shared memory array is allowed. You can call 
  SPS_FreeDataCopy to free the memory by the buffer but it is not necessary.
  If the buffer is not freed it is reused the next time you call the
  function SPS_GetDataCol, SPS_GetDataRow or SPS_GetData for fullname-array. 
  Input:
         name : The name of the specversion
	 array : The name of the array in this SPEC version
	 my_type : The data type how you would like to have the buffer. This
	           data type does not have to be the same as the data type
		   of the shared memory array. In this case it will be 
		   transformed in the data type you asked for (fast)
	 col : integer. Which column do you want.
	 row : integer. How many rows you want. If this is 0 then all
	       the rows are copied. 
	 act_rows : pointer to int. If not NULL then this is filled with the
               number of rows actually copied.
  Returns: Pointer to the data array (DO NOT USE free() to free this buffer)
*/

void * SPS_GetDataCol (char *name, char *array, int my_type, 
			int col, int row, int *act_rows);

/* 
   Tells you if the data in the shared memory array changed since you last 
   called this function with this opac pointer or since SPS_ConnectToArray 
   if this is your first call to this function.
   If you are not currently attached to the memory the function will
   attach, read and detach after.

   Input: fullname and array
   Returns: 0 No change
           -1 Error - memory not longer updated
            1 Contents changed - new valid data
*/
int SPS_IsUpdated (char *fullname, char *array); 
  
/* 
   Returns the current update counter for this array or -1 
   if not longer updated. This function is used in some cases where the
   simpler function SPS_IsUpdated is not sufficient (For example if you
   have multiple subroutines in your program which access the shared memory
   independently. The SPS routines can not distinguish between calls from
   the same subroutine twice or single calls from different subroutines)
   You should keep the old counter value and compare the returned value to
   the old. If it is changed you have to update your array.
   Input: fullname and array
   Returns: update counter value
           -1 Error - memory not longer updated
*/

int SPS_UpdateCounter (char *fullname, char *array); 

/* 
   Tells all the readers of the SPEC array that you updated its contents
   If you are not currently attached to the memory the function will
   attach, write and detach after.
   Input: fullname : Spec version
          array : Array in SPEC version
   Returns: 1 Error - memory not longer updated or no write permission
            0 success
*/

int SPS_UpdateDone (char *fullname, char *array);

/*
  Input: version : name of SPEC version. 
         array_name : Name of this spec array
         rows : Pointer to integer to return number of rows in this array 
         cols : Pointer to integer to return number of cols in this array 
         type : Pointer - Of which data type is the data in this array. 
	        Possible values are:
         flag : Pointer More information about the contents of this array 
	      (does it contain data for MCA, CCD cameras other info .. )
	      Returns: Error code : 0 == no error
*/

int 
SPS_GetArrayInfo (char * spec_version, char * array_name, int *rows, 
		  int *cols, int *type, int *flag);

/* 
   Creates a shared memory array and the shared memory structure for
   the spec version. You can only create new arrays for specversions
   which either do not exist or have been created by this process.
   After the call you are automatically attached to the shared memory.
   The process always stay attached to the shared memories he created.
   You can call all the other functions like Connect or Deconnect but 
   you will always stay attached to the shared memory
   Input: spec_version: The specversion you will create the new array in
          arrayname: The name of the array
	  rows: number of rows 
	  cols: number of columns
	  type : The data type of the data stored in this array.
	  flag: A flag to indicate the type of the array MCA data or 
	        CCD camera.
   Returns: 0 OK
            1 Error 
*/

int
SPS_CreateArray (char * spec_version, char *arrayname,
		 int rows, int cols, int type, int flags) ;

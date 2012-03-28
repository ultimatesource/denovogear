/*  File: arraysub.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1989-
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Arbitrary length arrays
 *              These functions are declared in array.h
 * Exported functions:
 *              See Header file: array.h (includes lots of macros)
 * HISTORY:
 * Last edited: Nov 21 10:54 2008 (rd)
 * Created: Thu Dec 12 15:43:25 1989 (mieg)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "array.h"

/********** Array : class to implement variable length arrays **********/

static int totalAllocatedMemory = 0 ;
static int totalNumberCreated = 0 ;
static int totalNumberActive = 0 ;
static Array reportArray = 0 ;

#define arrayExists(a) ((a) && (a)->magic == ARRAY_MAGIC)

Array uArrayCreate (int n, int size)
{
  Array a = (Array) calloc (1, sizeof (struct ArrayStruct)) ;
  static int isFirst ;

  if (isFirst)
    { isFirst = 0 ;
      if (ARRAY_REPORT_MAX) reportArray = arrayCreate (512, Array) ;
    }
  if (size <= 0)
    { fprintf (stderr, "negative size %d in uArrayCreate\n", size) ; exit (-1) ; }
  if (n < 1)
    n = 1 ;
  totalAllocatedMemory += n * size ;

  a->magic = ARRAY_MAGIC ;
  if (!(a->base = calloc (n, size)))
    { fprintf (stderr, "alloc failure in arrayCreate - requesting %d bytes\n", n*size) ;
      exit (-1) ;
    }
  a->dim = n ;
  a->max = 0 ;
  a->size = size ;
  a->id = ++totalNumberCreated ;
  ++totalNumberActive ;
  if (reportArray)
    { if (a->id < ARRAY_REPORT_MAX)
	array (reportArray, a->id, Array) = a ;
      else
	{ arrayDestroy (reportArray) ;
	  reportArray = 0 ;
	}
    }
  return a ;
}

/**************/

Array uArrayReCreate (Array a, int n, int size)
{
  if (!arrayExists(a))
    return  uArrayCreate (n, size) ;

  if (a->size != size)
    { fprintf (stderr, "type size mismatch in arrayReCreate: size %d != a->size %d\n", size, a->size) ;
      exit (-1) ;
    }

  if (n < 1) n = 1 ;

  if (a->dim < n || (a->dim - n)*size > (1 << 20) ) /* free if save > 1 MB */
    { totalAllocatedMemory -= a->dim * size ;
      free (a->base) ;
      a->dim = n ;
      totalAllocatedMemory += a->dim * size ;
      /* base-mem isn't alloc'd on handle, it's free'd by finalisation */
      if (!(a->base = calloc (n, size)))
	{ fprintf (stderr, "alloc failure in arrayCreate - requesting %d bytes\n", n*size) ;
	  exit (-1) ;
	}
    }
  else
    memset (a->base, 0, n*size) ;

  a->max = 0 ;
  return a ;
}

/**************/

void arrayDestroy (Array a)
{
  if (!arrayExists (a)) { fprintf (stderr, "arrayDestroy called on bad array %x\n", a) ; exit (-1) ; }

  totalAllocatedMemory -= a->dim * a->size ;
  totalNumberActive-- ;
  if (reportArray)
    arr(reportArray, a->id, Array) = 0 ;
  a->magic = 0 ;
  free (a->base) ;
  free (a) ;
}

/**************/

Array arrayCopy (Array a) 
{
  Array new ;

  if (!arrayExists (a)) { fprintf (stderr, "arrayCopy called on bad array %x\n", a) ; exit (-1) ; }
 
  new = uArrayCreate (a->dim, a->size) ;
  memcpy (new->base, a->base, a->dim * a->size) ;
  new->max = a->max ;
  return new;
}

/******************************/

void arrayExtend (Array a, int n) 
{
  char *new ;

  if (!arrayExists (a)) { fprintf (stderr, "arrayExtend called on bad array %x\n", a) ; exit (-1) ; }

  if (n < a->dim)
    return ;

  totalAllocatedMemory -= a->dim * a->size ;
  if (a->dim*a->size < 1 << 23)	/* 8MB */
    a->dim *= 2 ;
  else
    a->dim += 1024 + ((1 << 23) / a->size) ;
  if (n >= a->dim)
    a->dim = n + 1 ;

  totalAllocatedMemory += a->dim * a->size ;

  if (!(new = calloc (a->dim, a->size))) { fprintf (stderr, "calloc failure in arrayExtend\n") ; exit (-1) ; }
  memcpy (new,a->base,a->size*a->max) ;
  free (a->base) ;
  a->base = new ;

  return;
}

/***************/

char *uArray (Array a, int i)
{
  if (!arrayExists (a)) { fprintf (stderr, "array() called on bad array %x\n", a) ; exit (-1) ; }

  if (i < 0) { fprintf (stderr, "referencing array element %d < 0", i) ; exit (-1) ; }

  if (i >= a->max)
    { if (i >= a->dim)
        arrayExtend (a,i) ;
      a->max = i+1 ;
    }
  return a->base + i*a->size ;
}

/***********************************************/
       /* Finds Entry s from Array  a
        * sorted in ascending order of order()
        * If found, returns TRUE and sets *ip
        * if not, returns FALSE and sets *ip one step left
        */

BOOL arrayFind(Array a, void *s, int *ip, int (* order)(void*, void*))
{
  int ord ;
  int i = 0 , j, k ;

  if (!arrayExists (a)) { fprintf (stderr, "arrayFind called on bad array %x\n", a) ; exit (-1) ; }

  j = arrayMax(a) ;
  if (!j || (ord = order(s,uArray(a,0))) < 0)
    { if (ip) *ip = -1 ;
      return FALSE ;
    }   /* not found */

  if (ord == 0)
    { if (ip) *ip = 0 ;
      return TRUE ;
    }

  if ((ord = order(s,uArray(a,--j))) > 0)
    { if (ip) *ip = j ;
      return FALSE ;
    }
  
  if (ord == 0)
    { if (ip) *ip = j ;
      return TRUE ;
    }

  while(TRUE)
    { k = i + ((j-i) >> 1) ; /* midpoint */
      if ((ord = order(s, uArray(a,k))) == 0)
	{ if (ip) *ip = k ;
	  return TRUE ;
	}
      if (ord > 0) i = k ;
      else j = k ;
      if (i == (j-1))
        break ;
    }
  if (ip)
    *ip = i ;
  return FALSE ;
}

/**************************************************************/
       /* Removes Entry s from Array  a
        * sorted in ascending order of order()
        */

BOOL arrayRemove (Array a, void * s, int (* order)(void*, void*))
{
  int i;

  if (!arrayExists (a)) { fprintf (stderr, "arrayRemove called on bad array %x\n", a) ; exit (-1) ; }

  if (arrayFind(a, s, &i,order))
    {
      /* memcpy would be faster but regions overlap
       * and memcpy is said to fail with some compilers
       */
      char *cp = uArray(a,i), *cq = cp + a->size ;
      int j = (arrayMax(a) - i)*(a->size) ;
      while (j--)
	*cp++ = *cq++ ;

      arrayMax(a)-- ;
      return TRUE ;
    }
  else

    return FALSE ;
}

/**************************************************************/
       /* Insert Segment s in Array  a
        * in ascending order of s.begin
        */

BOOL arrayInsert(Array a, void * s, int (*order)(void*, void*))
{
  int i, j, arraySize;

  if (!arrayExists (a)) { fprintf (stderr, "arrayInsert called on bad array %x\n", a) ; exit (-1) ; }

  if (arrayFind(a, s, &i,order))
    return FALSE ;  /* no doubles */
  
  arraySize = arrayMax(a) ;
  j = arraySize + 1 ;

  uArray(a,j-1) ; /* to create space */

	/* avoid memcpy for same reasons as above */
  {
    char *cp, *cq ;
    int k ;
    
    if (arraySize > 0)
      { cp = uArray(a,j - 1) + a->size - 1 ;
	cq = cp - a->size ;
	k = (j - i - 1)*(a->size) ;
	while (k--)
	  *cp-- = *cq-- ;
      }
    
    cp = uArray(a,i+1) ; 
    cq = (char *) s ; 
    k = a->size ;
    while (k--)
      *cp++ = *cq++ ;
  }
  return TRUE ;
}

/**************/

void arrayCompress(Array a)
{
  int i, j, k , as ;
  char *x, *y, *ab ;

  if (!arrayExists (a)) { fprintf (stderr, "arrayCompress called on bad array %x\n", a) ; exit (-1) ; }

  if (arrayMax(a) < 2)
    return ;

  ab = a->base ; 
  as = a->size ;
  for (i = 1, j = 0 ; i < arrayMax(a) ; i++)
    { x = ab + i * as ; y = ab + j * as ;
      for (k = a->size ; k-- ;)		
	if (*x++ != *y++) 
	  goto different ;
      continue ;
      
    different:
      if (i != ++j)
	{ x = ab + i * as ; y = ab + j * as ;
	  for (k = a->size ; k-- ;)	 
	    *y++ = *x++ ;
	}
    }
  arrayMax(a) = j + 1 ;
}

/**************/

int arrayReportMark (void)
{
  if (reportArray)
    return arrayMax (reportArray) ;
  else
    return 0 ;
}

/**************/

void arrayReport (int j)
{
  int i ;
  Array a ;

  fprintf(stderr, "Array report: %d created, %d active, %d MB allocated\n",   
	  totalNumberCreated, totalNumberActive, totalAllocatedMemory/(1024*1024)) ;

  if (reportArray)
    { i = arrayMax (reportArray) ;
      while (i-- && i > j)
	{ a = arr (reportArray, i, Array) ;
	  if (arrayExists(a))
	    fprintf (stderr, " array %d  size %d, max %d\n", i, a->size, a->max) ;
	}
    }
}

/**************/

void arrayStatus (int *nmadep, int *nusedp, int *memAllocp, int *memUsedp)
{ 
  int i ;
  Array a ;

  *nmadep = totalNumberCreated ; 
  *nusedp = totalNumberActive ;
  *memAllocp = totalAllocatedMemory ;
  *memUsedp = 0 ;

  if (reportArray)
    for (i = 0 ; i < arrayMax(reportArray) ; ++i)
      if (arrayExists (a = arr(reportArray, i, Array)))
	*memUsedp += a->max * a->size ;
}

/************************  end of file ********************************/
/**********************************************************************/
 

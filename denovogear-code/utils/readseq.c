/*  File: readseq.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 1993-2006
 *-------------------------------------------------------------------
 * Description: generic code to read Pearson format files (fasta)
 		>header line
		conv[x] is the internal code for char 'x'
		conv[x] == -1 or -2 means ignore. 
                conv[x] < -2 means error.
		will work on fil == stdin
 * Exported functions: readSequence, writeSequence, seqConvert, readFastq, writeFastq
 * HISTORY:
 * Last edited: Oct 22 17:53 2008 (rd)
 * * Oct 22 17:30 2008 (rd): added desc to read/writeFastq
 * * Jul 29 22:42 2008 (rd): added writeFastq
 * * Oct 22 06:45 2006 (rd): added readFastq
 * * Dec 29 23:35 1993 (rd): now works off FILE*, returns id and desc
 * Created: Tue Jan 19 21:14:35 1993 (rd)
 * CVS info: $Id: readseq.c,v 1.3 2001/10/26 09:35:03 searle Exp $
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#pragma inline(add)
static char *messalloc (int n)
{
  char *result ;

  if (!(result = (char*) malloc (n)))
    { fprintf (stderr, "MALLOC failure reqesting %d bytes - aborting\n", n) ;
      exit (-1) ;
    }
  return result ;
}

#define messfree(x) free(x)

static void add (char c, char* *buf, int *buflen, int n)
{
  if (n >= *buflen)
    { 
      int blen = *buflen;
      if (blen < 0)
	{ blen = -blen ;
	  *buf = (char*) messalloc (blen) ;
	}
      else
	{ blen *= 2 ;
	  if ((*buf = realloc(*buf,blen)) == NULL)
             { fprintf (stderr, "REALLOC failure reqesting %d bytes - aborting\n", blen) ;
               exit (-1) ;
             }
	}
      *buflen = blen; 
    }
  (*buf)[n] = c;
}

int readSequence (FILE *fil, int *conv, 
                  char **seq, char **id, char **desc, int *length)
{
  char c ;
  int n ;
  static FILE *oldFil = 0 ;
  static int line ;
  int buflen ;

  if (fil != oldFil)
    { line = 1 ;
      oldFil = fil ;
    }
  
/* get id, descriptor */
  c = getc (fil) ;
  if (c == '>')			/* header line */
    { c = getc(fil) ;

      n = 0 ;			/* id */
      buflen = -32;
      while (!feof (fil) && c != ' ' && c != '\n' && c != '\t')
	{ if (id) add (c, id, &buflen, n++) ;
	  c = getc (fil) ;
	}
      if (id) add (0, id, &buflen, n) ;

				/* white space */
      while (!feof (fil) && (c == ' ' || c == '\t'))
	c = getc (fil) ;

      n = 0 ;			/* desc */
      buflen = -32 ;
      while (!feof (fil) && c != '\n')
	{ if (desc) add (c, desc, &buflen, n++) ;
	  c = getc (fil) ;
	}
      if (desc) add (0, desc, &buflen, n) ;

      ++line ;
    }
  else
    { ungetc (c, fil) ;		/* no header line */
      if (id) 
	*id = "" ;
      if (desc)
	*desc = "" ;
    }

  /* ensure whitespace ignored */

  conv[' '] = conv['\t'] = -1 ;
  conv['\n'] = -3 ;

  n = 0 ;			/* sequence */
  buflen = -1024 ;

  while (!feof (fil))
    { 
      c = getc (fil) ;
      if (c == '>')
	{ ungetc (c, fil) ;
	  break ;
	}

      if (c == EOF || c == EOF + 256) /* satisfies all compilers */
	break ;

    
      switch (conv[c]) {
        case -2:
          { if (id) 
              fprintf (stderr, "Bad char 0x%x = '%c' at line %d, base %d, sequence %s\n",
          	     c, c, line, n, *id) ;
            else
              fprintf (stderr, "Bad char 0x%x = '%c' at line %d, base %d\n",
          	     c, c, line, n) ;
            return 0 ;
          }
          break;
        case -3:
          ++line;
        case -1:
          break;
        default:
	  if (seq) add (conv[c], seq, &buflen, n++) ;
	  else ++n ;
      }
    }
  if (seq) add (0, seq, &buflen, n) ;

  if (length)
    *length = n ;

  return n ;
}

/*****************************************************/

/* a fastq block looks like:

@slxa_0020_2_0001_85
CCTTCTGCCNGCCCCCAGGAAACATTCTCATAAT
+slxa_0020_2_0001_85
???????6%"??2&?7/??'1875%?7:%)1$%2

*/

int readFastq (FILE *fil, int *conv, 
	       char **seq, char **qval, char **id, char **desc, int *length)
{
  char c, *cp ;
  int n, m ;
  static FILE *oldFil = 0 ;
  static int line ;
  int buflen ;

  if (fil != oldFil)
    { line = 1 ;
      oldFil = fil ;
    }
  
/* get id */
  c = getc (fil) ;
  if (!feof(fil) && c != '@')			/* header line */
    { fprintf (stderr, "bad @ header line %d\n", line) ; return 0 ; }
  
  c = getc(fil) ;
  n = 0 ;			/* id */
  buflen = -32;
  while (!feof (fil) && c != ' ' && c != '\n' && c != '\t')
    { if (id) add (c, id, &buflen, n++) ;
      c = getc (fil) ;
    }
  if (id) add (0, id, &buflen, n) ;

				/* white space */
  while (!feof (fil) && (c == ' ' || c == '\t'))
    c = getc (fil) ;

  n = 0 ;			/* desc */
  buflen = -32 ;
  while (!feof (fil) && c != '\n')
    { if (desc) add (c, desc, &buflen, n++) ;
      c = getc (fil) ;
    }
  if (desc) add (0, desc, &buflen, n) ;
  
  if (!feof (fil) && c != '\n')
    { fprintf (stderr, "bad @ identifier line %n\n", line) ; return 0 ; }
  ++line ;

  /* ensure whitespace ignored */

  conv[' '] = conv['\t'] = -1 ;
  conv['\n'] = -3 ;

  n = 0 ;			/* sequence */
  buflen = -1024 ;

  while (!feof (fil))
    { 
      c = getc (fil) ;
      if (c == '+')
	{ ungetc (c, fil) ;
	  break ;
	}

      if (c == EOF || c == EOF + 256) /* satisfies all compilers */
	break ;
    
      switch (conv[c]) {
        case -2:
          { if (id) 
              fprintf (stderr, "Bad char 0x%x = '%c' at line %d, base %d, sequence %s\n",
          	     c, c, line, n, *id) ;
            else
              fprintf (stderr, "Bad char 0x%x = '%c' at line %d, base %d\n",
          	     c, c, line, n) ;
            return 0 ;
          }
          break;
        case -3:
          ++line;
        case -1:
          break;
        default:
	  if (seq) add (conv[c], seq, &buflen, n++) ;
	  else ++n ;
      }
    }
  if (seq) add (0, seq, &buflen, n) ;

  if (length)
    *length = n ;
  
/* second block with same identifier and matching quality values */
  c = getc (fil) ;
  if (!feof(fil) && c != '+')			/* header line */
    { fprintf (stderr, "bad + header line %d\n", line) ; return 0 ; }
  
  c = getc(fil) ; if (id) cp = *id ;
  while (!feof (fil) && c != ' ' && c != '\n' && c != '\t')
    { if (id && c != *cp++) 
	{ fprintf (stderr, "mismatching + identifier line %d\n", line) ; return 0 ; }
      c = getc (fil) ;
    }

  while (!feof (fil) && c != '\n')		/* ignore rest of line - may have desc or not */
    c = getc (fil) ;

  ++line ;

  m = 0 ;			/* qualities */
  buflen = -(n+1) ;

  while (!feof (fil))
    { 
      c = getc (fil) ;
      if (c == '@' && m == n)
	{ ungetc (c, fil) ;
	  break ;
	}
      else if (m > n)
	break ;

      if (c == EOF || c == EOF + 256) /* satisfies all compilers */
	break ;

      if (c == '\n')
	++line;
      else if (qval) 
	add (c-33, qval, &buflen, m++) ;
      else 
	++m ;
    }
  if (qval) add (0, qval, &buflen, m) ;

  if (m != n)
    { fprintf (stderr, "mismatching seq, q length line %n\n", line) ; return 0 ; }

  return n ;
}

/*****************************************************/

int seqConvert (char *seq, int *length, int *conv)
{
  int i, n = 0 ;
  int c ;

  for (i = 0 ; seq[i] ; ++i)
    { c = seq[i] ;
      if (length && i >= *length)
	break ;
      if (conv[c] < -2)
	{ fprintf (stderr, "Bad char 0x%x = '%c' at base %d in seqConvert\n", c, c, n) ;
	  return 0 ;
	}
      if (conv[c] >= 0)
	seq[n++] = conv[c] ;
    }
  if (n < i)
    seq[n] = 0 ;

  if (length)
    *length = n ;
  return n ;
}

/*****************************************************/

int writeSequence (FILE *fil, int *conv, 
		   char *seq, char *id, char *desc, int len)
{
  int i ;

  if (!fil)
    { fprintf (stderr, "ERROR: writeSequence requires a file\n") ;
      return 0 ;
    }
    
  if (id && *id)
    { fprintf (fil, ">%s", id) ;
      if (desc && *desc)
	fprintf (fil, " %s", desc) ;
      fprintf (fil, "\n") ;
    }

  for (i = 0 ; i < len ; ++i)
    { if (i && !(i%60))
	fputc ('\n', fil) ;
      if (conv[seq[i]] > 0)
	fputc (conv[seq[i]], fil) ;
      else
	{ fprintf (stderr, "ERROR in writeSequence: %s[%d] = %d does not convert\n",
		   id, i, seq[i]) ;
	  return 0 ;
	}
    }
  fputc ('\n', fil) ;
  return len ;
}

/*****************************************************/

int writeFastq (FILE *fil, int *conv, 
		char *seq, unsigned char *qual, char *id, char *desc, int len)
{
  int i ;

  if (!fil) { fprintf (stderr, "ERROR: writeFastq requires a file\n") ; return 0 ; }
  if (!id || !*id) { fprintf (stderr, "ERROR: writeFastq requires an id\n") ; return 0 ; }

  if (desc && *desc) 
    fprintf (fil, "@%s %s\n", id, desc) ;
  else
    fprintf (fil, "@%s\n", id) ;
  for (i = 0 ; i < len ; ++i)
    { if (i && !(i%60))
	fputc ('\n', fil) ;
      if (conv[seq[i]] > 0)
	fputc (conv[seq[i]], fil) ;
      else
	{ fprintf (stderr, "ERROR in writeSequence: %s[%d] = %d does not convert\n",
		   id, i, seq[i]) ;
	  return 0 ;
	}
    }
  fputc ('\n', fil) ;
    
  fprintf (fil, "+%s\n", id) ;
  for (i = 0 ; i < len ; ++i)
    { if (i && !(i%60))
	fputc ('\n', fil) ;
      fputc (qual[i]+33, fil) ;
    }
  fputc ('\n', fil) ;

  return len ;
}

/*********** read a matrix, using conv ************/

int readMatrix (char *name, int *conv, int** *mat)
{
  char matdirname[256] ;
  char fullname[512] ;
  FILE *fil ;
  char line[1024] = "#", *p;
  int i, j, nsymb, smax = 0 ;
  int symb[128] ;
  extern char* strtok (char*, const char*) ;

  if (getenv ("BLASTMAT")) 
    strcpy (matdirname, getenv ("BLASTMAT")) ;
  else
    strcpy (matdirname, "/nfs/disk100/pubseq/blastdb/") ;
  strcpy (fullname, matdirname) ;
  strcat (fullname, name) ;

  if (!(fil = fopen (name, "r")) && !(fil = fopen (fullname, "r")))
    { fprintf (stderr, "ERROR in readMatrix: could not open %s or %s\n",
	       name, matdirname) ;
      return 0 ;
    }
    
  while (!feof(fil) && *line == '#') /* comments */
    fgets (line, 1023, fil) ;

				/* character set */
  p = line ; while (*p == ' ' || *p == '\t' || *p == '\n') ++p ;
  for (i = 0 ; *p && i < 128 ; ++i)
    { symb[i] = conv[*p] ;
      if (symb[i] < -2)
	{ fprintf (stderr, "ERROR in readMatrix: illegal symbol %c\n", *p) ;
	  fclose (fil) ;
	  return 0 ;
	}
      if (symb[i] > smax)
	smax = symb[i] ;
      ++p ; while (*p == ' ' || *p == '\t' || *p == '\n') ++p ;
    }
  nsymb = i ;

  ++smax ;
  *mat = (int**) messalloc (smax * sizeof (int*)) ;
  for (i = 0 ; i < smax ; ++i)
    (*mat)[i] = (int*) messalloc (smax * sizeof (int)) ;

  for (i = 0 ; fgets(line, 1023, fil) && i < nsymb ; ++i)
    { p = line ; while (*p == ' ' || *p == '\t' || *p == '\n') ++p ;
      if (p && conv[*p] == symb[i])
	{ ++p ; while (*p == ' ' || *p == '\t' || *p == '\n') ++p ; }
      for (j = 0 ; *p && j < 128 ; ++j)
	{ if (symb[i] >= 0 && symb[j] >= 0) 
	    (*mat)[symb[i]][symb[j]] = atoi(p) ;
	  if (*p == '-') ++p ;
	  while (p && *p >= '0' && *p <= '9') ++p ;
	  while (*p == ' ' || *p == '\t' || *p == '\n') ++p ;
	}
      if (j != nsymb)
	{ fprintf (stderr, "ERROR in readMatrix: bad line: %s\n", line) ;
	  fclose (fil) ;
	  return 0 ;
	} 
    }

  fclose (fil) ;
  return 1 ;
}

/*********** standard conversion tables **************/

int dna2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, '-',  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbig2NConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int aa2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
} ;

int aa2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
} ;

int noConv[] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
} ;
/**************** end of file ***************/

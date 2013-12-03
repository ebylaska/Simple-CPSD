/*
*************************************************************************
*									*
* compressed_io.c							*
*									*
* These routines implement a simple input/output package for writing to	*
* and reading from gzipped files.					*
*                                                                       *
*									*
* Fortran usage:							*
*	call openfile(unit,'file', 'X', n)	X is either 'r' or 'w'	*
*					n is the length of the filename	*
*	call cwrite(unit,c, n)		c is n element character array	*
*	call iwrite(unit,i, n)		i is n element integer array	*
*	call dwrite(unit,d, n)		d is n element double array	*
*	call cread(unit,c, n)		c is n element character array	*
*	call iread(unit,i, n)		i is n element integer array	*
*	call dread(unit,d, n)		d is n element double array	*
*	call closefile(unit)		close the datafile		*
*									*
* Author:      Scott Kohn (skohn@chem.ucsd.edu)				*
* modified by: Eric Bylaska (ebylaska@chem.ucsd.edu)			*
*									*
*************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**** You need to define the locations of gzip and zcat ****/
//  #define GZIP  "/bin/gzip"
//  #define GZCAT "/bin/zcat"

#define GZIP  ""
#define GZCAT ""



#define FORTRAN_NAME(X) X
extern FILE *popen(const char *, const char *);


#define MAX_UNIT	10

static FILE* fd[MAX_UNIT];	/* the file descriptor of the pipe */

#define BAIL(X) { fprintf(stderr, X); exit(-1); }

/*
*************************************************************************
*									*
* Define the Xwrite and Xread routines using the Fortran bindings.	*
*									*
*************************************************************************
*/

void FORTRAN_NAME(cwrite_)(const int *unit, 
                          const char *c,
                          const int *n)
{
   (void) fwrite(c, sizeof(char), *n, fd[*unit]);
}

void FORTRAN_NAME(cread_)(const int *unit, 
                               char *c, 
                         const int *n)
{
   (void) fread(c, sizeof(char), *n, fd[*unit]);
}

void FORTRAN_NAME(iwrite_)(const int *unit, const int *i, const int *n)
{
   (void) fwrite(i, sizeof(int), *n, fd[*unit]);
}

void FORTRAN_NAME(iread_)(const int *unit, int *i, const int *n)
{
   (void) fread(i, sizeof(int), *n, fd[*unit]);
}

void FORTRAN_NAME(dwrite_)(const int *unit, const double *d, const int *n)
{
   (void) fwrite(d, sizeof(double), *n, fd[*unit]);
}

void FORTRAN_NAME(dread_)(const int *unit, double *d, const int *n)
{
   (void) fread(d, sizeof(double), *n, fd[*unit]);
}

/*
*************************************************************************
*									*
* void openfile(char *filename, char *mode, int *n)			*
* void closefile()							*
*									*
* Function openfile opens a pipe to either gzip (to compress a stream)	*
* or zcat (to uncompress a stream).  Function closefile() closes the	*
* pipe stream created by openfile().					*
*									*
*************************************************************************
*/

#define FUDGE_FACTOR (8)

void FORTRAN_NAME(openfile_)(const int *unit, 
                                  char *filename, 
       				  char *mode,    
				  int *n)
{
   //const int buffersize = strlen(GZIP)+strlen(GZCAT)+(*n)+FUDGE_FACTOR;
   //char *command = (char *) malloc(buffersize);
   char *file = (char *) malloc(*n+1);

   (void) strncpy(file, filename, *n);
   file[*n] = '\0';

   if ((*mode == 'r') || (*mode == 'R')) {
      if (!(fd[*unit] = fopen(filename, "r")))
         BAIL("ERROR:  Could not open pipe from input file\n");
   } else {
      if (!(fd[*unit] = fopen(filename, "w")))
         BAIL("ERROR:  Could not open pipe to output file\n");
   }

   free(file);
   //free(command);
}

void FORTRAN_NAME(closefile_)(const int *unit)
{
   (void) fclose(fd[*unit]);
}


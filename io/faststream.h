#ifndef INCLUDED_FASTSTREAM_H
#define INCLUDED_FASTSTREAM_H 1
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpi.h"

typedef struct TextFileData *TextFile;

/* Input routines */
int TextFileReadOpen(const char filename[], MPI_Comm comm, MPI_Info info,
		     TextFile *fh_p);
int TextFileReadClose(TextFile fh);
int TextFileReadFree(TextFile *fh_p);
int TextFileReadIterate(TextFile fh,
			int (*processline)(const char *, int len, void *),
			void *extra_data);
int TextFileReadPrintStats(TextFile fh, FILE *fp);
int TextFileReadGetStats(TextFile fh, int tsize, double t[],
			 MPI_Offset *bytesRead);
int TextFileReadGetInfo(TextFile fh, char *str, int *len);
int TextFileReadDebug(int);

/* Output routines */
int TextFileWriteOpen(MPI_Comm comm, int root, const char *fname, MPI_Info info,
		      TextFile *fh_p);
int TextFileWriteOpenFromFH(MPI_Comm comm, int root, FILE *fp, MPI_Info info,
			    TextFile *fh);
int TextFileWriteOrdered(TextFile fh, const char *buf, size_t len);
int TextFileWritePrintf(TextFile fh, const char *restrict format, ...);
int TextFileWriteFlush(TextFile fh);
int TextFileWriteClose(TextFile fh);
int TextFileWriteFree(TextFile *fh_p);
int TextFileWritePrintStats(TextFile fh, FILE *fp);
int TextFileWriteGetStats(TextFile fh, int *tsize, double t[],
			  MPI_Offset *bytesWritten);
int TextFileWriteNameStats(int idx, int len, char label[]);

#endif

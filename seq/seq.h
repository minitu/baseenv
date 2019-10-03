#ifndef SEQ_H_INCLUDED
#define SEQ_H_INCLUDED 1
void seqBegin(MPI_Comm comm);
void seqEnd(MPI_Comm comm);
void seqChangeOrder(MPI_Comm, int, int);
#endif

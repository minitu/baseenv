#ifndef NODECART_H_INCLUDED
#define NODECART_H_INCLUDED

/* These routines are replacements for MPI_Cart_xxx routines */
int MPIX_Nodecart_create(MPI_Comm comm_old, int ndims, const int dims[],
			 const int periods[], int reorder,
			 MPI_Comm *comm_cart);
int MPIX_Nodecart_shift(MPI_Comm comm, int direction, int disp,
			int *rank_sources, int *rank_dest);
int MPIX_Nodecart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
int MPIX_Nodecart_rank(MPI_Comm comm, const int coords[], int *rank);
int MPIX_Nodecart_sub(MPI_Comm comm, const int remain[], MPI_Comm *subcomm);
int MPIX_Nodecart_get(MPI_Comm comm, int maxdims, int dims[], int periods[],
		      int coords[]);
int MPIX_Nodecart_dim_get(MPI_Comm comm, int *ndims);
/* End of routines that can replace MPI_Cart_xxx routines */

int MPIX_Nodecomm_create(MPI_Comm incomm, MPI_Comm *nodecomm,
			 MPI_Comm *leadercomm, int *noderank, int *nnodes);

int MPIX_Comm_dims_create(MPI_Comm comm, int nsize, int ndims, int dims[]);

/* Utility routine for information about the node */
int MPIX_GetSocketAndCPU(int *nsocket, int *socketrank, int *rankinsocket);
int MPIX_GetNodeTopoInfo(MPI_Comm comm, int *nnodes, int *noderank,
			 int *nodesize, int *rankonnode,
			 int *nchips, int *chiprank, int *rankinchip);

/* Utility routine based on having node decomposition */
int MPIX_NodeSubset(MPI_Comm incomm, MPI_Comm leadercomm, MPI_Comm nodecomm,
		    int nnodes, MPI_Comm *outcommPtr);

/* Utility routine to get information on decomposition for inter/intra node */
int MPIX_Nodecart_get_decomp(MPI_Comm comm, int *nlevel, int *ndim, int *dims);

/* Utility routine to provide cvar-like access */
void MPIX_Nodecart_cvar_set(const char *, int value);

/* Utility routine for provinding information output */
int MPIX_PrintNodeCommCounts(FILE *fp, MPI_Comm comm, int nr, const int ranks[],
			     MPI_Comm nodecomm);
int MPIX_PrintNodeCommCounts_X(FILE *fp, MPI_Comm comm, int nr, const int ranks[],
			       MPI_Comm nodecomm, MPI_Comm socketcomm);
int MPIX_PrintNodeDecomp(FILE *fp, MPI_Comm comm);
int MPIX_GetOffNodeCounts(MPI_Comm comm, int nr, const int ranks[],
			  MPI_Comm nodecomm,
			  int *onnode, int *offnode);
int MPIX_GetOffNodeCounts_X(MPI_Comm comm, int nr, const int ranks[],
			    int ncomm, MPI_Comm subcomm[],
			    int *onnode, int *offnode);
int MPIX_GetNestedCounts(MPI_Comm comm, MPI_Comm nodecomm, MPI_Comm socketcomm,
			 int nr, const int ranks[],
			 int *onsock, int *onnode, int *offnode);


/* Temp routine to enable use of socket info */
int MPIX_SetUseSocket(int flag);
int MPIX_Socketcomm_create(MPI_Comm nodecomm, MPI_Comm *socketcomm);

#endif /* NODECART_H_INCLUDED */

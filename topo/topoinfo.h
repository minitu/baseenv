#ifndef TOPOINFO_H_INCLUDED
#define TOPOINFO_H_INCLUDED

#define MAX_TOPO_DIM 5
typedef enum { TOPO_UNKNOWN, TOPO_NODELIST,
	       TOPO_TORUS, TOPO_MESH, TOPO_BUS, TOPO_TREE,
	       TOPO_NODE, TOPO_SOCKET, TOPO_CORE }
	topoType_t;

/*E
  topoMethod_t - Available methods for determining topology
  E*/
typedef enum { TOPO_DEFAULT, TOPO_GENERIC, TOPO_CRAY_RCA, TOPO_BGQ_MPIX,
	       TOPO_MPISPLIT, TOPO_NODENAME,
	       TOPO_HWLOC, TOPO_GETCPU, TOPO_GETAFFINITY,
	       TOPO_DEBUG } topoMethod_t;

/*S
  topoCoord_t - Structure containing the coordinates for a topology
  S*/
typedef struct { int coords[MAX_TOPO_DIM]; } topoCoord_t;

/*S
  topoentry_t - Structure containing information about a level of the system topology

 Notes:
 This structure may be examined by the user to determine properties of the
 system topology.  Most systems will have at least two levels, with 'topoType'
 'TOPO_CORE' and 'TOPO_CHIP'.  Additional levels may include 'TOPO_NODE' and
 for systems with a high performance mesh interconnect, 'TOPO_MESH'.  Other
 topology types are defined.

 Each level may have siblings (e.g., cores on the same chip), children
 (e.g., cores of the chip), and parents (e.g., the node hosting a chip).
 Coordinates are provided to identify siblings.  For simple arrangements,
 such as cores on a chip, this is just a numbering starting from zero;
 specifically, a one-tuple.  For a 'TOPO_MESH' or 'TOPO_TORUS' topology,
 the coordinates will be an n-tuple, for an n-Dimensional mesh.

 The source code for the routine 'topoPrint' shows how to use some of the
 fields in this structure.
  S*/
typedef struct topoentry_t {
    topoType_t  topoType;      /* e.g., unknown, torus, mesh, bus, tree */
    int         dim;           /* number of dimensions in topology */
    topoCoord_t coords;        /* coordinates in topology */
    topoCoord_t maxtopocoords; /* maximum value of valid coordinate in
				  topology */
    topoCoord_t maxcoords;     /* maximum value of valid coordinates */
    topoCoord_t mincoords;     /* minimum value of valid coordinates */
    struct topoentry_t *next;
} topoentry_t;

typedef struct {
    int         numLevels;     /* */
    /* The node is a special thing.  Maintain information just about the
       nodes */
    int         nranks;        /* Number of processes on this node */
    int         *ranks;        /* Ranks (in MPI_COMM_WORLD) of processes
				  on the same node */
    int         nodenum;       /* Number of node, from 0..nnodes-1 */
    int         nnodes;        /* Number of nodes.  If -1, node info has
				  not been initialized */
    MPI_Comm    nodeComm,      /* Communicator of processes on the same node */
	nodeLeaderComm;        /* Communicator of leaders of the nodes */
    topoentry_t *info, *infoTail;
} topoinfo_t;

int topoInit(int, topoinfo_t **);
int topoDebug(int);
int topoFinalize(topoinfo_t **);
int topoPrint(FILE *, const char *, topoinfo_t *);
int topoToStr(topoinfo_t *, int, char *, int);

int topoNodeInfoBasic(topoinfo_t *, int *, int *, int *, int *);
int topoMeshCoords(topoinfo_t *, int *, int [], int []);
int topoToArray(topoinfo_t *, int [], int [], int *, int *, int);
int topoGetRingFromMesh(topoinfo_t *topo, int ndim, const int meshcoords[],
			const int qtorus[], int *nodenum);
int topoMeshContainer(topoinfo_t *topo, int *ndim,
		      int mindim[], int maxdim[], int qtorus[]);
int TopoFindPartnerCliquenum(int nranks, const int ranks[], MPI_Comm comm,
			     int myClique, int cliqueArray[]);
int TopoComputeRanksSameClique(int nranks, const int ranks[], MPI_Comm comm,
			       int myClique, int *nSameClique );

/* Data structures and routines for getting distance information */

#define MAX_TOPO_TOTAL_DIM 10
/* This contains information on the coordinates of a process with the
   noted rank */
typedef struct topoarray_t {
    int rank;                          /* Rank of the process */
    int nlevels;                       /* Number of levels in topology */
    int nodeidx;                       /* Index of node (see topoToArray) */
    int mycoords[MAX_TOPO_TOTAL_DIM];  /* coordinates at each level */
    int maxcoords[MAX_TOPO_TOTAL_DIM]; /* Maximum coordinates at each level */
} topoarray_t;

typedef struct topodist_t {
    int         ntargets;        /* Number of processes sent to (number of
				    elements in tarray */
    topoarray_t myarray;
    topoarray_t *tarray;         /* Array of coordinates for each target
				    process */
} topodist_t;

int topodistInit(MPI_Comm comm, int nsend, int sendrank[],
		 int nrecv, int recvrank[],
		 topoinfo_t *ti, topodist_t **td);
void topodistFree(topodist_t *dt);
int topodistNodeCommInfo(topodist_t *td, MPI_Comm comm,
			 int *partnerOnNode, int *partnerOffNode,
			 MPI_Comm *nodecomm_p, int *totalOffNode);
int topoMeshHopDistance(topodist_t *dt, int rank, int *dist);

int topoAvailMethods(int maxnum, int *nnode, int *nnet,
		     topoMethod_t nodeMethod[], topoMethod_t netMethod[]);
int topoSetMethods(topoMethod_t nodeMethod, topoMethod_t netMethod);
int topoGetMethodsDesc(const char **nodestr, const char **netstr);
int topoNodeEnumeration(topoinfo_t *topo, int *numnodes, int *mynodenum,
			int *nranks, int noderanks[]);
const char *topoTypeStr(topoType_t t);
#endif

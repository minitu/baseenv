#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "topoinfo.h"
#include "seq.h"

#define MAX_DIM 6
#define MAX_METHODS 10

void topodistTest(FILE *);
int topodistPrint(topodist_t *, FILE *);

int main(int argc, char **argv)
{
  topoinfo_t *topoinfo;
  int  wrank, verbose=0, nodenum, rc;
  char leader[10];
  char topostr[256];
  int  ndim=MAX_DIM, coords[MAX_DIM], qtorus[MAX_DIM];
  int  mycoords[MAX_DIM], maxcoords[MAX_DIM], ncoords, nodeidx, i;
  topoMethod_t nodemeths[MAX_METHODS], netmeths[MAX_METHODS];
  int          n1, n2, nnodemeths, nnetmeths;
  int          debugtopo = 0;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&wrank);
  for (i=1; i<argc; i++) {
      if (argv[i]) {
	  if (strcmp(argv[i], "-v") == 0) verbose = 1;
	  else if (strcmp(argv[i], "-d") == 0) debugtopo = 1;
	  else {
	      if (wrank == 0) {
		  fprintf(stderr, "Unrecognized argument %s\n", argv[i]);
		  fflush(stderr);
	      }
	      MPI_Abort(MPI_COMM_WORLD, 1);
	  }
      }
  }

  /* Test all available methods in combinations */
  topoAvailMethods(MAX_METHODS, &nnodemeths, &nnetmeths, nodemeths, netmeths);
  for (n1=0; n1<nnodemeths; n1++) {
      if (!debugtopo && nodemeths[n1] == TOPO_DEBUG) continue;
      for (n2=0; n2<nnetmeths; n2++) {
	  const char *nodestr, *netstr;
	  if (!debugtopo && netmeths[n2] == TOPO_DEBUG) continue;
	  topoSetMethods(nodemeths[n1], netmeths[n2]);
	  topoGetMethodsDesc(&nodestr, &netstr);
	  if (wrank == 0) {
	      printf("Node method: %s\nNet method %s\n", nodestr, netstr);
	      printf("Rank:Type:Dimension:Coordinates\n");
	      fflush(stdout);
	  }

	  snprintf(leader,sizeof(leader),"%d:",wrank);
	  topoInit(verbose,&topoinfo);
	  MPI_Barrier(MPI_COMM_WORLD);
          seqBegin(MPI_COMM_WORLD);
	  topoPrint(stdout,leader,topoinfo);
          fflush(stdout);
          seqEnd(MPI_COMM_WORLD);

	  topoToStr(topoinfo,1,topostr,sizeof(topostr));
	  rc = topoMeshCoords(topoinfo, &ndim, coords, qtorus);
	  if (rc == 0) {
	      topoGetRingFromMesh(topoinfo, ndim, coords, qtorus, &nodenum);
	      if (topostr[0])
		  printf("%d:%s:ringnode %d\n", wrank, topostr, nodenum);
	  }
	  MPI_Barrier(MPI_COMM_WORLD);

	  /* Now, check the conversion of topo information into the array used
	     in some of the other tools, such as nodecomm */
	  if (wrank == 0) {
	      printf("Topo To array; Node elements (one line per level) then interconnect\n");
	      fflush(stdout);
	  }
	  topoToArray(topoinfo, mycoords, maxcoords, &ncoords, &nodeidx,
		      MAX_DIM);
	  seqBegin(MPI_COMM_WORLD);
	  printf("%d:Total coords: %d Nodeidx: %d\n",
		 wrank, ncoords, nodeidx);
	  for (i=0; i<nodeidx; i++) {
	      printf("%d:Node:[%d] %d in %d\n",
		     wrank, i, mycoords[i], maxcoords[i]);
	  }
	  for (i=nodeidx; i<ncoords; i++) {
	      printf("%d:Interconnect:[%d] %d in %d\n",
		     wrank, i, mycoords[i], maxcoords[i]);
	  }
	  fflush(stdout);
	  seqEnd(MPI_COMM_WORLD);

	  MPI_Barrier(MPI_COMM_WORLD);

	  topodistTest(stdout);

	  MPI_Barrier(MPI_COMM_WORLD);
	  {
	      int numnodes, mynodenum, nranks, noderanks[128];
              int err;
	      nranks = 128;
	      err = topoNodeEnumeration(topoinfo, &numnodes, &mynodenum,
		                        &nranks, noderanks);
              if (wrank == 0) {
	          if (err == 0) {
		      printf("Node Enumeration: with %d nodes of %d processes\n",
                             numnodes, nranks);
                  }
                  else {
                      printf("Node Enumeration failed.\n");
                      if (n2 == TOPO_GENERIC)
                          printf("Expected with generic net\n");
                  }
		  fflush(stdout);
              }
	      seqBegin(MPI_COMM_WORLD);
              if (err == 0) {
	          printf("%d: mynode %d: ", wrank, mynodenum);
	          for (i=0; i<nranks; i++) {
		      printf("%d%s", noderanks[i], (i < nranks-1) ? ",":"");
	          }
	          printf("\n");
	          fflush(stdout);
              }
	      seqEnd(MPI_COMM_WORLD);
	  }

	  topoFinalize(&topoinfo);
      }
  }

  MPI_Finalize();
  return 0;
}

#include <string.h>
int topodistPrint(topodist_t *td, FILE *fp)
{
    int nOffNode,
	nOnNode,
	totalOffNode,
	nodeOnNode,
	nodeOffNode,
	nrank, nsize, rank, size,
	err;
    MPI_Comm nodecomm;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
	printf("Test of on/off node determination\n");
	fflush(stdout);
    }

    err = topodistNodeCommInfo(td, MPI_COMM_WORLD, &nOnNode, &nOffNode,
			       &nodecomm, &totalOffNode);

    MPI_Allreduce(&nOnNode, &nodeOnNode, 1, MPI_INT, MPI_SUM, nodecomm);
    MPI_Allreduce(&nOffNode, &nodeOffNode, 1, MPI_INT, MPI_SUM, nodecomm);

    MPI_Comm_rank(nodecomm, &nrank);
    MPI_Comm_size(nodecomm, &nsize);
    if (nrank == 0) {
	printf("rank = %d, nsize = %d, Sends: Node On = %d, Off = %d; procs sending off = %d\n",
	       rank, nsize, nodeOnNode, nodeOffNode, totalOffNode);
	fflush(stdout);
    }
    MPI_Comm_free(&nodecomm);
    MPI_Barrier(MPI_COMM_WORLD);

    {
	char coordStr[1024];
	int j;
	coordStr[0] = 0;
	for (j=td->myarray.nodeidx; j<td->myarray.nlevels; j++) {
	    char numstr[20];
	    snprintf(numstr, 20, "%d", td->myarray.mycoords[j]);
	    strncat(coordStr, numstr, 1024-1);
	    strncat(coordStr, ",", 1024-1);
	}
	if (rank > 0) {
	    MPI_Recv((void*)0, 0, MPI_INT, rank-1, 1, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	}
	printf("%d:[%s] OnNode = %d, OffNode = %d\n", rank, coordStr, nOnNode,
	       nOffNode); fflush(stdout);
	if (rank < size-1) {
	    MPI_Send((void*)0, 0, MPI_INT, rank+1, 1, MPI_COMM_WORLD);
	}
    }
    return 0;
}

#include <math.h>
void  topodistTest(FILE *fp)
{
    topodist_t *dt;
    topoinfo_t *ti;
    int        rank, size, sends[4], m, ns;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    m = sqrt((double)size);
    ns = 0;
    if (rank > 0)           sends[ns++] = rank - 1;
    if (rank < size-1)      sends[ns++] = rank + 1;
    if (rank - m >= 0)      sends[ns++] = rank - m;
    if (rank + m <= size-1) sends[ns++] = rank + m;
    topoInit(1,&ti);
    topodistInit(MPI_COMM_WORLD, ns, sends, ns, sends, ti, &dt);
    topodistPrint(dt, fp);
    topodistFree(dt);
}

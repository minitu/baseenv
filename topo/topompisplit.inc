/* -*- Mode: C; c-basic-offset:4 ; -*- */
static int topoiGetNodeTopoInfoMPISPLIT(topoinfo_t *topo);

static int topoiGetNodeTopoInfoMPISPLIT(topoinfo_t *topo)
{
    topoentry_t *e;
    MPI_Comm nodecomm, leadercomm;
    int      nrank, nsize, color, nodenums[2];
    int      i, *range;
    MPI_Group gn, gw;

    /* Note that with a constant keys, ties are broken wrt the
       input communicator/group */
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
			MPI_INFO_NULL, &nodecomm);

    /* We need to also enumerate the different comms (nodes).
       First, determine how many there are */
    MPI_Comm_size(nodecomm, &nsize);
    MPI_Comm_rank(nodecomm, &nrank);
    color = 0;
    if (nrank != 0) color = MPI_UNDEFINED;
    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &leadercomm);
    if (color == 0) {
	MPI_Comm_rank(leadercomm, &nodenums[0]);
	MPI_Comm_size(leadercomm, &nodenums[1]);
	MPI_Comm_free(&leadercomm);
    }
    MPI_Bcast(&nodenums, 2, MPI_INT, 0, nodecomm);

#if 0
    e = topoiAllocEntry(topo);
    if (!e) return 1;
    e->dim = 1;
    e->maxcoords.coords[0] = nsize;
    e->topoType            = TOPO_NODE;
    e->coords.coords[0]    = nrank;
    e->coords.coords[1]    = -1;
#endif

    
    e = topoiAllocEntry(topo);
    if (!e) return 1;
    e->dim = 1;
    e->maxcoords.coords[0] = nodenums[1];
    e->topoType            = TOPO_NODELIST;
    e->coords.coords[0]    = nodenums[0];
    e->coords.coords[1]    = -1;

    /* Compute the ranks of the processes on the same node */
    if (topo->nnodes == -1) {
	topo->nnodes  = nodenums[1];
	topo->nodenum = nodenums[0];
	topo->nranks  = nsize;
	topo->ranks   = (int *)malloc(nsize * sizeof(int));
	range         = (int *)malloc(nsize * sizeof(int));
	if (!range || !topo->ranks) return 1;

	for (i=0; i<nsize; i++) range[i] = i;
	MPI_Comm_group(nodecomm, &gn);
	MPI_Comm_group(MPI_COMM_WORLD, &gw);
	MPI_Group_translate_ranks(gn, nsize, range, gw, topo->ranks);
	MPI_Group_free(&gn);
	MPI_Group_free(&gw);
	free(range);
    }

#if 0
    e = topoiAllocEntry(topo);
    if (!e) return 1;
    e->dim = 1;
    e->maxcoords.coords[0] = nsize;
    e->topoType            = TOPO_CORE;
    e->coords.coords[0]    = nrank;
    e->coords.coords[1]    = -1;

    e = topoiAllocEntry(topo);
    if (!e) return 1;
    e->dim = 1;
    e->maxcoords.coords[0] = nodenums[1];
    e->topoType            = TOPO_NODE;
    e->coords.coords[0]    = nodenums[0];
    e->coords.coords[1]    = -1;
#endif

    MPI_Comm_free(&nodecomm);
    if (leadercomm != MPI_COMM_NULL)
	MPI_Comm_free(&leadercomm);


    return 0;
}

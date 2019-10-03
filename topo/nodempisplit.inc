/* -*- Mode: C; c-basic-offset:4 ; -*- */
static int topoiGetNodeInfo(int isMultithreaded, topoinfo_t *topo)
{
    topoentry_t *e;
    MPI_Comm nodecomm, leadercomm;
    int      nrank, nsize, color, nodenums[2];
    MPI_Group g1, gw;

    /* Note that with a constant keys, ties are borken wrt the
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

    MPI_Comm_size(nodecomm, &e->nranks);

    MPI_Comm_group(nodecomm, &g1);
    MPI_Comm_group(MPI_COMM_WORLD, &gw);
    e->ranks = (int *)malloc(e->nranks * sizeof(int));
    MPI_Group_translate_ranks(g1, e->nranks, gw, e->ranks);
    MPI_Group_free(&g1);
    MPI_Group_free(&gw);

    MPI_Comm_free(&nodecomm);

    return 0;
}

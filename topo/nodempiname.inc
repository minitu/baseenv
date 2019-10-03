/* -*- Mode: C; c-basic-offset:4 ; -*- */
#error 'Not yet implemented'

static int topoiGetNodeInfo(int isMultithreaded, topoinfo_t *topo)
{
    topoentry_t *e;
    MPI_Comm nodecomm, leadercomm;
    int      nrank, nsize, color, nodenums[2];

    err = gethostname(myname, sizeof(myname));
    if (err) {
	perror( "gethostname: " );
	MPI_Abort( MPI_COMM_WORLD, 1 );
    }
    err = MPE_FindCliqueFromName( myname, comm, 16, MAX_CLIQUE_SIZE, 
				  cliqueNum, cliqueSize, cliqueRanks );
    
    ---- old code ----
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

    MPI_Comm_free(&nodecomm);

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

    return -1;
}

/*@
  MPE_FindCliqueFromName - Given a name that is unique to a node, find the
                           ranks of processes on the same node

  Input Parameters:
+ name - name of the node on which the process is running
. comm - communicator 
. hintSize - Used to control initial hash size
- maxSize - Maximum number of ranks in cliqueRanks

  Output Parameters:
+ cliqueNum - Id of this clique, starting from zero
. cliqueSize - Number of processes in this clique
- cliqueRanks - Ranks of processes in this clique

  Note:
  This is the general routine.  A number of convenience routines
  are provided to create communicators from cliques, and to return
  cliques in a communicator by using the node name.
  @*/
int MPE_FindCliqueFromName( const char *name, MPI_Comm comm, 
			    int hintSize, int maxSize, 
			    int *cliqueNum, int *cliqueSize, int cliqueRanks[] )
{
    int         csize, hsize, i, hash, goalVal, goalMaxVal, maxVal;
    int         rank, namelen, maxlen, *hashArray = 0;
    const char *p = 0;
    int        *cliqueInfo = 0;
    MPI_Status  status;
    MPI_Request req;

    /* Initial hash size and create hash array */
    ttime = MPI_Wtime();
    MPI_Comm_size( comm, &csize );
    if (csize == 1) {
	*cliqueSize = 1;
	*cliqueNum  = 1;
	cliqueRanks[0] = 0;
	return 0;
    }
    MPI_Comm_rank( comm, &rank );

    hsize = (hintSize + csize) / hintSize;
    hashArray = (int *) calloc( hsize, sizeof(int) );
    if (!hashArray) {
	return errReturn( "Unable to allocate hash table" );
    }

    goalVal    = (3 + hintSize) / 4;
    goalMaxVal = hintSize * 4;
    do {
	int nEmpty;
	/* Compute hash */
	p    = name;
	hash = 0;
	while (*p) hash += (*p++ - ' ');
	hash = hash % hsize;
	hashArray[hash] = 1;
	
	MPI_Allreduce( MPI_IN_PLACE, hashArray, hsize, MPI_INT, MPI_SUM, comm );

	/* */
	maxVal = 0;
	nEmpty = 0;
	for (i=0; i<hsize; i++) {
	    if (hashArray[i] > maxVal) maxVal = hashArray[i];
	    else if (hashArray[i] == 0) nEmpty++;
	}
	/* hashOccupancy is % occupancy. */
	hashOccupancy = (nEmpty * 100)/hsize;
	
	if (maxVal < goalVal && hsize > 16) hsize = hsize / 2;
	else if (maxVal > goalMaxVal && hsize < csize / 4) {
	    free( hashArray );
	    hsize *= 2;
	    hashArray = (int *)calloc( hsize, sizeof(int) );
	    if (!hashArray) {
		return errReturn( "Unable to reallocate hash table" );
	    }
	}
	else
	    /* Good enough */
	    break;
    } while (1);

    if (verbose && rank == 0) printf( "[0] Hash size = %d\n", hsize );

    /* At this point, we know that at most maxVal messages will be sent to any
       process */
    namelen = strlen( name ) + 1;
    MPI_Allreduce( &namelen, &maxlen, 1, MPI_INT, MPI_MAX, comm );

    /* Send to my target the info */
    /* Should avoid sending to myself */
    MPI_Isend( (char *)name, namelen, MPI_CHAR, hash, 0, comm, &req );
    if (verbose) printf( "[%d] sent to %d\n", rank, hash );

    /* Allocate space for the clique info */
    cliqueInfo = (int *)calloc( maxVal + 2, sizeof(int) );
    if (!cliqueInfo) {
	return errReturn( "Unable to allocate clique info" );
    }

    /* If I'm a target, receive data on the names */
    if (rank < hsize && hashArray[rank]) {
	NameInfo *names;
	int nclique, myIndex=0, j;
	char *curName;

	/* Allocate an array for the data */
	names = (NameInfo *) calloc( hashArray[rank], sizeof(NameInfo) );
	if (!names) {
	    /* Really need to panic */
	    return errReturn( "Unable to allocate names array" );
	}
	for (i=0; i<hashArray[rank]; i++) {
	    names[i].name = (char *)calloc( maxlen, sizeof(char) );
	    MPI_Recv( names[i].name, maxlen, MPI_CHAR, MPI_ANY_SOURCE, 0, comm, 
		      &status );
	    names[i].srcRank = status.MPI_SOURCE;
	    if (verbose & 0x2) printf( "[%d] received %s from %d\n", 
				       rank, names[i].name, status.MPI_SOURCE );
	}
	if (verbose) printf( "[%d] Received info\n", rank );
	
	/* Sort by name (note that there may be multiple names that hashed
	   to the same value */
	qsort( names, hashArray[rank], sizeof(NameInfo), nameCompare );

	/* Compute the number of cliques */
	nclique = 0;
	curName = (char *)"\0";
	for (i=0; i<hashArray[rank]; i++) {
	    if (strcmp(names[i].name, curName)) {
		curName = names[i].name;
		nclique++;
	    }
	}

	/* Compute the clique index*/
	MPI_Exscan( &nclique, &myIndex, 1, MPI_INT, MPI_SUM, comm );

	if (verbose) printf( "[%d] Computed clique index = %d\n", rank, myIndex);
	/* Return the information to each sender.  Note that we could 
	   parallelize this step by sending just to the root; let each 
	   root notify its members.
	   Each process is sent an integer array containing:
	   cliqueIndex, cliqueSize, ranks in clique
	*/
	j = 0;
	for (i=0; i<nclique; i++) {
	    int nInClique = 0, k;

	    cliqueInfo[0] = myIndex++;
	    curName       = names[j].name;
	    k             = j;
	    for (;j<hashArray[rank] && strcmp(curName,names[j].name) == 0; j++) {
		cliqueInfo[2+j-k] = names[j].srcRank;
		if (verbose & 0x2) printf( "[%d] Adding src %d to clique[%d]\n", 
				     rank, names[j].srcRank, j-k );
		nInClique++;
	    }
	    cliqueInfo[1] = nInClique;
	    if (verbose) printf( "[%d] clique size = %d\n", rank, nInClique );
	    for (k=0; k<nInClique; k++) {
		if (verbose & 0x2) printf( "[%d] Sending back to %d\n", rank, 
				     cliqueInfo[2+k] );
		if (rank == cliqueInfo[2+k]) {
		    int kk;
		    *cliqueNum = cliqueInfo[0];
		    *cliqueSize = cliqueInfo[1]; 
		    if (*cliqueSize > maxSize) {
			if (verbose) {
			    fprintf( stderr, "Cliquesize of %d too large\n", 
				     *cliqueSize );
			}
			MPI_Abort( MPI_COMM_WORLD, 1 );
			return 1;
		    }
		    for (kk=0; kk<*cliqueSize; kk++) {
			cliqueRanks[kk] = cliqueInfo[2+kk];
		    }
		}
		else {
		    MPI_Send( cliqueInfo, nInClique+2, MPI_INT, 
			      cliqueInfo[2+k], 1, comm );
		}
	    }
	}

	/* Free all allocated memory */
	for (i=0; i<hashArray[rank]; i++) {
	    free( names[i].name );
	}
	free( names );
    }
    else {
	/* Participate in the Exscan */
	int zero = 0, dummy;
	MPI_Exscan( &zero, &dummy, 1, MPI_INT, MPI_SUM, comm );
	if (verbose) printf( "[%d] Computed clique index = %d\n", rank, dummy);
    }
    MPI_Wait( &req, MPI_STATUS_IGNORE );
    
    /* Receive the clique info (if we need to) */
    if (hash != rank) {
	MPI_Recv( cliqueInfo, 2+maxVal, MPI_INT, hash, 1, comm,
		  MPI_STATUS_IGNORE );

	/* Unpack */
	*cliqueNum = cliqueInfo[0];
	*cliqueSize = cliqueInfo[1]; 
	for (i=0; i<*cliqueSize; i++) {
	    cliqueRanks[i] = cliqueInfo[2+i];
	}
    }

    ttime = MPI_Wtime() - ttime;
    return 0;
}


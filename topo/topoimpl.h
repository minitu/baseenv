#ifndef TOPOIMPL_H_INCLUDED
#define TOPOIMPL_H_INCLUDED

/* This structure holds the available functions to determine node and
   interconnect information.  By using an array of these structures,
   we can provide multiple methods for determining node and interconnect
   information.  This is mostly useful for debugging, but can also be
   helpful in the case where only partial information is available */
typedef struct topoiNodeFns {
    int (*getNodeInfo)(int, topoinfo_t *);
    const char *descString;
    topoMethod_t id;
} topoiNodeFns_t;


typedef struct topoiNetFns {
    int (*getNodeTopoInfo)(topoinfo_t *);
    const char *descString;
    topoMethod_t id;
} topoiNetFns_t;

/*int topoiSetupDummyTopo(int kind, topoinfo_t *topo);*/
int topoiSetDummyTopo(int kind);
topoentry_t *topoiAllocEntry(topoinfo_t *);
int topoiGetNodeTopoInfoDEBUG(topoinfo_t *topo);
int topoiGetNodeInfoDEBUG(int isMultithreaded, topoinfo_t *topo);

#endif

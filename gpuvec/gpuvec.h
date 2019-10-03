#ifndef HAVE_GPUVEC_H
#define HAVE_GPUVEC_H 1

#ifndef BASETYPE
#define BASETYPE float
#endif
/* Stringize basetype */
#define _STR1(x) #x
#define _STR2(x) _STR1(x)
#define BASETYPE_NAME _STR2(BASETYPE)

#ifndef LEN
#define LEN 32000
#endif
#ifndef LEN2
#define LEN2 256
#endif

#endif

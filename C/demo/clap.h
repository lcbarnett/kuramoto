#ifndef CLAP_H
#define CLAP_H

// CLAP: a simple command-line parser
//
// argc, argv are the standard parameters for your `main' routine.

#define CLAP_ARGTYPE_MAX 100
#define CLAP_ARGNAME_MAX 100

typedef unsigned char     uchar;
typedef unsigned int      uint;
typedef unsigned long     ulong;
typedef long double       ldouble;
typedef const char* const cstr;

int clap_find_arg(int argc, char* argv[], cstr argname, cstr argtype, void* parg);

#define CLAP_ARG(argname,argtype,argdef,argdesc) \
const argtype argname = (clap_find_arg(argc,argv,#argname,#argtype,(void*)&argname) ? argname : argdef); \
void* clap_v##argname = (void*)&argname; \
if      (strncmp(#argtype,"char",    CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24c %s\n",  #argname, *(char*)    clap_v##argname, argdesc); \
else if (strncmp(#argtype,"int",     CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24d %s\n",  #argname, *(int*)     clap_v##argname, argdesc); \
else if (strncmp(#argtype,"uint",    CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24u %s\n",  #argname, *(uint*)    clap_v##argname, argdesc); \
else if (strncmp(#argtype,"long",    CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24ld %s\n", #argname, *(long*)    clap_v##argname, argdesc); \
else if (strncmp(#argtype,"ulong",   CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24lu %s\n", #argname, *(ulong*)   clap_v##argname, argdesc); \
else if (strncmp(#argtype,"size_t",  CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24zu %s\n", #argname, *(size_t*)  clap_v##argname, argdesc); \
else if (strncmp(#argtype,"float",   CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24g %s\n",  #argname, *(float*)   clap_v##argname, argdesc); \
else if (strncmp(#argtype,"double",  CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24lg %s\n", #argname, *(double*)  clap_v##argname, argdesc); \
else if (strncmp(#argtype,"ldouble", CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24Lg %s\n", #argname, *(ldouble*) clap_v##argname, argdesc); \
else if (strncmp(#argtype,"cstr",    CLAP_ARGTYPE_MAX) == 0) printf("%-12s = %-24s %s\n",  #argname, *(cstr*)    clap_v##argname, argdesc); \
else {fprintf(stderr,"CLAP ERROR: unhandled arg type \"%s\"\n",#argtype); exit(EXIT_FAILURE);} \
fflush(stdout)

#endif // CLAP_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "clap.h"

int clap_find_arg(int argc, char* argv[], cstr argname, cstr argtype, void* parg)
{
	--argc;
	++argv;
// fprintf(stderr,"argc = %d\n",argc);
	if (argc == 0) return 0;
	if (argc%2) {fprintf(stderr,"CLAP ERROR: must be an even number of args\n"); exit(EXIT_FAILURE);}
	for (int i = 0; i < argc; i += 2)
	{
// fprintf(stderr,"argv[%d] = %s\n",i,argv[i]);
		if (argv[i][0] != '-')  {fprintf(stderr,"CLAP ERROR: argument %d is not a switch!\n",i+1); exit(EXIT_FAILURE);}
		if (strncmp(argname,argv[i]+1, CLAP_ARGNAME_MAX) == 0)
		{
			if      (strncmp(argtype,"char",    CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%c",  (char*)   parg); return 1;}
			else if (strncmp(argtype,"int",     CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%d",  (int*)    parg); return 1;}
			else if (strncmp(argtype,"uint",    CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%u",  (uint*)   parg); return 1;}
			else if (strncmp(argtype,"long",    CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%ld", (long*)   parg); return 1;}
			else if (strncmp(argtype,"ulong",   CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%lu", (ulong*)  parg); return 1;}
			else if (strncmp(argtype,"size_t",  CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%zu", (size_t*) parg); return 1;}
			else if (strncmp(argtype,"float",   CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%g",  (float*)  parg); return 1;}
			else if (strncmp(argtype,"double",  CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%lg", (double*) parg); return 1;}
			else if (strncmp(argtype,"ldouble", CLAP_ARGTYPE_MAX) == 0) {sscanf(argv[i+1],"%Lg", (ldouble*)parg); return 1;}
			// this is so sneaky... we point *parg at the actual argv
			// string, so we don't need to allocate memory for it :)
			// [but we do break constness of cstr type (you gotta love C :) ... could this have repercussions?]
			else if (strncmp(argtype,"cstr",    CLAP_ARGTYPE_MAX) == 0) {*(char**)parg = argv[i+1];	return 1;}
		}

	}
	return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int demo      (int argc, char* argv[]);
int audio     (int argc, char* argv[]);
int stulan    (int argc, char* argv[]);
int chaos     (int argc, char* argv[]);
int lorenz96  (int argc, char* argv[]);
//int scratch (int argc, char* argv[]);

int main(int argc, char* argv[])
{
	if (!(argc>1)) {
		fprintf(stderr,"%s: must be at least one program argument\n",argv[0]);
		return EXIT_FAILURE;
	}

	if (strcmp(argv[1], "demo"     ) == 0) return demo     (argc-2,argv+2);
	if (strcmp(argv[1], "audio"    ) == 0) return audio    (argc-2,argv+2);
	if (strcmp(argv[1], "stulan"   ) == 0) return stulan   (argc-2,argv+2);
	if (strcmp(argv[1], "chaos"    ) == 0) return chaos    (argc-2,argv+2);
	if (strcmp(argv[1], "lorenz96" ) == 0) return lorenz96 (argc-2,argv+2);
//	if (strcmp(argv[1], "scratch"  ) == 0) return scratch  (argc-2,argv+2);

	fprintf(stderr,"%s: unknown simulation '%s'\n",argv[0],argv[1]);
	return EXIT_FAILURE;
}

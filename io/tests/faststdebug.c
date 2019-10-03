#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    char *filename = 0;
    int i, nline = 0;
    FILE *fp;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s filename linecount\n", argv[0]);
        return 1;
    }

    filename = argv[1];
    nline    = atoi(argv[2]);

    fp = fopen(filename, "w");
    for (i=0; i<nline; i++) {
	fprintf(fp, "%d: This is a line of text for line %d\n", i, i);
    }

    fclose(fp);
    return 0;
}

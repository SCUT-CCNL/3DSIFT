#include "../../Include/Util/matrixIO3D.h"

#include <cstring>
#include <cctype>

// writes out the m x n matrix header onto the stream
int WriteMatrixHeaderToStream(FILE * file, int m, int n, int p)
{
	if (fwrite(&m, sizeof(int), 1, file) < 1)
		return 1;
	if (fwrite(&n, sizeof(int), 1, file) < 1)
		return 1;
	if (fwrite(&p, sizeof(int), 1, file) < 1)
		return 1;

	return 0;
}

int ReadMatrixSizeFromStream(FILE * file, int * m, int * n, int *p)
{
	if (fread(m, sizeof(int), 1, file) < 1)
		return 1;
	if (fread(n, sizeof(int), 1, file) < 1)
		return 1;
	if (fread(p, sizeof(int), 1, file) < 1)
		return 1;

	return 0;
}


#ifndef _MATRIXIO3D_H_
#define _MATRIXIO3D_H_

#include <cstdio>
#include <cstdlib>

// === 3D Matrix Input and Output Rountines ===
/*
* Read and write image from and to disk binary files. *
* File Format: m n p, data...                         *
* NOTE: Matrix are saved as row-major                 *
*/


int ReadMatrixSizeFromStream(FILE * file, int * m, int * n, int *p);

int ReadMatrixSizeFromDisk(const char * filename, int * m, int * n, int *p);

int WriteMatrixHeaderToStream(FILE * file, int m, int n, int p);

template <class T>
int ReadMatrixFromDisk(const char * filename, int * m, int * n, int *p, T ** matrix)
{
	FILE * file;
	file = fopen(filename, "rb");
	if (!file)
	{
		printf("Can't open input matrix file: %s.\n", filename);
		return 1;
	}

	if (ReadMatrixSizeFromStream(file, m, n, p) != 0)
	{
		printf("Error reading matrix header from disk file: %s.\n", filename);
		return 1;
	}

	//int size = (*m) * (*n) * sizeof(T) + 2 * sizeof(int);
	*matrix = (T *)malloc(sizeof(T)*(*m)*(*n)*(*p));

	if (ReadMatrixFromStream(file, *m, *n, *p, *matrix) != 0)
	{
		printf("Error reading matrix data from disk file: %s.\n", filename);
		return 1;
	}

	fclose(file);

	return 0;
}

// read the m x n x p matrix from the stream, in binary format
template <class T>
int ReadMatrixFromStream(FILE * file, int M, int N, int P, T * matrix)
{
	unsigned int readBytes;
	if ((readBytes = fread(matrix, sizeof(T), M*N*P, file)) < (unsigned int)M*N*P)
	{
		printf("Error: I have only read %u bytes. sizeof(T)=%zu\n", readBytes, sizeof(T));
		return 1;
	}

	return 0;
}

//Writes matrix to disk binary files
template <class T>
int WriteMatrixToDisk(const char * filename, int m, int n, int p, T * matrix)
{
	FILE * file;
	file = fopen(filename, "wb");
	if (!file)
	{
		printf("Can't open output file: %s.\n", filename);
		return 1;
	}

	if (WriteMatrixHeaderToStream(file, m, n, p) != 0)
	{
		printf("Error writing the matrix header to disk file: %s.\n", filename);
		return 1;
	}

	if (WriteMatrixToStream(file, m, n, p, matrix) != 0)
	{
		printf("Error writing the matrix to disk file: %s.\n", filename);
		return 1;
	}

	fclose(file);

	return 0;
}

// writes out the m x n x p matrix onto the stream, in binary format
// Row-major
template <class T>
int WriteMatrixToStream(FILE * file, size_t m, size_t n, size_t p, T * matrix)
{
	if (fwrite(matrix, sizeof(T), m*n*p, file) < m*n*p)
		return 1;
	return 0;
}

// prints the matrix to standard output in Matlab format
//Print Matrix to command line
template <class T>
void PrintMatrixInMatlabFormat(int m, int n, int p, T * U)
{
	for (int k = 0; k < p; k++) {
		for (int i = 0; i<n; i++)
		{
			printf("{");
			for (int j = 0; j < m; j++)
			{
				printf("%f", U[ELT(n, p, j, i, k)]);
				if (j != m - 1)
					printf(", ");
			}
			printf("}");

			if (i != n - 1)
				printf(",\n");
		}

		printf("\n}\n");
	}
}
#endif // !_MATRIXIO3D_H_

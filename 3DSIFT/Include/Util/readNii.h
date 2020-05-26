#ifndef __READ_NII_H__
#define __READ_NII_H__

#include "common.h"

SIFT_LIBRARY_API float* readNiiFile(const char *fname, int &nx, int &ny, int &nz);

#endif
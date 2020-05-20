#ifndef __SIFT_TIME_H__
#define __SIFT_TIME_H__

#include "common.h"

struct SIFT_LIBRARY_API SIFT_TimerPara
{
	double d_Allocation = 0;
	double d_BuildGSS = 0;
	double d_BuildDOG = 0;
	double d_Detect = 0;
	double d_AssignOrientation = 0;
	double d_Extraction = 0;
	double d_release = 0;
};

struct SIFT_LIBRARY_API SIFT_PROCESS
{
	SIFT_TimerPara REF;
	SIFT_TimerPara TAR;
	double d_RegTime = 0;
};

#endif

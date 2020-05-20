#ifndef  __COMMON_H__
#define __COMMON_H__

#ifdef SIFT_LIBRARY_EXPORTS
#define SIFT_LIBRARY_API __declspec(dllexport)
#else
#define SIFT_LIBRARY_API __declspec(dllimport)
#endif

//#define CHECK_ENABLE
//#define JR_DEBUG
#include <vector>
#include <iostream>
#include <chrono>

template <class T>
T getDoubleMill(decltype(std::chrono::high_resolution_clock::now()) start, decltype(std::chrono::high_resolution_clock::now()) end) {
	std::chrono::duration<T, std::milli> totalMs = end - start;
	return totalMs.count();
}

struct SIFT_LIBRARY_API SIFT_TimerPara {
	double d_TotalTime = 0;

	double d_Allocation = 0;
	double d_BuildGSS = 0;
	double d_BuildDOG = 0;
	double d_Detect = 0;
	double d_AssignOrientation = 0;
	double d_Extraction = 0;
	double d_release = 0;

	double d_memoryOverhead = 0;

	std::vector<double> vD_octaveTime;
	std::vector<double> vD_octaveCompute;

	double getAllComputeTime() {
		return (d_BuildGSS + d_BuildDOG + d_Detect + d_AssignOrientation + d_Extraction);
	}
};

struct SIFT_LIBRARY_API SIFT_Dev_TimerPara {

	int devId = 0;

	double d_Allocation = 0;
	double d_hostCopy = 0;
	double d_BuildGSS = 0;
	double d_Detect = 0;
	double d_KeypointCompute = 0;
	double d_downsample = 0;
};

struct SIFT_LIBRARY_API SIFT_PROCESS{
	SIFT_TimerPara REF;
	SIFT_TimerPara TAR;
	double d_RegTime = 0;
};

SIFT_LIBRARY_API std::ostream & operator << (std::ostream &os, const SIFT_TimerPara &st);

SIFT_LIBRARY_API std::ostream & operator << (std::ostream &os, const SIFT_PROCESS &sp);


#endif // ! __COMMON_H__

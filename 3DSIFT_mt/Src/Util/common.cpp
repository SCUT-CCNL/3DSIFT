#include "../../Include/Util/common.h"

using namespace std;

std::ostream & operator << (std::ostream &os, const SIFT_TimerPara &st) {
	os << "\tTotal Time:" << st.d_TotalTime << "\n"
		<< "\tAllocation time:" << st.d_Allocation << "\n"
		<< "\tBuild-GSS time:" << st.d_BuildGSS << "\n"
		<< "\tBuild-DOG time:" << st.d_BuildDOG << "\n"
		<< "\tDetection time:" << st.d_Detect << "\n"
		<< "\tOrientation time:" << st.d_AssignOrientation << "\n"
		<< "\tExtract time:" << st.d_Extraction << "\n"
		<< "\tRelease time:" << st.d_release << "\n"
		<< "\tMemory ovhead time:" << st.d_memoryOverhead << std::endl;

	if (st.vD_octaveTime.size() == st.vD_octaveCompute.size() && !st.vD_octaveTime.empty()) {
		for (int i = 0; i < static_cast<int>(st.vD_octaveTime.size()); ++i) {
			os << "\t----Octave " << i << "\ttime:" << st.vD_octaveTime[i] << "\tcompute-time:" << st.vD_octaveCompute[i] << std::endl;
		}
	}

	return os;
}

std::ostream & operator << (std::ostream &os, const SIFT_PROCESS &sp) {
	os << "Reference:\n"
		<< sp.REF << endl;

	os << "Target:\n"
		<< sp.TAR << endl;

	os << "Match:\n"
		<< sp.d_RegTime;

	return os;
}
#include "Include/cSIFT3D.h"
#include "Include/Util/matrixIO3D.h"

using namespace std;

int main(){

	const char *Ref = R"(D:\SIFT_EXPERI_DATA\Hole_intersect\624_500-Ori440Margin92_92_30-m12_std8_r8-intersect_6000-GN.bin)";
	//const char *path_ = R"(F:\SIFT_EXPERI_LOG\DEBUG\Voxel_384_ori_GN_.txt)";

	int m, n, p;
	float* SIFT_Vol = nullptr;
	ReadMatrixFromDisk(Ref, &m, &n, &p, &SIFT_Vol);
	cout << m << " " << n << " " << p << endl;

	auto CSIFT3D = CPUSIFT::CSIFT3DFactory::CreateCSIFT3D(SIFT_Vol, m, n, p);

	CSIFT3D->KpSiftAlgorithm();
	auto vKP = CSIFT3D->GetKeypoints();


	return 0;
}
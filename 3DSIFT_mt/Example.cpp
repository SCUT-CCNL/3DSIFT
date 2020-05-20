#include "Include/cSIFT3D.h"
#include "Include/Util/matrixIO3D.h"
#include "Include/Util/readNii.h"

using namespace std;

int main(){

    const char *Ref = R"(D:\SIFT_EXPERI_DATA\Hole_intersect\624_500-Ori440Margin92_92_30-m12_std8_r8-intersect_6000-GN.bin)";
    const char *RefNii = R"(D:\SIFT_EXPERI_DATA\Hole_intersect\624_500-Ori440Margin92_92_30-m12_std8_r8-intersect_6000-GN.bin)";
	//const char *path_ = R"(F:\SIFT_EXPERI_LOG\DEBUG\Voxel_384_ori_GN_.txt)";

	int m, n, p;
	float* SIFT_Vol = nullptr;
	ReadMatrixFromDisk(Ref, &m, &n, &p, &SIFT_Vol);
    cout << m << " " << n << " " << p << endl;

    int m1, n1, p1;
    float* SIFT_VolNii = nullptr;
    SIFT_VolNii = readNiiFile(RefNii, m1, n1, p1);
    cout << m1 << " " << n1 << " " << p1 << endl;

    float bias = 0;
    int64_t count = 0;
    const int64_t nx = m;
    const int64_t nxy = m*n;
    for (int64_t x = 0; x < m; ++x) {
        for (int64_t y = 0; y < n; ++y) {
            for (int64_t z = 0; z < p; ++z) {
                int64_t idx = z* nxy + y*nx + x;
                float vBin = SIFT_Vol[idx];
                float vNii = SIFT_VolNii[idx];

                float tmp = abs(vNii - vBin);
                if (tmp > 0.01f) {
                    count += 1;
                    std::cerr << "warning, high bias at (" << x << "," << y << "," << z << ") : " << tmp << std::endl;
                }
                bias += tmp;
            }
        }
    }

    std::cout << "total bias:" << bias<<std::endl;
    std::cout << "bias num:" << count << std::endl;

	auto CSIFT3D = CPUSIFT::CSIFT3DFactory::CreateCSIFT3D(SIFT_Vol, m, n, p);

	CSIFT3D->KpSiftAlgorithm();
	auto vKP = CSIFT3D->GetKeypoints();


	return 0;
}
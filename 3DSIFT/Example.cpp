#include "Include/cSIFT3D.h"
#include "Include/Util/readNii.h"
#include "Include/cMatcher.h"
#include <vector>

using namespace std;

int main() {

    const char *RefNiiPath = R"(../Example/Torus_Ref.nii.gz)";
    const char *TarNiiPath = R"(../Example/Torus_Def.nii.gz)";

    
    //read reference image
    int nx = 0, ny = 0, nz = 0;
    float* refVolNii = nullptr;
    refVolNii = readNiiFile(RefNiiPath, nx, ny, nz);
    cout << "Dimensions of reference image:" << nx << " " << ny << " " << nz << endl;

    //extract keypoint from ref image
    auto SIFT_ref = CPUSIFT::CSIFT3DFactory::CreateCSIFT3D(refVolNii, nx, ny, nz);
    SIFT_ref->KpSiftAlgorithm();
    auto vRefKp = SIFT_ref->GetKeypoints();
    

    
    //read target image
    int nxTar = 0, nyTar = 0, nzTar = 0;
    float* tarVolNii = nullptr;
    tarVolNii = readNiiFile(TarNiiPath, nxTar, nyTar, nzTar);
    cout << "Dimensions of target image:" << nxTar << " " << nyTar << " " << nzTar << endl;
    //extract keypoint from tar image
    auto SIFT_tar = CPUSIFT::CSIFT3DFactory::CreateCSIFT3D(tarVolNii, nxTar, nyTar, nzTar);
    SIFT_tar->KpSiftAlgorithm();
    auto vTarKp = SIFT_tar->GetKeypoints();
    

    //match procedure
    CPUSIFT::muBruteMatcher matcher;
    //the coordinates of matched keypoint pairs in reference image and target image.
    vector<CPUSIFT::Cvec> matchRefCoor, matchTarCoor;
    //threshold of filtering matches
    const float threshold = 0.85f;
    matcher.enhancedMatch(matchRefCoor, matchTarCoor, vRefKp, vTarKp, threshold);

    cout << "Matched Points: reference coordinate(x,y,z);target coordinate(x,y,z)" << endl;
    for (int i = 0; i < matchRefCoor.size(); ++i) {
        cout <<
            matchRefCoor[i].x << "," <<
            matchRefCoor[i].y << "," <<
            matchRefCoor[i].z << ";" <<

            matchTarCoor[i].x << "," <<
            matchTarCoor[i].y << "," <<
            matchTarCoor[i].z << endl;
    }

    delete[] refVolNii;
    delete[] tarVolNii;
    delete SIFT_ref;
    delete SIFT_tar;

	return 0;
}
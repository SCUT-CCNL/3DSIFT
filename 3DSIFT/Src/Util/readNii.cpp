#include "../../Include/Util/readNii.h"

#include "../../3party/layNii/dep/laynii_lib.h"

float* readNiiFile(const char *fname, int &nx, int &ny, int &nz) {

    nifti_image* nii = nifti_image_read(fname, 1);

    nx = nii->nx;
    ny = nii->ny;
    nz = nii->nz;

    const int64_t nxy = int64_t(nx)*int64_t(ny);
    const int64_t nxyz = int64_t(nx)*int64_t(ny)*int64_t(nz);
    
    //if the input file is not FLOAT32, then convert it 
    if (nii->datatype != DT_FLOAT32) {
        auto niiFloat = copy_nifti_as_float32(nii);
        nifti_image_free(nii);
        nii = niiFloat;
    }
    float *niiData = static_cast<float*>(nii->data);

    //malloc data
    float *data = new float[nxyz];

    //fetch voxels
    for (int64_t z = 0; z < nz; ++z) {
        for (int64_t y = 0; y < ny ; ++y) {
            for (int64_t x = 0; x < nx; ++x) {
                int64_t voxel_idx = nxy * z + nx * y + x;
                data[voxel_idx] = *(niiData + voxel_idx);
            }
        }
    }

    nifti_image_free(nii);
    return data;
}
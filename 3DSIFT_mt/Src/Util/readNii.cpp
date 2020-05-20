#include "../../Include/Util/readNii.h"

#include <laynii_lib.h>

float* readNiiFile(const char *fname, int &nx, int &ny, int &nz) {

    nifti_image* nii = nifti_image_read(fname, 1);

    nx = nii->nx;
    ny = nii->ny;
    nz = nii->nz;

    const int64_t nxy = int64_t(nx)*int64_t(ny);
    const int64_t nxyz = int64_t(nx)*int64_t(ny)*int64_t(nz);
    
    auto niiFloat = copy_nifti_as_float32(nii);
    float* niiData = static_cast<float*>(niiFloat->data);

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
    nifti_image_free(niiFloat);

    return data;
}
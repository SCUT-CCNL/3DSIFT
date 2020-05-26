/*
BSD 3-Clause License

Copyright (c) 2020, Laurentius Huber
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//
// some small modifications:
// 1.remove "using namespace std;" in the header filebuf
// 2.remove "#include <math.h>; " and "#include <cmath>", and add the latter one to the source file
// 3.add include guard

#ifndef __LAYNII_LIB_H__
#define __LAYNII_LIB_H__

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <tuple>

#include "./nifti2_io.h"

//using namespace std;

// ============================================================================
// Declarations
// ============================================================================
double ren_average(double arr[], int size);
double ren_stdev(double arr[], int size);
double ren_correl(double arr1[], double arr2[], int size);
double ren_skew(double arr[], int size);
double ren_kurt(double arr[], int size);
double ren_autocor(double arr[], int size);

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ);
float dist2d(float x1, float y1, float x2, float y2);
float angle(float a, float b, float c);
float gaus(float distance, float sigma);

void log_welcome(const char* programname);
void log_output(const char* filename);
void log_nifti_descriptives(nifti_image* nii);

void save_output_nifti(std::string filename, std::string prefix, nifti_image* nii,
                       bool log = true, bool use_outpath = false);

nifti_image* copy_nifti_as_double(nifti_image* nii);
nifti_image* copy_nifti_as_float32(nifti_image* nii);
nifti_image* copy_nifti_as_float16(nifti_image* nii);
nifti_image* copy_nifti_as_int32(nifti_image* nii);
nifti_image* copy_nifti_as_int16(nifti_image* nii);

std::tuple<uint32_t, uint32_t, uint32_t> ind2sub_3D(
    const uint32_t linear_index, const uint32_t size_x, const uint32_t size_y);

uint32_t sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                    const uint32_t size_x, const uint32_t size_y);

std::tuple<float, float> simplex_closure_2D(float x, float y);
std::tuple<float, float> simplex_perturb_2D(float x, float y, float a, float b);

// ============================================================================
// Preprocessor macros.
// ============================================================================
#define PI 3.14159265;

#endif //__LAYNII_LIB_H__
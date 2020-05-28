# 3DSIFT
A multi-thread CPU implementation of 3D SIFT, the Scale invariant feature transform (SIFT) for 3D image (volumetric image). The feature matching of 3D SIFT features is also provided.

The API of reading NIFTI images is included in this program.



The program is written in C++ language and is parallelized using OpenMp.



## Contents

- **Read nifti image:**
  - Header file:`Include/Util/readNii.h`
  - Function:`readNiiFile(path, nx, ny, nz);`
  - The function read the `path`, and return a pointer with the data of volume image, and the dimensions of the image are written on `nx`, `ny` and `nz` .
- **3D SIFT feature extraction.**
  - Header file:`Include/cSIFT3D.h`
  - Object create function:`CPUSIFT::CSIFT3DFactory::CreateCSIFT3D(volData, nx, ny, nz)`
    - This is a factory function belongs to namespace `CPUSIFT`, it creates a 3D SIFT object with the datapointer `volData` and dimensions `nx`,`ny` and `nz` as input. Also, the parameters of 3D SIFT can be indicated, they are set as default values, as seen in the declaration of the function.
    - The created 3D SIFT object should be deleted.
  - Extract class function:`    CSIFT3D->KpSiftAlgorithm();`
    - perform feature extraction procedure
  - Obtain keypoint function:`CSIFT3D->GetKeypoints()`
    - return vector of extracted keypoints.
- **Feature matching on two sets of features.**
  - Header file:`Include/cMatcher.h`
  - Class: `CPUSIFT::muBruteMatcher`
  - This class provides several functions of brute-force matching. 
    - **Overall**: The `enhancedMatch` is recommended.
    - **Parameters:** `(matchRefCoor, matchTarCoor, vRefKp, vTarKp, thrshold)`
      - The first two parameters are coordinates of returned matched keypoints, representing coordinates in reference image and target image.
      - `vRefKp` and `vTarKp` represents the two sets of keypoints extracted from the image
      - `threshold` is the threshold for filtering ambiguous matches. The smaller the value, the more strict the filtering. `0.85` is recommended and used as default parameter.
    - `injectMatch()`: perform feature matching from reference set of keypoint to target set of keypoint. Like mostly used in 2D SIFT.
    - `bijectMatch()`: perform matching from reference to target, and also from target to reference, the valid matches must be appeared in the both two matching, as used in Ref[1].
    - `enhancedMatch()`: perform match from the reference to target first, and then for target features that are matched by multiple reference feature the matching from target to reference is performed to remove ambiguity, as used in Ref[2].



## How to use

Currently, the visual studio project is provided, it's recommended to use windows PC to run the codes.

A example of feature extraction and matching is provided in `example.cpp`, which is included in the visual studio project too. The steps of running the example data is as follows:

1. Open the visual studio project `.sln` file.
2. Set the Windows SDK version of the visual studio project to that are installed on your PC. PS: version below  10.0.17763.0 is non-tested.
3. Modify the path of data files in `example.cpp`
4. Build and run the codes



If the readers want to use this library in another programs, two methods are recommended:

1. Put all of the source codes of this project (besides `example.cpp`) into the another program, and the dependencies of different files and third-party libraries should be handled correctly. An example can be seen in our another project [3D SIFT PiDVC]( https://github.com/ParallelCCNL/3DSIFT_PiDVC_mt ).
2. Build `dll` and `lib` files using the visual studio projects by setting the `Configuration Types` of visual studio project to be `Dynamic Library (.dll)` instead of `Application (.exe)`. And in the program importing 3D SIFT, it's required to include only the header files of 3D SIFT and import the `.lib` flies for compilation, and the `.dll` files are required for running the generated programs.



## Example data

Several nifti files are provided:

https://drive.google.com/open?id=1IlOBp1hu-KUh648lXR3YjgNih71Mj8ke



## Dependency

Several third-party libraries are used in this work, and are also included in this project.

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), put in path `3DSIFT\3party\Eigen`, used to perform matrix  computation, such as eigendecompostion.
- [layNii](https://github.com/layerfMRI/LAYNII), put in path `3DSIFT\3party\layNii`, used to read Nifti format files.
- [Zlib](https://zlib.net/), put in path `3DSIFT\3party\zlib`, used by layNii. It should be noticed that the provided lib files are built on x64 Windows 10 with Visual studio 2017 and WINDOWS SDK 10.0.17763.0 from the source code. 



## About 3D SIFT

This program is based on the work of [Rister](https://github.com/bbrister/SIFT3D). He has proposed a well-performed 3D SIFT with rotation-invariance, as seen in Reference [1].

There are several different details in our  implement:

- The Gaussian Smoothing is fixed in this work, by removing the zooming of gaussian kernel in higher octave. A series of gaussian kernels are used in different octaves without zooming of kernel radius. The produced  gaussian images are verified by comparing with that produced by matlab with the same sigma parameter.
- In descriptor construction, the keypoint is located at the center of the sub-regions as described in the paper (reference [1]) , while in this [implementation](https://github.com/bbrister/SIFT3D) the keypoint is coincide with one of the sub-region center.



## Reference and Citations

If you want to cite this work, please refer to the papers as follows:

[1] B. Rister, M. A. Horowitz and D. L. Rubin, "Volumetric Image Registration From Invariant Keypoints," in *IEEE Transactions on Image Processing*, vol. 26, no. 10, pp. 4900-4910, Oct. 2017. doi: 10.1109/TIP.2017.2722689

[2] Our paper in the manuscript.
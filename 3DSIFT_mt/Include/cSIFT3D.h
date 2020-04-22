#ifndef __CSIFT3D_H__
#define __CSIFT3D_H__

#include <string>
#include <vector>
#include "Util/cTexImage.h"
#include "Util/common.h"

namespace CPUSIFT {
	//define some constant
	//#define INIT_SIGMA 0.5 // Use for 2D blurring first time

#define SIGMA_DEFAULT 1.6 //Scale of the base octave
#define SIGMA_N_DEFAULT 1.15 //Nomial scale,即视为正常尺度
#define NUM_KP_LEVELS 3 //Number of levels per octave in which keypoints are found


#define PEAK_THRESH 0.1 //DoG peak threshold
#define EIG_THRES 0.9
#define CORNER_THRESH 0.4 // Minimum corner score

//internel params
#define IMG_BORDER 1 //Consideration
#define NHIST_PER_DIM 4 // 4x4x4 for the cube
#define ICOS_NFACES 20 // Number of faces in an icosahedron
#define ICOS_NVERT 12 // Number of vertices in an icosahedron
#define DESC_NUMEL (NHIST_PER_DIM * NHIST_PER_DIM * NHIST_PER_DIM * ICOS_NVERT)

/*
!------------------------------------------------------------------
Following struct needs to be reorganized the element order for efficiency
!------------------------------------------------------------------
*/
	//global for
	extern int sift_thread_num;

	//Struct defining a vector in Cartesian coordinates
	typedef struct _cCvec
	{
		float x;
		float y;
		float z;
		_cCvec(float x_ = 0, float y_ = 0, float z_ = 0) {
			x = x_;
			y = y_;
			z = z_;
		}

	} SIFT_LIBRARY_API Cvec;

	//Primary structure, keypoint
	typedef struct _cKeypoint
	{
		float x, y, z;// location coordinate, float maybe better here for sub-pixel
		float scale;//scale
		int octave, level; //pyramid index

		float rx, ry, rz;//Original Image location

		//for debug
		Cvec win;
		float eigvalue[3];
		float eigvector[9];

		//Descriptor and Rotation Matrix, structure tensor should be added following
		float Rotation[9]; //3x3 matrix for rotation
		float str_tensor[9]; //3X3 matrix for structure_tensor
		float *desc = nullptr; //Corresponding to the 4x4x4x12 bins

	} SIFT_LIBRARY_API Keypoint;

	bool cmp_kp(const Keypoint &a, const Keypoint &b);

	bool cmp_kp_orig(const Keypoint &a, const Keypoint &b);

	//Triangle definition
	typedef struct _cTri
	{
		Cvec v[3]; // Vertices
		int idx[3]; // Index of each vertex in the solid

	}Tri;

	// Triangle mesh
	typedef struct _Mesh
	{
		Tri *tri; //Triangles
		int num; //Number of triangles

	}Mesh;

	//struct for 3D image, self-defined, channels default to be 1
	//Base element for Gaussian and DOG pyramid
	typedef struct _cImage
	{
		float* data;// from raw data, don't use unsigned char here
		//float* data_f; // convert to float data, maybe better

		int nx, ny, nz;// Image dimension in x-axis, y-axis, z-axis
		size_t xs, ys, zs;//Stride in x,y,z: xs = 1, ys = nx, zs = nx * ny

		//Not very clear about the float world dimensions
		float ux, uy, uz;

		//size_t -> unsigned long long
		size_t size; // Total size in pixels
		float s; // scale-space location

	}Image;

	typedef struct _cEigenVal
	{
		float val;// eigen value
		float vec[3]; // eigen vector

	}EigenVal;

	class CSIFT3D
	{
	protected:
		std::vector<Keypoint> filter;
		std::vector<TexImage> Gss_Pyramid;
		std::vector<TexImage> DoG_Pyramid;
		TexImage Host_Im;
		Mesh mesh;
		int octave_num;
		int num_kp_levels;
		float sigma_default;     // Scale of the base octave
		float sigma_n_default;   // Nomial scale
		float peak_thresh;       // DoG peak threshold
		float max_eig_thres;	 //Maximum ratio of eigenvalue magnitudes
		float corner_thresh;     // Minimum corner score
		std::vector<Keypoint> extre;

		//for debug
		std::vector<std::vector<Keypoint>> level_extrema;

		//global descriptor
		// DESC_NUMEL * N
		float *global_descriptor = nullptr;

	public:
		SIFT_TimerPara m_timer;

		SIFT_LIBRARY_API CSIFT3D();
		SIFT_LIBRARY_API CSIFT3D(float* volume, int x_dim, int y_dim, int z_dim,
			int num_kp_levels_, float sigma_default_,
			float sigma_n_default_, float peak_thresh_, float max_eig_thres_, float corner_thresh_);
		SIFT_LIBRARY_API ~CSIFT3D();

		SIFT_LIBRARY_API void KpSiftAlgorithm();

		//for running parameter settings
		SIFT_LIBRARY_API void SetNumThreads(int t_num);
		SIFT_LIBRARY_API std::vector<Keypoint> GetKeypoints();

		//!- Ordinary Algorithm process
		void Initialize();
		void Build_Gaussian_Scale_Space();
		void Build_DOG_Scale_Space();
		void Detect_KeyPoints();
		void Assign_Orientation();
		void Extract_Description();
		void Release_SIFT();
		void SetHostImNull() { Host_Im.SetImageDataPt(nullptr); }

		//added by JRYoung for checking
		//MAY 2019
		std::vector<TexImage>* GET_GSS() {
			return &Gss_Pyramid;
		};
		std::vector<TexImage>* GET_DOG() {
			return &DoG_Pyramid;
		};
		std::vector<std::vector<Keypoint> >* GET_LEVEL() {
			return &level_extrema;
		}

		//test funs
		void test_build2sigma(TexImage &img);

	};

	class SIFT_LIBRARY_API CSIFT3DFactory
	{
	public:
		static CSIFT3D* CreateCSIFT3D(float* volume,
			int x_dim, int y_dim, int z_dim,
			int num_kp_levels = NUM_KP_LEVELS,
			float sigma_default = SIGMA_DEFAULT,
			float sigma_n_default = SIGMA_N_DEFAULT,
			float peak_thresh = PEAK_THRESH,
			float max_eigo_thres = EIG_THRES,
			float corner_thresh = CORNER_THRESH);

		static CSIFT3D* CreateCSIFT3D(std::string path_,
			int num_kp_levels = NUM_KP_LEVELS,
			float sigma_default = SIGMA_DEFAULT,
			float sigma_n_default = SIGMA_N_DEFAULT,
			float peak_thresh = PEAK_THRESH,
			float max_eigo_thres = EIG_THRES,
			float corner_thresh = CORNER_THRESH);
	
	};

	// Interface for implementing the algorithm, refactor code

	void DownSample_3D(TexImage * src, TexImage * dst);

	void GaussianSmooth_3D(TexImage * src, TexImage * dst, float sigma);

	void GaussianSmooth_3D_Imp(TexImage * src, TexImage * dst, int dim, float unit, float* weight, int width);

	void Im_permute(TexImage * src, TexImage * dst, int dim1, int dim2);

	void Sub(TexImage * prev, TexImage * cur, TexImage * dog);

	bool IsExtrema_neighbor(TexImage* prev, TexImage* cur, TexImage* next, int x, int y, int z);

	int Assign_Orientation_Imp(Keypoint& kp, TexImage* gaussian, const float sigma, const float max_eig_ratio, const float corner_thresh);

	bool DistinctEig(float a, float b, float c);

	void Extract_Descriptor_Imp(Keypoint& kp, TexImage * gaussian, Mesh* mesh);

	void Trilinear_interpolation_over_desc(Mesh * mesh, Keypoint &kp, Cvec &vbins, Cvec &grad, int loop_idx);

	int Check_intersect_faces(Mesh *mesh, Cvec* grad, Cvec* bary);

	void Transpose_Matrix(float* Rot);

	void Swap_Element(float &a, float &b);

	void normailize_desc(float* desc);

	int cart2bary(Cvec * cart, const Tri * const tri, Cvec * const bary, float * const k);

	void Trilinear_interpolation_over_desc_debug(Mesh * mesh, Keypoint &kp, Cvec &vbins, Cvec &grad, int loop_idx,
		float* host_dvbins, int* host_intersect_id, float* host_bary, int* host_offset, float* host_desc_accum, int debug);

}

#endif // !__SIFT3D_H__

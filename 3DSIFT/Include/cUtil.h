#ifndef __CUTIL_H__
#define __CUTIL_H__

#include <vector>
#include "cSIFT3D.h"
#include "Util/common.h"

namespace CPUSIFT {

//Need to modify to be the stride, better to understand
#ifndef ELT
#define ELT(height,width,x,y,z) ((((long)z)*((long)height)+((long)y))*((long)width)+(long)x)

#define ELT2(x, y, z, x_stride, y_stride, z_stride) (long)((x) * x_stride + (y) * y_stride + (z) * z_stride)

//#define CMP_PREV(im, x, y, z, CMP, val) ( \
//        (val) CMP im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)] \
//)
//#define CMP_CUR(im, x, y, z, CMP, val) ( \
//	(val) CMP im->data[ELT2(x + 1, y, z, im->xs, im->ys, im->zs)] && \
//	(val) CMP im->data[ELT2(x - 1, y, z, im->xs, im->ys, im->zs)] && \
//	(val) CMP im->data[ELT2(x, y + 1, z, im->xs, im->ys, im->zs)] && \
//	(val) CMP im->data[ELT2(x, y - 1, z, im->xs, im->ys, im->zs)] && \
//	(val) CMP im->data[ELT2(x, y, z + 1, im->xs, im->ys, im->zs)] && \
//	(val) CMP im->data[ELT2(x, y, z - 1, im->xs, im->ys, im->zs)] \
//)
//#define CMP_NEXT(im, x, y, z, CMP, val) \
//        CMP_PREV(im, x, y, z, CMP, val)
#endif

//int Initialize(SIFT3D *sift3D);

//!-------------------------- refactor function in new version -----------------------------------------
	float SIFT3D_CVEC_L2_NORM(Cvec vec);

	float SIFT3D_CVEC_L2_NORM_SQ(Cvec vec);

	void SIFT3D_CVEC_SCALE(Cvec& vec, double sca);

	float SIFT3D_CVEC_DOT(Cvec& vec1, Cvec& vec2);

	void SIFT3D_CVEC_CROSS(Cvec& vec1, Cvec& vec2, Cvec& out);

	void Im2TexIm(Image* im, TexImage* tex_im);

	int Initialize_geometry(Mesh * mesh);

	void Initialize_Pyramid(std::vector<TexImage> &Pyramid, TexImage* Base_Im,
		int num_oct, int interval, int num_kp_levels, float sigma_default);

	// Refactor debug function
	void new_write_Im(TexImage* im);

	void new_write_Im(TexImage* im, int idx);

	void new_write_GSS(std::vector<TexImage> &Gss_Pyramid, int octave_num, int gss_level);

	void new_write_DOG(std::vector<TexImage> &DoG_Pyramid, int octave_num, int dog_level);

	void new_write_Desc(std::vector<Keypoint> &filter);

	void new_write_Keypoint(std::vector<Keypoint> &kp);

	void new_write_level_Keypoint(std::vector<std::vector<Keypoint>>& level_extrema);

	SIFT_LIBRARY_API void write_sift_kp(std::vector<Cvec> &kp, const char* file_name);

	SIFT_LIBRARY_API void read_sift_kp(const char* file_name, std::vector<Cvec> &kp);

	void write_debug_middle(float* loop_point, float* vrot, float* vbins, int* intersect_id, float* bary,
		float* dvbins, int* offset, float* desc_accum, int batchsize, Keypoint& kp, float* Rot);

	//!-----------------------------------------------------------------------------------------------------
	//int Initialize_geometry(SIFT3D *sift3D);

	//int Initialize_Pyramid(SIFT3D * sift3D);

	//int Initialize_Parameter(SIFT3D *sift3D);

	void Initialize_Keypoint(Keypoint& kp);

	int get_Min(int row, int col, int height);

	void CopyVolumeImage(Image * src, Image * dst);//check whether to use & or *

	void init_im(Image *const im);

	void im_free(Image * im);

	void im_scale(Image * im);

	void data_scale(float * data, int nx, int ny, int nz, size_t xs, size_t ys, size_t zs);

	float im_max_abs(Image * im);

	float im_max_abs(TexImage * im);

	//for debug and compare
	//void write_Im(SIFT3D * sift3D);
	//
	//void write_Im(Image* im, int idx);
	//
	//void write_GSS(SIFT3D * sift3D);
	//
	//void write_DOG(SIFT3D * sift3D);
	//
	//void write_Keypoint(SIFT3D * sift3D);
	//
	//void write_Desc(SIFT3D * sift3D);

	void write_Str_tensor(std::vector<Keypoint> extrema);

	void write_Eig(std::vector<Keypoint> extrema);

	void write_Rot(std::vector<Keypoint> filter);

	void write_Desc(std::vector<Keypoint> filter);

	SIFT_LIBRARY_API void read_desc_all_level(char* file_name, std::vector<Keypoint>& kp);

}

#endif // !__UTIL_H__

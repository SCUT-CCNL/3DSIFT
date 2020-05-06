
#include "../Include/cSIFT3D.h"

//#include "dicom.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>

#include "../Include/Util/matrixIO3D.h"
#include "../Include/cutil.h"
#include "../Include/Util/cMemManager.h"

using namespace std;
using namespace Eigen;

namespace CPUSIFT{

//#define __JR_DEBUG

const float ori_grad_thresh = 1E-10; //Minimum norm of average gradient
const float bary_eps = FLT_EPSILON * 1E1; //Error tolerance for barycentric coordinates
const float trunc_thresh = 0.2 * 128 / DESC_NUMEL; //Descriptor truncation threshold

//Some details need to understand
const float ori_sig_fctr = 1.5; //Ratio of window parameter to keypoint scale
const float ori_rad_fctr = 3.0; //Ratio of window radius to parameter

const float desc_sig_fctr = 7.071067812;//Ratio of window parameter to keypoint scale, 5 * sqrt(2)
const float desc_rad_fctr = 2.0; //Ratio of window radius to parameter

int sift_thread_num = omp_get_max_threads();

bool cmp(EigenVal a, EigenVal b)
{
	return a.val < b.val;
}

bool cmp_kp(const Keypoint &a, const Keypoint &b) {
	if (a.rz < b.rz)
		return true;
	else if (a.rz > b.rz)
		return false;
	else {
		if (a.ry < b.ry)
			return true;
		else if (a.ry > b.ry)
			return false;

		if (a.rx < b.rx)
			return true;
		else
			return false;
	}
}


bool cmp_kp_orig(const Keypoint &a, const Keypoint &b) {
	if (a.z < b.z)
		return true;
	else if (a.z > b.z)
		return false;
	else {
		if (a.y < b.y)
			return true;
		else if (a.y > b.y)
			return false;

		if (a.x < b.x)
			return true;
		else
			return false;
	}
}


void time_info(const char* info)
{
	static double complete_start = omp_get_wtime();
	static double total_time = -1.0f;
	static double init = -1.0f;

	//checkerror(info);
	double now = omp_get_wtime();
	if (init < 0){
		init = now;
		total_time = now;
	}

	//image finish
	if (info[0] == '@'){
		cout << "\ttotal time:" << (now - total_time) << "s  ----" << info + 1 << endl;
		init = now;
		total_time = now;
		return;
	}

	cout << "\t\ttime:" << 1000 * (now - init) << "ms  ----" << info << endl;
	init = now;
}

CSIFT3D* CSIFT3DFactory::CreateCSIFT3D(
	float* volume, int x_dim, int y_dim, int z_dim,
	int num_kp_levels, float sigma_default, float sigma_n_default, float peak_thresh, float max_eig_thres, float corner_thresh )
{
	return (new CSIFT3D(volume, x_dim, y_dim, z_dim,
	num_kp_levels, sigma_default, sigma_n_default,
	peak_thresh, max_eig_thres, corner_thresh));
}

CSIFT3D* CSIFT3DFactory::CreateCSIFT3D(std::string path_,
	int num_kp_levels, float sigma_default, float sigma_n_default, float peak_thresh, float max_eig_thres, float corner_thresh) {

	int m, n, p;
	float* volume = nullptr;
	ReadMatrixFromDisk(path_.c_str(), &m, &n, &p, &volume);
	cout << m << " " << n << " " << p << endl;

	CSIFT3D *CSIFT = new CSIFT3D(volume, m, n, p,
		num_kp_levels, sigma_default, sigma_n_default,
		peak_thresh, max_eig_thres, corner_thresh);
	free(volume);
	return CSIFT;
}

//--------------------------CPU SIFT 3D--------------------------------------

CSIFT3D::CSIFT3D()
{
	octave_num = 0;
	num_kp_levels = 0;
	sigma_default = 0.0;     // Scale of the base octave
	sigma_n_default = 0.0;   // Nomial scale
	peak_thresh = 0.0;       // DoG peak threshold
	max_eig_thres = 0.0;
	corner_thresh = 0.0;
}

CSIFT3D::~CSIFT3D() {
	Release_SIFT();
	if (global_descriptor != nullptr)
		free(global_descriptor);
}

CSIFT3D::CSIFT3D(float* volume, int x_dim, int y_dim, int z_dim,
	int num_kp_levels_, float sigma_default_, float sigma_n_default_, float peak_thresh_, float max_eig_thres_, float corner_thresh_)
{
	this->num_kp_levels = num_kp_levels_;
	this->sigma_default = sigma_default_;
	this->sigma_n_default = sigma_n_default_;
	this->peak_thresh = peak_thresh_;
	this->max_eig_thres = max_eig_thres_;
	this->corner_thresh = corner_thresh_;

	Host_Im.SetImageSize(x_dim, y_dim, z_dim);
	Host_Im.SetImageUnit(1.0, 1.0, 1.0);//Careful here
	Host_Im.SetImageScale(1.0);
	Host_Im.MallocArrayMemory();

	memcpy(Host_Im._Data, volume, sizeof(float)*x_dim*y_dim*z_dim);
	data_scale(Host_Im._Data, x_dim, y_dim, z_dim, 1, x_dim, x_dim * y_dim);
}

void CSIFT3D::KpSiftAlgorithm()
{
	//!- 1.Initialization
	time_info("----start");
	double pre_start = omp_get_wtime();
	Initialize();
	double pre_end = omp_get_wtime();
	m_timer.d_Allocation = pre_end - pre_start;
	cout << "Initialization is OK" << endl;
	time_info("----Init done");

	//!- 2.Build Gss scale space
	double start_GSS = omp_get_wtime();
	Build_Gaussian_Scale_Space();
	//new_write_GSS(Gss_Pyramid, octave_num, num_kp_levels + 3);
	//cout << "Build GSS is OK" << endl;
	time_info("----Build GSS");

	//!- 3.Build Dog scale space
	double start_DOG = omp_get_wtime();
	Build_DOG_Scale_Space();
	//new_write_DOG(DoG_Pyramid, octave_num, num_kp_levels + 2);
	//cout << "Build DOG is OK" << endl;
	time_info("----Build DOG");

	////!- 4.Detect Keypoint
	double start_DECT = omp_get_wtime();
	Detect_KeyPoints();
	//new_write_level_Keypoint(level_extrema);
	////new_write_Keypoint(extre);
	//cout << "Detect Keypoint is OK" << endl;
	time_info("----Detect keypoint");

	cout << "After detecting keypoints, kp size is : " << extre.size() << endl;

	////!- 5.Assign Orientation
	double start_Orientation = omp_get_wtime();
	Assign_Orientation();

	cout << "After Orientation, kp size is : " << filter.size() << endl;
	//write_Rot(extre);
	//write_Rot(filter);
	//write_Str_tensor(extre);
	//write_Eig(extre);
	////new_write_Keypoint(filter);
	//cout << "Assign orientation is OK" << endl;
	time_info("----Orientation");

	////!- 6.Extract Description
	double start_Extract = omp_get_wtime();
	Extract_Description();
	time_info("----Description");
	//write_Desc(filter);
	//cout << "extract description is OK" << endl;
	//new_write_Desc(filter);

	double start_release = omp_get_wtime();

#ifndef CHECK_ENABLE
	Release_SIFT();
#endif
	double end_release = omp_get_wtime();

	m_timer.d_BuildGSS = start_DOG - start_GSS;
	m_timer.d_BuildDOG = start_DECT - start_DOG;
	m_timer.d_Detect = start_Orientation - start_DECT;
	m_timer.d_AssignOrientation = start_Extract - start_Orientation;
	m_timer.d_Extraction = start_release - start_Extract;
	m_timer.d_release = end_release - start_release;
	time_info("@finish");
}

void CSIFT3D::Initialize()
{
	//Be careful here

	//Image* im_from_dicom;
	//CMemManager<Image>::hCreatePtr(im_from_dicom, 1);
	//im_from_dicom->data = NULL;
	//init_im(im_from_dicom);

	//read_dcm_dir_cpp("D:\\sift_dvc_data\\2dicom", im_from_dicom);

	//im_scale(im_from_dicom);

	//Im2TexIm(im_from_dicom, &Host_Im);

	//Set number of octaves
	//int last_octave = (int)log2((float)get_Min(im_from_dicom->nx, im_from_dicom->ny, im_from_dicom->nz)) - 3;
	int last_octave = (int)log2((float)get_Min(this->Host_Im.GetDimX(), this->Host_Im.GetDimY(), this->Host_Im.GetDimZ())) - 3;
	octave_num = last_octave + 1;// num = last - first + 1, first = 0 

	Initialize_geometry(&mesh);

	Gss_Pyramid.resize(octave_num * (num_kp_levels + 3));
	DoG_Pyramid.resize(octave_num * (num_kp_levels + 2));

	Initialize_Pyramid(Gss_Pyramid, &Host_Im, octave_num, num_kp_levels + 3, num_kp_levels, sigma_default);

	Initialize_Pyramid(DoG_Pyramid, &Host_Im, octave_num, num_kp_levels + 2, num_kp_levels, sigma_default);

}

void CSIFT3D::Build_Gaussian_Scale_Space()
{
	int gss_level = num_kp_levels + 3;

	float *sigmas = new float[num_kp_levels + 3];//多层octave通用，暂且这么理解
	float k = pow(2.0, 1.0 / num_kp_levels); // 2^(1/s)

	float sigma_default_base = sigma_default * pow(2.0, -1.0 / 3.0);//careful
	sigmas[0] = sigma_default_base;

	float sig_prev, sig_total;
	for (int i = 1; i < num_kp_levels + 3; i++)
	{
		//每一层在上一幅图的尺度上再模糊sigma[i]的尺度
		//sig_prev 为 k^n * sigma 的尺度
		//sig_total为 k^(n+1) * sigma 的尺度
		sig_prev = pow(k, i - 1) * sigma_default_base;
		sig_total = sig_prev * k;
		sigmas[i] = sqrt(sig_total * sig_total - sig_prev * sig_prev);
	}

	//Each octave
	for (int o = 0; o < octave_num; o++)
	{
		//Each level inside the octave
 		for (int i = 0; i < num_kp_levels + 3; i++)
		{
			//Image Volume
			//第一张图，直接拷贝 or 做初始sigma的模糊
			if (o == 0 && i == 0)
			{
				float base_sigma = sqrt(sigmas[i] * sigmas[i] - sigma_n_default * sigma_n_default);
				GaussianSmooth_3D(&Host_Im, &Gss_Pyramid[0], base_sigma);
			}
			else if (i == 0)
			{
				//前一组高斯图像的倒数第三层
				//如图像下标为：
				// 0 1 2 3 4 5    octave = 0      
				// 6 7 8 9 10 11  octave = 1
				// 12 13 14 15 16 17 octave = 2
				// 第一组下标为6的图是由第0组下标为3的图降采样得到的
				// (current_octave - 1) * (intervals + 3) + intervals
				DownSample_3D(&Gss_Pyramid[(o - 1) * gss_level + num_kp_levels], &Gss_Pyramid[o * gss_level]);
			}
			else
			{
				GaussianSmooth_3D(&Gss_Pyramid[o * gss_level + i - 1], &Gss_Pyramid[o * gss_level + i], sigmas[i]);
			}
		}
	}
}


void CSIFT3D::test_build2sigma(TexImage &img) {

	int factor = 4;
	float sig_prev = 1.60, sig_total = 1.6 * factor;
	float sig = sqrtf(sig_total*sig_total - sig_prev*sig_prev);

	TexImage tmp(Gss_Pyramid[1].GetDimX(), Gss_Pyramid[1].GetDimY(), Gss_Pyramid[1].GetDimZ());
	tmp.SetImageUnit(Gss_Pyramid[1].GetUnitX(), Gss_Pyramid[1].GetUnitY(), Gss_Pyramid[1].GetUnitZ());
	tmp.MallocArrayMemory();
	GaussianSmooth_3D(&Gss_Pyramid[1], &tmp, sig);


	TexImage tmp2(Gss_Pyramid[1].GetDimX() / 2 , Gss_Pyramid[1].GetDimY() / 2, Gss_Pyramid[1].GetDimZ() / 2);
	tmp2.SetImageUnit(Gss_Pyramid[1].GetUnitX() * 2, Gss_Pyramid[1].GetUnitY() * 2, Gss_Pyramid[1].GetUnitZ() * 2);
	tmp2.MallocArrayMemory();
	DownSample_3D(&tmp, &tmp2);

	img.ReSetImageSize(Gss_Pyramid[1].GetDimX() / factor, Gss_Pyramid[1].GetDimY() / factor, Gss_Pyramid[1].GetDimZ() / factor);
	img.SetImageUnit(Gss_Pyramid[1].GetUnitX() * factor, Gss_Pyramid[1].GetUnitY() * factor, Gss_Pyramid[1].GetUnitZ() * factor);
	img.MallocArrayMemory();
	DownSample_3D(&tmp2, &img);
}


void CSIFT3D::Build_DOG_Scale_Space()
{
	int gss_idx = -1;
	int dog_idx = -1;

	for (int o = 0; o < octave_num; o++)
	{
		for (int i = 1; i < num_kp_levels + 3; i++)
		{
			gss_idx = o * (num_kp_levels + 3) + i;
			dog_idx = o * (num_kp_levels + 2) + i - 1;
			Sub(&Gss_Pyramid[gss_idx - 1], &Gss_Pyramid[gss_idx], &DoG_Pyramid[dog_idx]);
		}
	}
}

void CSIFT3D::Detect_KeyPoints()
{
	int dog_interval = num_kp_levels + 2; //Range 0-4, select 1-3

	float dog_max = 0.0f;

	std::vector<Keypoint> level;
	printf("Detect Keypoint initialization done\n\n");

	int detect_num = 0;

	for (int o = 0; o < octave_num; o++)
	{
		//第一层和最后一层忽略
		for (int i = 1; i < dog_interval - 1; i++)
		{
			int idx = o * dog_interval + i;

			level.clear();

			TexImage* cur = &DoG_Pyramid[idx];
			//Search max
			dog_max = im_max_abs(&DoG_Pyramid[idx]);
			float thres = this->peak_thresh * dog_max;
			std::cerr << "[CPU] thresh*peak of octave:" << o << " level:" << i << " is " << thres << std::endl;

			for (int z = IMG_BORDER; z < cur->GetDimZ() - IMG_BORDER; z++)
			{
				for (int y = IMG_BORDER; y < cur->GetDimY() - IMG_BORDER; y++)
				{
					for (int x = IMG_BORDER; x < cur->GetDimX() - IMG_BORDER; x++)
					{
						float val = cur->GetImageDataWithIdx(x, y, z);
						//float val = cur->data[ELT2((x), (y), (z), x_stride, y_stride, z_stride)];
						
						if ((val > thres || val < -thres) &&
							IsExtrema_neighbor(&DoG_Pyramid[idx - 1], &DoG_Pyramid[idx],
								&DoG_Pyramid[idx + 1], x, y, z))
						{
							Keypoint new_kp;
							new_kp.x = x;
							new_kp.y = y;
							new_kp.z = z;
							new_kp.octave = o;
							new_kp.level = i;
							new_kp.scale = cur->GetScale();
							Initialize_Keypoint(new_kp);

							extre.push_back(new_kp);
							//for debug
							level.push_back(new_kp);
							detect_num++;
						}
					}
				}
			}

			level_extrema.push_back(level);
		}
	}

	std::cout << "DOG detected keypoint num:" << detect_num <<std::endl;

}

void CSIFT3D::Assign_Orientation()
{
	int num_week_gra = 0, num_non_dinstinct = 0, num_large_cornor = 0;
	int *RET = new int[extre.size()];

#pragma omp parallel for schedule(dynamic) num_threads(sift_thread_num)
	for (int i = 0; i < extre.size(); i++)
	{
		int level_idx;
		int res;

		Keypoint& kp = extre[i];
		level_idx = kp.octave * (num_kp_levels + 3) + kp.level;
		TexImage* level = &Gss_Pyramid[level_idx];

		const float sigma = ori_sig_fctr * kp.scale;

		res = Assign_Orientation_Imp(kp, level, sigma, max_eig_thres, corner_thresh);
		RET[i] = res;
		if (res < 1)
		{
			//sift3D->filter.push_back(kp);
			kp.x = kp.y = kp.z = -1.0;
			/*for (int j = 0; j < 9; j++)
			{
				kp.Rotation[j] = 0;
			}*/
		}
	}
#pragma omp barrier


	for (int i = 0; i < extre.size(); i++)
	{
		Keypoint& kp = extre[i];

		if (kp.x < 0) 
			continue;
		filter.push_back(kp);
	}

	for (int i = 0; i < extre.size(); ++i) {

		if (RET[i] == -1)
			num_week_gra++;
		else if (RET[i] == -2)
			num_non_dinstinct++;
		else if (RET[i] == -3)
			num_large_cornor++;
	}

	std::cout << "Num of KP filtered out: " << num_week_gra << "(weak gradient)\t" << num_non_dinstinct << "(non distinct eigen)\t"<< num_large_cornor <<"(large cornor angel)" << std::endl;
	std::cout << "After Orientation and Cornor Detection, Key points num:" << filter.size() << std::endl;

	delete RET;
}

void CSIFT3D::Extract_Description()
{
	global_descriptor = (float *)calloc(DESC_NUMEL*filter.size(),sizeof(float));

#pragma omp parallel for schedule(dynamic) num_threads(sift_thread_num)
	for (int i = 0; i < filter.size(); i++)
	{
		int level_idx = 0;
		int oct, lev;

		Keypoint& kp = filter[i];
		kp.desc = global_descriptor + i*DESC_NUMEL ;
		oct = kp.octave;
		lev = kp.level;
		level_idx = oct * (num_kp_levels + 3) + lev;
		TexImage* level = &Gss_Pyramid[level_idx];
		Extract_Descriptor_Imp(kp, level, &mesh);
	}
}


//Interface function for algorithm implementation
void DownSample_3D(TexImage * src, TexImage * dst)
{
	//Problem may cause here about the memory usage
	//隔点采样，三个维度均为原来的一半
	int nx = dst->GetDimX();
	int ny = dst->GetDimY();
	int nz = dst->GetDimZ();

#pragma omp parallel for
	for (int k = 0; k < nz; k++)
	{
		for (int m = 0; m < ny; m++)
		{
			for (int n = 0; n < nx; n++)
			{
				int src_x, src_y, src_z;
				float tmp;

				src_x = n << 1;
				src_y = m << 1;
				src_z = k << 1;

				tmp = src->GetImageDataWithIdx(src_x, src_y, src_z);
				dst->SetImageDataWithIdx(tmp, n, m, k);
			}
		}
	}
}

void GaussianSmooth_3D(TexImage * src, TexImage * dst, float sigma)
{
	//-----------------------------------------------------------------------------------
	//生成所需要的高斯核权重数组

	//确保sigma为正数
	sigma = sigma > 0 ? sigma : 0;

	//高斯核矩阵的大小为(6*sigma+1)*(6*sigma+1)
	//ksize为奇数
	//int ksize = round(sigma * 3) * 2 + 1;
	const int half_width = sigma > 0 ? max((int)ceil(sigma * 3.0), 1) : 1;
	const int width = 2 * half_width + 1;

	//if (width == 1)
	//{
	//	CopyVolumeImage(src, dst);
	//	return;
	//}

	//计算高斯核矩阵， 1D高斯核
	float *kernel = new float[width];
	float acc = 0;
	int i;
	float x;
	for (i = 0; i < width; i++)
	{
		x = i - half_width;
		x /= sigma + DBL_EPSILON;
		*(kernel + i) = exp(-0.5 * x * x);
		acc += *(kernel + i);
	}

	//归一化
	for (i = 0; i < width; i++)
	{
		*(kernel + i) /= acc;
	}

	//---------------------------------------------------------------------------------------------------------
	//！-利用1D高斯核分别进行x,y,z方向的高斯模糊
	//! -Be careful for the unit
	// Temp Variable for gaussian process
	//Image* Temp;

	TexImage tmp;
	tmp.SetImageSize(src->GetDimX(), src->GetDimY(), src->GetDimZ());
	tmp.SetImageUnit(src->GetUnitX(), src->GetUnitY(), src->GetUnitZ());
	tmp.MallocArrayMemory();
	if (tmp._Data == NULL)
		std::cerr << "Wrong At Allocating Temp Memory For SIFT" << std::endl;

	//vector<TexImage> Temp;
	//for (int i = 0; i < 3; i++)
	//{
	//	TexImage tmp;
	//	Temp.push_back(tmp);
	//}

	//for (int i = 0; i < 3; i++)
	//{
	//	Temp[i].SetImageSize(src->GetDimX(), src->GetDimY(), src->GetDimZ());
	//	Temp[i].SetImageUnit(src->GetUnitX(), src->GetUnitY(), src->GetUnitZ());
	//	Temp[i].MallocArrayMemory();
	//	if (Temp[i]._Data == NULL)
	//		std::cerr<< "Wrong At Allocating Temp Memory For SIFT" << std::endl;
	//}


	//TexImage tmp(src->GetDimX(), src->GetDimY(), src->GetDimZ());
	//tmp.MallocArrayMemory();
	//tmp.SetImageUnit(src->GetUnitX(), src->GetUnitY(), src->GetUnitZ());

	// X -> Y_transpose -> Y -> Y_transpose_back -> Z_transpose -> Z -> Z_transpose_back
	GaussianSmooth_3D_Imp(src, dst, 0, src->GetUnitX(), kernel, width);
	Im_permute(dst, &tmp, 0, 1);

	GaussianSmooth_3D_Imp(&tmp, dst, 0, src->GetUnitY(), kernel, width);
	Im_permute(dst, &tmp, 0, 1);

	Im_permute(&tmp, dst, 0, 2);
	GaussianSmooth_3D_Imp(dst, &tmp, 0, src->GetUnitZ(), kernel, width);
	Im_permute(&tmp, dst, 0, 2);

	delete[] kernel;
	//Temp.clear();
	
}

void GaussianSmooth_3D_Imp(TexImage * src, TexImage * dst, int dim, float unit, float* weight, int width)
{
	//姑且认为unit是有用的
	int half_width = width / 2;
	int nx = src->GetDimX();
	int ny = src->GetDimY();
	int nz = src->GetDimZ();

	dst->SetImageSize(nx, ny, nz);
	dst->SetImageUnit(src->GetUnitX(), src->GetUnitY(), src->GetUnitZ());
	//dst->SetImageScale(src->GetScale());

	float conv_eps = 0.1f;

	/*----------------------------------------------------------------------------------*/
	//\original implementations, from Rister
	//\this unit_factor is physical coordinate based, assuming 
	//float unit_factor = 1.0 / unit;
	/*----------------------------------------------------------------------------------*/
	

	/*----------------------------------------------------------------------------------*/
	//\corrected implementations(I think)
	//\Due to same gaussian kernels is used in multiple octaves, actually different octave has different sigma parameters stack
	//\When using same kernels, should guarantee it's performed on Image Coordinates rather than physical coordinates.
	//\Hence, the 1.0 (i.e image coordinate based) is used.
	float unit_factor = 1.0;
	//-----------------------------------------------//

	int unit_half_width = static_cast<int>(ceilf(static_cast<float>(half_width) * unit_factor));

	int dim_end = nx - 1;

	int start[] = { 0, 0, 0 };
	int end[] = { nx - 1, ny - 1, nz - 1 };

	start[dim] += unit_half_width;
	end[dim] -= (unit_half_width + 1);

	//int x_stride = 1;
	//int y_stride = nx;
	//int z_stride = nx * ny;

	//Set units, strides and dimensions for the dst, the same as src
	//dst->nx = src->nx;
	//dst->ny = src->ny;
	//dst->nz = src->nz;

	//dst->ux = src->ux;
	//dst->uy = src->uy;
	//dst->uz = src->uz;

	//dst->xs = 1;
	//dst->ys = src->nx;
	//dst->zs = src->nx * src->ny;


	//First pass: process the interior
#pragma omp parallel for
	for (int z = start[2]; z <= end[2]; z++)
	{
		for (int y = start[1]; y <= end[1]; y++)
		{
			for (int x = start[0]; x <= end[0]; x++)
			{
				int d;
				float frac;
				float tmp_low, tmp_high;
				float coords[] = { x,y,z };
				float tmp_res = 0.0;
				for (d = -half_width; d <= half_width; d++)
				{
					const float tap = weight[d + half_width];
					const float step = d * unit_factor;

					coords[dim] -= step;

					int idx_lo[] = { coords[0], coords[1], coords[2] };
					int idx_hi[] = { idx_lo[0], idx_lo[1], idx_lo[2] };
					idx_hi[dim] += 1;
					frac = coords[dim] - (float)idx_lo[dim];

					tmp_low = src->GetImageDataWithIdx(idx_lo[0], idx_lo[1], idx_lo[2]);
					tmp_high = src->GetImageDataWithIdx(idx_hi[0], idx_hi[1], idx_hi[2]);
					tmp_res += tap * ((1.0f - frac) * tmp_low + frac* tmp_high);

					//dst->data[ELT2(x, y, z, x_stride, y_stride, z_stride)] +=
					//	tap * ((1.0 - frac) * src->data[ELT2(idx_lo[0], idx_lo[1], idx_lo[2], x_stride, y_stride, z_stride)] +
					//		frac * src->data[ELT2(idx_hi[0], idx_hi[1], idx_hi[2], x_stride, y_stride, z_stride)]);

					coords[dim] += step;
				}
				dst->SetImageDataWithIdx(tmp_res, x, y, z);
			}
		}
	}

	//Second pass: process the boundaries
	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				int d;
				float frac;
				float tmp_low, tmp_high;
				const int i_coords[] = { x,y,z };
				float tmp_res = 0.0f;
				if (i_coords[dim] >= start[dim] && i_coords[dim] <= end[dim])
					continue;

				for (d = -half_width; d <= half_width; d++)
				{
					float coords[] = { x,y,z };
					const float tap = weight[d + half_width];
					const float step = d * unit_factor;

					coords[dim] -= step;

					//Mirror coordinates
					//Bug fix for using (int), leading to cutting to zero-side for negative points.
					//example: int(-0.5) === 0
					if (coords[dim] < 0)
					{
						coords[dim] = -1 * coords[dim];
					}
					else if (coords[dim] >= dim_end)
					{
						//dim_end - (coords[dim] - dim_end) - conv_eps
						coords[dim] = 2 * dim_end - coords[dim] - conv_eps;
					}
					//Explicit cast, more readable
					int idx_lo[] = { static_cast<int>(coords[0]), static_cast<int>(coords[1]), static_cast<int>(coords[2]) };
					int idx_hi[] = { idx_lo[0], idx_lo[1], idx_lo[2] };
					idx_hi[dim] += 1;
					frac = coords[dim] - (float)idx_lo[dim];

					tmp_low = src->GetImageDataWithIdx(idx_lo[0], idx_lo[1], idx_lo[2]);
					tmp_high = src->GetImageDataWithIdx(idx_hi[0], idx_hi[1], idx_hi[2]);
					tmp_res += tap * ((1.0f - frac) * tmp_low + frac* tmp_high);

#ifdef __JR_DEBUG
					if (x == 969 && y == 251 && z == 290) {
						printf("(%d,%d,%d : %f):%f\n", x, y, z, coords[dim], tap * ((1.0f - frac) * tmp_low + frac* tmp_high));
					}
#endif

					//dst->data[ELT2(x, y, z, x_stride, y_stride, z_stride)] +=
					//	tap * ((1.0 - frac) * src->data[ELT2(idx_lo[0], idx_lo[1], idx_lo[2], x_stride, y_stride, z_stride)] +
					//		frac * src->data[ELT2(idx_hi[0], idx_hi[1], idx_hi[2], x_stride, y_stride, z_stride)]);

					//Attention!

				}
#ifdef __JR_DEBUG
				if (x == 969 && y == 251 && z == 290) {
					printf("(%d,%d,%d):%f\n", x, y, z, tmp_res);
				}
#endif 

				dst->SetImageDataWithIdx(tmp_res, x, y, z);
			}
		}
	}

}

void Im_permute(TexImage * src, TexImage * dst, int dim1, int dim2)
{
	//reset the dimension and units, strides
	int dim[3] = { src->GetDimX(), src->GetDimY(), src->GetDimZ() };
	float unit[3] = { src->GetUnitX(), src->GetUnitY(), src->GetUnitZ() };

	int old_dim[3] = { src->GetDimX(), src->GetDimY(), src->GetDimZ() };
	float old_unit[3] = { src->GetUnitX(), src->GetUnitY(), src->GetUnitZ() };

	dim[dim1] = old_dim[dim2];
	dim[dim2] = old_dim[dim1];

	unit[dim1] = old_unit[dim2];
	unit[dim2] = old_unit[dim1];

	dst->SetImageSize(dim[0], dim[1], dim[2]);
	dst->SetImageUnit(unit[0], unit[1], unit[2]);
	//dst->SetImageScale(src->GetScale());

	//dst->nx = dim[0];
	//dst->ny = dim[1];
	//dst->nz = dim[2];

	//dst->ux = unit[0];
	//dst->uy = unit[1];
	//dst->uz = unit[2];

	//dst->xs = 1;
	//dst->ys = dst->nx;
	//dst->zs = dst->nx * dst->ny;

	//Transpose data
#pragma omp parallel for
	for (int z = 0; z < dst->GetDimZ(); z++)
	{
		for (int y = 0; y < dst->GetDimY(); y++)
		{
			for (int x = 0; x < dst->GetDimX(); x++)
			{
				float tmp;

				int src_coords[] = { x, y, z };
				int temp;

				temp = src_coords[dim1];
				src_coords[dim1] = src_coords[dim2];
				src_coords[dim2] = temp;

				tmp = src->GetImageDataWithIdx(src_coords[0], src_coords[1], src_coords[2]);
				dst->SetImageDataWithIdx(tmp, x, y, z);
				//dst->data[ELT2(x, y, z, dst->xs, dst->ys, dst->zs)] =
				//	src->data[ELT2(src_coords[0], src_coords[1], src_coords[2], src->xs, src->ys, src->zs)];
			}
		}
	}
}

void Sub(TexImage * prev, TexImage * cur, TexImage * dog)
{
	if (prev->GetDimX() != cur->GetDimX() || prev->GetDimX() != dog->GetDimX() ||
		prev->GetDimY() != cur->GetDimY() || prev->GetDimY() != dog->GetDimY() ||
		prev->GetDimZ() != cur->GetDimZ() || prev->GetDimZ() != dog->GetDimZ() )
	{
		std::cout << "Wrong in Sub function" << std::endl;
		return;
	}

	int x_dim = prev->GetDimX();
	int y_dim = prev->GetDimY();
	int z_dim = prev->GetDimZ();

	//int x_stride = prev->xs;
	//int y_stride = prev->ys;
	//int z_stride = prev->zs;

#pragma omp parallel for
	for (int z = 0; z < z_dim; z++)
	{
		for (int y = 0; y < y_dim; y++)
		{
			for (int x = 0; x < x_dim; x++)
			{
				float tmp;
				tmp = (cur->GetImageDataWithIdx(x, y, z) - prev->GetImageDataWithIdx(x, y, z)) * (-1);
				dog->SetImageDataWithIdx(tmp, x, y, z);
				//dog->data[ELT2(x, y, z, x_stride, y_stride, z_stride)] =
				//	(cur->data[ELT2(x, y, z, x_stride, y_stride, z_stride)] - prev->data[ELT2(x, y, z, x_stride, y_stride, z_stride)]) * (-1);
			}
		}
	}
}

bool IsExtrema_neighbor(TexImage* prev, TexImage* cur, TexImage* next, int x, int y, int z)
{
	float val = cur->GetImageDataWithIdx(x, y, z);
	float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;

	tmp0 = prev->GetImageDataWithIdx(x, y, z);
	tmp1 = cur->GetImageDataWithIdx(x - 1, y, z);
	tmp2 = cur->GetImageDataWithIdx(x + 1, y, z);
	tmp3 = cur->GetImageDataWithIdx(x, y + 1, z);
	tmp4 = cur->GetImageDataWithIdx(x, y - 1, z);
	tmp5 = cur->GetImageDataWithIdx(x, y, z + 1);
	tmp6 = cur->GetImageDataWithIdx(x, y, z - 1);
	tmp7 = next->GetImageDataWithIdx(x, y, z);

	if ((val < tmp0 && val < tmp1 && val < tmp2 && val < tmp3 &&
		val < tmp4 && val < tmp5 && val < tmp6 && val < tmp7)
		||
		(val > tmp0 && val > tmp1 && val > tmp2 && val > tmp3 &&
			val > tmp4 && val > tmp5 && val > tmp6 && val > tmp7))
	{
		return true;
	}
	else
	{
		return false;
	}

}

int Assign_Orientation_Imp(Keypoint& kp, TexImage* gaussian,const float sigma, const float max_eig_ratio, const float corner_thresh)
{
	const float win_radius = sigma * ori_rad_fctr; //not keypoint scale in original implementation, need to be aware here, see paper

	Cvec vcenter = { kp.x, kp.y, kp.z };
	Cvec vd, vd_win, vdisp, vr;//vd->gradient, vd_win-> window gradient / direction directive
	float sq_disp, weight;
	float* str_tensor = kp.str_tensor;

	int nx = gaussian->GetDimX();
	int ny = gaussian->GetDimY();
	int nz = gaussian->GetDimZ();

	//int x_stride = gaussian->xs;
	//int y_stride = gaussian->ys;
	//int z_stride = gaussian->zs;

	vd_win.x = 0.0;
	vd_win.y = 0.0;
	vd_win.z = 0.0;

	//Remain sphere size
	const float uxf = (float)gaussian->GetUnitX();
	const float uyf = (float)gaussian->GetUnitY();
	const float uzf = (float)gaussian->GetUnitZ();

	int x_start = floorf(vcenter.x - win_radius / uxf);//scale the radius
	x_start = x_start > 1 ? x_start : IMG_BORDER;

	int x_end = ceilf(vcenter.x + win_radius / uxf);
	x_end = x_end < (nx - 2) ? x_end : nx - 1 - IMG_BORDER;

	int y_start = floorf(vcenter.y - win_radius / uyf);//scale the radius
	y_start = y_start > 1 ? y_start : IMG_BORDER;

	int y_end = ceilf(vcenter.y + win_radius / uyf);
	y_end = y_end < (ny - 2) ? y_end : ny - 1 - IMG_BORDER;

	int z_start = floorf(vcenter.z - win_radius / uzf);//scale the radius
	z_start = z_start > 1 ? z_start : IMG_BORDER;

	int z_end = ceilf(vcenter.z + win_radius / uzf);
	z_end = z_end < (nz - 2) ? z_end : nz - 1 - IMG_BORDER;

	//loop sphere, construct the structure tensor
	for (int z = z_start; z <= z_end; z++)
	{
		for (int y = y_start; y <= y_end; y++)
		{
			for (int x = x_start; x <= x_end; x++)
			{
				vdisp.x = ((float)x - vcenter.x) * uxf;
				vdisp.y = ((float)y - vcenter.y) * uyf;
				vdisp.z = ((float)z - vcenter.z) * uzf;
				sq_disp = SIFT3D_CVEC_L2_NORM_SQ(vdisp);
				if (sq_disp > win_radius * win_radius)
					continue;

				weight = expf(-0.5 * sq_disp / (sigma * sigma));

				vd.x = 0.5 * (gaussian->GetImageDataWithIdx(x + 1, y, z) - gaussian->GetImageDataWithIdx(x - 1, y, z));
				vd.y = 0.5 * (gaussian->GetImageDataWithIdx(x, y + 1, z) - gaussian->GetImageDataWithIdx(x, y - 1, z));
				vd.z = 0.5 * (gaussian->GetImageDataWithIdx(x, y, z + 1) - gaussian->GetImageDataWithIdx(x, y, z - 1));

				//vd.x = 0.5 * (gaussian->data[ELT2((x + 1), (y), (z), x_stride, y_stride, z_stride)] - gaussian->data[ELT2((x - 1), (y), (z), x_stride, y_stride, z_stride)]);
				//vd.y = 0.5 * (gaussian->data[ELT2((x), (y + 1), (z), x_stride, y_stride, z_stride)] - gaussian->data[ELT2((x), (y - 1), (z), x_stride, y_stride, z_stride)]);
				//vd.z = 0.5 * (gaussian->data[ELT2((x), (y), (z + 1), x_stride, y_stride, z_stride)] - gaussian->data[ELT2((x), (y), (z - 1), x_stride, y_stride, z_stride)]);

				vd.x *= 1.0f / uxf;
				vd.y *= 1.0f / uyf;
				vd.z *= 1.0f / uzf;

				str_tensor[0] += vd.x * vd.x * weight;
				str_tensor[1] += vd.x * vd.y * weight;
				str_tensor[2] += vd.x * vd.z * weight;
				str_tensor[4] += vd.y * vd.y * weight;
				str_tensor[5] += vd.y * vd.z * weight;
				str_tensor[8] += vd.z * vd.z * weight;

				vd_win.x += vd.x * weight;
				vd_win.y += vd.y * weight;
				vd_win.z += vd.z * weight;

			}
		}
	}

	//Fill in the remaining element
	str_tensor[3] = str_tensor[1];
	str_tensor[6] = str_tensor[2];
	str_tensor[7] = str_tensor[5];

	kp.win = vd_win;

	//Reject the weak gradient
	if (SIFT3D_CVEC_L2_NORM_SQ(vd_win) < ori_grad_thresh)
	{
		printf("Reject the weak gradient\n");
		//debug for check
		return -1;
	}

	//Make the matrix eigendecomposition, use the eigen library, structure_tensor = V * D * V^(-1)
	Matrix3d str_mat(3, 3);
	str_mat(0, 0) = str_tensor[0];
	str_mat(0, 1) = str_tensor[1];
	str_mat(0, 2) = str_tensor[2];
	str_mat(1, 0) = str_tensor[3];
	str_mat(1, 1) = str_tensor[4];
	str_mat(1, 2) = str_tensor[5];
	str_mat(2, 0) = str_tensor[6];
	str_mat(2, 1) = str_tensor[7];
	str_mat(2, 2) = str_tensor[8];

	EigenSolver<Matrix3d> es(str_mat);
	Matrix3d D = es.pseudoEigenvalueMatrix();
	Matrix3d V = es.pseudoEigenvectors();

	//Be aware of the eigen value in descending order or not, make a sorting here
	//Check the eigen value and see whether they are distinct or not
	EigenVal EV[3];
	EV[0].val = D(0, 0);
	EV[1].val = D(1, 1);
	EV[2].val = D(2, 2);

	EV[0].vec[0] = V(0, 0);
	EV[0].vec[1] = V(1, 0);
	EV[0].vec[2] = V(2, 0);

	EV[1].vec[0] = V(0, 1);
	EV[1].vec[1] = V(1, 1);
	EV[1].vec[2] = V(2, 1);

	EV[2].vec[0] = V(0, 2);
	EV[2].vec[1] = V(1, 2);
	EV[2].vec[2] = V(2, 2);

	sort(EV, EV + 3, cmp);

	//debug
	kp.eigvalue[0] = EV[0].val;
	kp.eigvalue[1] = EV[1].val;
	kp.eigvalue[2] = EV[2].val;

	kp.eigvector[0] = EV[0].vec[0];
	kp.eigvector[1] = EV[0].vec[1];
	kp.eigvector[2] = EV[0].vec[2];

	kp.eigvector[3] = EV[1].vec[0];
	kp.eigvector[4] = EV[1].vec[1];
	kp.eigvector[5] = EV[1].vec[2];

	kp.eigvector[6] = EV[2].vec[0];
	kp.eigvector[7] = EV[2].vec[1];
	kp.eigvector[8] = EV[2].vec[2];

	if (fabs(EV[0].val / EV[1].val) > max_eig_ratio || fabs(EV[1].val / EV[2].val) > max_eig_ratio)
	{
		//printf("Reject the keypoint in max_eig_ratio\n");
		//debug for check
		return -2;
	}

	if (!DistinctEig(EV[0].val, EV[1].val, EV[2].val))
	{
		//debug for check
		return -2;
	}

	//Now construct the Rotation Matrix
	//Assign signs to the first n - 1 vectors
	float d, cos_ang, abs_cos_ang, q_NORM, sgn, corner_score;
	float d_NORM = SIFT3D_CVEC_L2_NORM(vd_win);
	Cvec eig_val;

	corner_score = FLT_MAX;
	for (int i = 2; i > 0; i--)
	{
		eig_val.x = EV[i].vec[0];
		eig_val.y = EV[i].vec[1];
		eig_val.z = EV[i].vec[2];
		d = SIFT3D_CVEC_DOT(eig_val, vd_win);
		q_NORM = SIFT3D_CVEC_L2_NORM(eig_val);

		cos_ang = d / (d_NORM * q_NORM);
		abs_cos_ang = fabs(cos_ang);

		//seems rejection should be here
		corner_score = corner_score < abs_cos_ang ? corner_score : abs_cos_ang;

		sgn = d > 0.0 ? 1.0 : -1.0;

		EV[i].vec[0] *= sgn;
		EV[i].vec[1] *= sgn;
		EV[i].vec[2] *= sgn;
	}

	if (corner_score < corner_thresh) {
		//debug for check
		return -3;
	}

	// Take the cross product of the first two vectors
	Cvec v1, v2;
	v1.x = EV[2].vec[0];
	v1.y = EV[2].vec[1];
	v1.z = EV[2].vec[2];

	v2.x = EV[1].vec[0];
	v2.y = EV[1].vec[1];
	v2.z = EV[1].vec[2];
	SIFT3D_CVEC_CROSS(v1, v2, vr);

	//Finally construct the Rotation Matrix
	kp.Rotation[0] = v1.x;
	kp.Rotation[1] = v2.x;
	kp.Rotation[2] = vr.x;
	kp.Rotation[3] = v1.y;
	kp.Rotation[4] = v2.y;
	kp.Rotation[5] = vr.y;
	kp.Rotation[6] = v1.z;
	kp.Rotation[7] = v2.z;
	kp.Rotation[8] = vr.z;

	return 1;
}

bool DistinctEig(float a, float b, float c)
{
	if (fabs(a - b) < DBL_EPSILON)
		return false;
	else if (fabs(a - c) < DBL_EPSILON)
		return false;
	else if (fabs(c - b) < DBL_EPSILON)
		return false;

	return true;
}

void Extract_Descriptor_Imp(Keypoint& kp, TexImage * gaussian, Mesh* mesh)
{
	//Basic parameter, implementation detail
	const float sigma = kp.scale * desc_sig_fctr;
	const float win_radius = desc_rad_fctr * sigma;
	const float desc_hw = win_radius / sqrt(2); //half width
	const float desc_width = 2.0f * desc_hw;//whole width
	const float desc_bin_fctr = (float)NHIST_PER_DIM / desc_width; //???

	const float coord_factor = pow(2.0, kp.octave);


	Cvec vcenter = { kp.x, kp.y, kp.z };

	const int height = gaussian->GetDimY();
	const int width = gaussian->GetDimX();

	int nx = gaussian->GetDimX();
	int ny = gaussian->GetDimY();
	int nz = gaussian->GetDimZ();

	//int x_stride = gaussian->xs;
	//int y_stride = gaussian->ys;
	//int z_stride = gaussian->zs;

	//Remain sphere size
	const float uxf = (float)gaussian->GetUnitX();
	const float uyf = (float)gaussian->GetUnitY();
	const float uzf = (float)gaussian->GetUnitZ();

	int x_start = floorf(vcenter.x - win_radius / uxf);//scale the radius
	x_start = x_start > 1 ? x_start : IMG_BORDER;

	int x_end = ceilf(vcenter.x + win_radius / uxf);
	x_end = x_end < (nx - 2) ? x_end : nx - 1 - IMG_BORDER;

	int y_start = floorf(vcenter.y - win_radius / uyf);//scale the radius
	y_start = y_start > 1 ? y_start : IMG_BORDER;

	int y_end = ceilf(vcenter.y + win_radius / uyf);
	y_end = y_end < (ny - 2) ? y_end : ny - 1 - IMG_BORDER;

	int z_start = floorf(vcenter.z - win_radius / uzf);//scale the radius
	z_start = z_start > 1 ? z_start : IMG_BORDER;

	int z_end = ceilf(vcenter.z + win_radius / uzf);
	z_end = z_end < (nz - 2) ? z_end : nz - 1 - IMG_BORDER;

	//Calculate the distance to confirm whether the point is in the sphere
	float weight;
	float sq_disp;

	float* R = kp.Rotation;
	Cvec vrot, vbins;
	Cvec vdisp, grad, grad_rot;

	//Save keypoint location in the original image
	/*kp.rx = kp.x * coord_factor;
	kp.ry = kp.y * coord_factor;
	kp.rz = kp.z * coord_factor;*/

	//Invert the rotation matrix
	Transpose_Matrix(kp.Rotation);

	//debug mode
	const int debug = 0;

	/*if (sigma > 280)
	{
		debug = 1;
		printf("\n\nNow Debug: x:%lf, y:%lf, z:%lf, sigma:%lf\n", kp.x, kp.y, kp.z, kp.scale);
	}*/

	int batchsize = 1;
	if (debug > 0)
	{
		batchsize = (x_end - x_start + 1) * (y_end - y_start + 1) * (z_end - z_start + 1);
	}

	//record all the middle result
	float* host_loop_point = nullptr;
	float* host_vrot = nullptr;
	float* host_vbins = nullptr;
	int* host_intersect_id = nullptr;
	float* host_bary = nullptr;
	float* host_dvbins = nullptr;
	int* host_offset = nullptr;
	float* host_desc_accum = nullptr;


	CMemManager<float>::hCreatePtr(host_loop_point, batchsize * 3);
	CMemManager<float>::hCreatePtr(host_vrot, batchsize * 3);
	CMemManager<float>::hCreatePtr(host_vbins, batchsize * 3);
	CMemManager<int>::hCreatePtr(host_intersect_id, batchsize);
	CMemManager<float>::hCreatePtr(host_bary, batchsize * 3);
	CMemManager<float>::hCreatePtr(host_dvbins, batchsize * 3);
	CMemManager<int>::hCreatePtr(host_offset, batchsize * 3 * 8);
	CMemManager<float>::hCreatePtr(host_desc_accum, batchsize * 3 * 8);


	int loop_idx = 0;
	for (int z = z_start; z <= z_end; z++)
	{
		for (int y = y_start; y <= y_end; y++)
		{
			for (int x = x_start; x <= x_end; x++)
			{
				if (debug > 0)
				{
					host_loop_point[loop_idx * 3 + 0] = x;
					host_loop_point[loop_idx * 3 + 1] = y;
					host_loop_point[loop_idx * 3 + 2] = z;
				}

				vdisp.x = ((float)x - vcenter.x) * uxf;
				vdisp.y = ((float)y - vcenter.y) * uyf;
				vdisp.z = ((float)z - vcenter.z) * uzf;
				sq_disp = SIFT3D_CVEC_L2_NORM_SQ(vdisp);
				if (sq_disp > win_radius * win_radius)
					continue;

				//Rotate to keypoint space
				vrot.x = R[0] * vdisp.x + R[1] * vdisp.y + R[2] * vdisp.z;
				vrot.y = R[3] * vdisp.x + R[4] * vdisp.y + R[5] * vdisp.z;
				vrot.z = R[6] * vdisp.x + R[7] * vdisp.y + R[8] * vdisp.z;

				//floating point locations of index of cube, 4x4x4
				vbins.x = (vrot.x + desc_hw) * desc_bin_fctr;
				vbins.y = (vrot.y + desc_hw) * desc_bin_fctr;
				vbins.z = (vrot.z + desc_hw) * desc_bin_fctr;

				/*----------------------------------------------------------------------------------*/
				//\original implementations, from Rister
				//if (vbins.x < 0 || vbins.y < 0 || vbins.z < 0 ||
				//	vbins.x >= 4 || vbins.y >= 4 || vbins.z >= 4)
				//	continue;
				/*----------------------------------------------------------------------------------*/

				/*----------------------------------------------------------------------------------*/
				//\modification, to make sure same with the paper
				//let the keypoint located at [1.5, 1.5, 1.5]
				//the left most voxel located at [-0.5, -0.5, -0.5]
				//the right most voxel located at [3.5, 3.5, 3.5]
				//centers of sub-regions are at {0, 1, 2, 3} x {0, 1, 2, 3} x {0, 1, 2, 3}
				/*----------------------------------------------------------------------------------*/
				vbins.x -= 0.5f;
				vbins.y -= 0.5f;
				vbins.z -= 0.5f;
				if (vbins.x <= -0.5f || vbins.y <= -0.5f || vbins.z <= -0.5f ||
					vbins.x >= 3.5f || vbins.y >= 3.5f || vbins.z >= 3.5f )
					continue;

				if (debug > 0)
				{
					host_vbins[loop_idx * 3 + 0] = vbins.x;
					host_vbins[loop_idx * 3 + 1] = vbins.y;
					host_vbins[loop_idx * 3 + 2] = vbins.z;
				}

				//Take gradient
				weight = expf(-0.5f * sq_disp / (sigma * sigma)); //Guassian weight

				grad.x = 0.5 * (gaussian->GetImageDataWithIdx(x + 1, y, z) - gaussian->GetImageDataWithIdx(x - 1, y, z));
				grad.y = 0.5 * (gaussian->GetImageDataWithIdx(x, y + 1, z) - gaussian->GetImageDataWithIdx(x, y - 1, z));
				grad.z = 0.5 * (gaussian->GetImageDataWithIdx(x, y, z + 1) - gaussian->GetImageDataWithIdx(x, y, z - 1));

				grad.x *= 1.0f / uxf;
				grad.y *= 1.0f / uyf;
				grad.z *= 1.0f / uzf;

				SIFT3D_CVEC_SCALE(grad, weight);

				//Rotate the gradient to keypoint space
				grad_rot.x = R[0] * grad.x + R[1] * grad.y + R[2] * grad.z;
				grad_rot.y = R[3] * grad.x + R[4] * grad.y + R[5] * grad.z;
				grad_rot.z = R[6] * grad.x + R[7] * grad.y + R[8] * grad.z;

				if (debug > 0)
				{
					host_vrot[loop_idx * 3 + 0] = grad_rot.x;
					host_vrot[loop_idx * 3 + 1] = grad_rot.y;
					host_vrot[loop_idx * 3 + 2] = grad_rot.z;
				}

				//Make trilinear interpolation
				//Trilinear_interpolation_over_desc(mesh, kp, vbins, grad_rot, loop_idx);

				Trilinear_interpolation_over_desc_debug(mesh, kp, vbins, grad_rot, loop_idx,
					host_dvbins, host_intersect_id, host_bary, host_offset, host_desc_accum, debug);

				loop_idx++;

			}
		}
	}



	normailize_desc(kp.desc);

	//truncate
	for (int i = 0; i < DESC_NUMEL; i++)
	{
		kp.desc[i] = (kp.desc[i] < trunc_thresh ? kp.desc[i] : trunc_thresh);
	}

	normailize_desc(kp.desc);

	/*if (debug > 0)
	{
		write_debug_middle(host_loop_point, host_vrot, host_vbins, host_intersect_id, host_bary,
			host_dvbins, host_offset, host_desc_accum, batchsize, kp, R);
	}*/

	CMemManager<float>::hDestroyPtr(host_loop_point);
	CMemManager<float>::hDestroyPtr(host_vrot);
	CMemManager<float>::hDestroyPtr(host_vbins);
	CMemManager<int>::hDestroyPtr(host_intersect_id);
	CMemManager<float>::hDestroyPtr(host_bary);
	CMemManager<float>::hDestroyPtr(host_dvbins);
	CMemManager<int>::hDestroyPtr(host_offset);
	CMemManager<float>::hDestroyPtr(host_desc_accum);


	//Save keypoint location in the original image
	kp.rx = kp.x * coord_factor;
	kp.ry = kp.y * coord_factor;
	kp.rz = kp.z * coord_factor;

}

void Trilinear_interpolation_over_desc(Mesh * mesh, Keypoint &kp, Cvec &vbins, Cvec &grad, int loop_idx)
{
	Cvec bary, dvbins;
	int dx, dy, dz, x, y, z;
	int hist_idx, offset_v1, offset_v2, offset_v3;
	float mag;
	float weight;
	int intersect_idx;

	//!-------------------------------------------------------
	const int y_stride = NHIST_PER_DIM;
	const int z_stride = NHIST_PER_DIM * NHIST_PER_DIM;
	//!-------------------------------------------------------

	//Decimal part
	dvbins.x = vbins.x - floorf(vbins.x);
	dvbins.y = vbins.y - floorf(vbins.y);
	dvbins.z = vbins.z - floorf(vbins.z);

	intersect_idx = Check_intersect_faces(mesh, &grad, &bary);
	if (intersect_idx < 0)
		return;

	//Get the magnitude of the vector
	mag = SIFT3D_CVEC_L2_NORM(grad);
	int hist_stripe = ICOS_NVERT;

	//Make the trilinear interpolation and find the corresponding bins
	for (dx = 0; dx < 2; dx++)
	{
		for (dy = 0; dy < 2; dy++)
		{
			for (dz = 0; dz < 2; dz++)
			{
				x = (int)vbins.x + dx;
				y = (int)vbins.y + dy;
				z = (int)vbins.z + dz;

				if (x < 0 || y < 0 || z < 0 ||
					x >= NHIST_PER_DIM || y >= NHIST_PER_DIM || z >= NHIST_PER_DIM)
					continue;

				//Fill the descriptor histogram
				hist_idx = x + y * y_stride + z * z_stride;

				//Get the spatial interpolation weight
				weight = ((dx == 0) ? (1.0 - dvbins.x) : dvbins.x) *
					((dy == 0) ? (1.0 - dvbins.y) : dvbins.y) *
					((dz == 0) ? (1.0 - dvbins.z) : dvbins.z);

				//Get the bin and add the weighted magnitude
				Tri * intersect_tri = mesh->tri + intersect_idx;
				offset_v1 = hist_idx * hist_stripe + intersect_tri->idx[0];
				offset_v2 = hist_idx * hist_stripe + intersect_tri->idx[1];
				offset_v3 = hist_idx * hist_stripe + intersect_tri->idx[2];

				kp.desc[offset_v1] += mag * weight * bary.x;
				kp.desc[offset_v2] += mag * weight * bary.y;
				kp.desc[offset_v3] += mag * weight * bary.z;

			}
		}
	}

}


void Trilinear_interpolation_over_desc_debug(Mesh * mesh, Keypoint &kp, Cvec &vbins, Cvec &grad, int loop_idx,
	float* host_dvbins, int* host_intersect_id, float* host_bary, int* host_offset, float* host_desc_accum, int debug)
{
	Cvec bary, dvbins;
	int dx, dy, dz, x, y, z;
	int hist_idx, offset_v1, offset_v2, offset_v3;
	float mag;
	float weight;
	int intersect_idx;

	//!-------------------------------------------------------
	const int y_stride = NHIST_PER_DIM;
	const int z_stride = NHIST_PER_DIM * NHIST_PER_DIM;
	//!-------------------------------------------------------

	//Decimal part
	dvbins.x = vbins.x - floorf(vbins.x);
	dvbins.y = vbins.y - floorf(vbins.y);
	dvbins.z = vbins.z - floorf(vbins.z);

	intersect_idx = Check_intersect_faces(mesh, &grad, &bary);

	if (debug > 0)
	{
		host_dvbins[loop_idx * 3 + 0] = dvbins.x;
		host_dvbins[loop_idx * 3 + 1] = dvbins.y;
		host_dvbins[loop_idx * 3 + 2] = dvbins.z;
		host_bary[loop_idx * 3 + 0] = bary.x;
		host_bary[loop_idx * 3 + 1] = bary.y;
		host_bary[loop_idx * 3 + 2] = bary.z;
		host_intersect_id[loop_idx] = intersect_idx;
	}

	if (intersect_idx < 0)
		return;

	//Get the magnitude of the vector
	mag = SIFT3D_CVEC_L2_NORM(grad);
	int hist_stripe = ICOS_NVERT;

	//Make the trilinear interpolation and find the corresponding bins
	int loop = 0;
	for (dx = 0; dx < 2; dx++)
	{
		for (dy = 0; dy < 2; dy++)
		{
			for (dz = 0; dz < 2; dz++)
			{
				x = (int)vbins.x + dx;
				y = (int)vbins.y + dy;
				z = (int)vbins.z + dz;

				if (x < 0 || y < 0 || z < 0 ||
					x >= NHIST_PER_DIM || y >= NHIST_PER_DIM || z >= NHIST_PER_DIM)
					continue;

				//Fill the descriptor histogram
				hist_idx = x + y * y_stride + z * z_stride;

				//Get the spatial interpolation weight
				weight = ((dx == 0) ? (1.0 - dvbins.x) : dvbins.x) *
					((dy == 0) ? (1.0 - dvbins.y) : dvbins.y) *
					((dz == 0) ? (1.0 - dvbins.z) : dvbins.z);

				//Get the bin and add the weighted magnitude
				Tri * intersect_tri = mesh->tri + intersect_idx;
				offset_v1 = hist_idx * hist_stripe + intersect_tri->idx[0];
				offset_v2 = hist_idx * hist_stripe + intersect_tri->idx[1];
				offset_v3 = hist_idx * hist_stripe + intersect_tri->idx[2];

				kp.desc[offset_v1] += mag * weight * bary.x;
				kp.desc[offset_v2] += mag * weight * bary.y;
				kp.desc[offset_v3] += mag * weight * bary.z;

				if (debug > 0)
				{
					host_offset[loop_idx * 3 * 8 + loop * 3 + 0] = offset_v1;
					host_offset[loop_idx * 3 * 8 + loop * 3 + 1] = offset_v2;
					host_offset[loop_idx * 3 * 8 + loop * 3 + 2] = offset_v3;

					host_desc_accum[loop_idx * 3 * 8 + loop * 3 + 0] = mag * weight * bary.x;
					host_desc_accum[loop_idx * 3 * 8 + loop * 3 + 1] = mag * weight * bary.y;
					host_desc_accum[loop_idx * 3 * 8 + loop * 3 + 2] = mag * weight * bary.z;
				}

				loop++;

			}
		}
	}
}

int Check_intersect_faces(Mesh *mesh, Cvec* grad, Cvec* bary)
{
	if (SIFT3D_CVEC_L2_NORM_SQ(*grad) < bary_eps)
		return -1;

	//Iterate through the faces
	int i;
	float k;
	int res;
	for (i = 0; i < ICOS_NFACES; i++)
	{
		const Tri* const tri = mesh->tri + i;

		// Convert to barycentric coordinates
		res = cart2bary(grad, tri, bary, &k);
		if (res < 0)
			continue;

		if (bary->x < -bary_eps || bary->y < -bary_eps || bary->z < -bary_eps || k < 0)
		{
			continue;
		}
		else
		{
			return i;
			//break;
		}

	}

	return -1;
}

void Transpose_Matrix(float *Rot)
{
	//Since we know it's 3x3 matrix, make an easy implementation without loop
	//Swap three pairs of elements
	Swap_Element(Rot[1], Rot[3]);
	Swap_Element(Rot[2], Rot[6]);
	Swap_Element(Rot[5], Rot[7]);
}

void Swap_Element(float &a, float &b)
{
	float tmp;
	tmp = a;
	a = b;
	b = tmp;
}

int cart2bary(Cvec * cart, const Tri * const tri, Cvec * const bary, float * const k)
{
	Cvec e1, e2, t, p, q;
	float det, det_inv;

	const Cvec * const v = tri->v;

	//Pre-compute e1, e2, t, p, q as paper illustration
	// E1 = V1 - V0
	e1.x = (v + 1)->x - v->x;
	e1.y = (v + 1)->y - v->y;
	e1.z = (v + 1)->z - v->z;

	// E2 = V2 - V0
	e2.x = (v + 2)->x - v->x;
	e2.y = (v + 2)->y - v->y;
	e2.z = (v + 2)->z - v->z;

	// T = O - V0; O = (0,0,0)
	t.x = v->x * (-1.0);
	t.y = v->y * (-1.0);
	t.z = v->z * (-1.0);

	// P = D X E2 (D -> grad)
	SIFT3D_CVEC_CROSS(*cart, e2, p);

	// Q = T X E1
	SIFT3D_CVEC_CROSS(t, e1, q);

	det = SIFT3D_CVEC_DOT(e1, p);

	if (fabsf(det) < bary_eps)
	{
		return -1;
	}

	det_inv = 1.0 / det;

	bary->y = det_inv * SIFT3D_CVEC_DOT(p, t);
	bary->z = det_inv * SIFT3D_CVEC_DOT(*cart, q);
	bary->x = 1 - bary->y - bary->z;

	*k = det_inv * SIFT3D_CVEC_DOT(q, e2);

	return 0;
}

void normailize_desc(float* desc)
{
	float norm = 0.0;
	float tmp;
	for (int i = 0; i < DESC_NUMEL; i++)
	{
		tmp = *(desc + i);
		norm += tmp * tmp;
	}

	norm = sqrt(norm) + DBL_EPSILON;

	for (int i = 0; i < DESC_NUMEL; i++)
	{
		float norm_inv = 1.0 / norm;
		*(desc + i) *= norm_inv;
	}
}


void CSIFT3D::Release_SIFT()
{
	for (int i = 0; i < octave_num * (num_kp_levels + 3) && i < Gss_Pyramid.size(); i++)
	{
		CMemManager<float>::hDestroyPtr(Gss_Pyramid[i]._Data);
	}
	Gss_Pyramid.clear();
	Gss_Pyramid.shrink_to_fit();

	for (int i = 0; i < octave_num * (num_kp_levels + 2) && i < DoG_Pyramid.size(); i++)
	{
		CMemManager<float>::hDestroyPtr(DoG_Pyramid[i]._Data);
	}
	DoG_Pyramid.clear();
	DoG_Pyramid.shrink_to_fit();

	extre.clear();
	extre.shrink_to_fit();

}

void CSIFT3D::SetNumThreads(int t_num) 
{
	if (t_num > 0)
		sift_thread_num = t_num;
}

std::vector<Keypoint> CSIFT3D::GetKeypoints() {
	return filter;
}

}
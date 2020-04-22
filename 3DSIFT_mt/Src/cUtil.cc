
#include "../../Include/cUtil.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "../../Include/Util/cMemManager.h"
//#include "cTexImage.h"
//#include "dicom.h"
using namespace std;

namespace CPUSIFT{

/* Internal math constants */
const double gr = 1.6180339887; // Golden ratio

/* Verices of a regular icosahedron inscribed in the unit sphere. */
const double vert[] = { 0,  1,  gr,
						0, -1,  gr,
						0,  1, -gr,
						0, -1, -gr,
						1,  gr,  0,
						-1,  gr,  0,
						1, -gr,  0,
						-1, -gr,  0,
						gr,   0,  1,
						-gr,   0,  1,
						gr,   0, -1,
						-gr,   0, -1 };

/* Vertex triplets forming the faces of the icosahedron. */
const int faces[] = { 0, 1, 8,
						0, 8, 4,
						0, 4, 5,
						0, 5, 9,
						0, 9, 1,
						1, 6, 8,
						8, 6, 10,
						8, 10, 4,
						4, 10, 2,
						4, 2, 5,
						5, 2, 11,
						5, 11, 9,
						9, 11, 7,
						9, 7, 1,
						1, 7, 6,
						3, 6, 7,
						3, 7, 11,
						3, 11, 2,
						3, 2, 10,
						3, 10, 6 };

//int Initialize(SIFT3D *sift3D)
//{
//	//Get the input data
//	Image * im = &sift3D->Im_origin;
//
//	read_dcm_dir_cpp("D:\\sift_dvc_data\\2dicom", im);
//
//	Initialize_Parameter(sift3D);
//
//	Initialize_geometry(sift3D);
//
//	Initialize_Pyramid(sift3D);
//
//	return 0;
//}

//!-------------------------------------- Refactor ------------------------------------------------------
float SIFT3D_CVEC_L2_NORM(Cvec vec)
{
	return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

float SIFT3D_CVEC_L2_NORM_SQ(Cvec vec)
{
	return (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

float SIFT3D_CVEC_DOT(Cvec& vec1, Cvec& vec2)
{
	return (vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
}

void SIFT3D_CVEC_CROSS(Cvec& vec1, Cvec& vec2, Cvec& out)
{
	out.x = vec1.y * vec2.z - vec1.z * vec2.y;
	out.y = vec1.z * vec2.x - vec1.x * vec2.z;
	out.z = vec1.x * vec2.y - vec1.y * vec2.x;
}

void SIFT3D_CVEC_SCALE(Cvec& vec, double sca)
{
	vec.x = vec.x * sca;
	vec.y = vec.y * sca;
	vec.z = vec.z * sca;
}


void Im2TexIm(Image* im, TexImage* tex_im)
{
	tex_im->SetImageSize(im->nx, im->ny, im->nz);
	tex_im->SetImageScale(im->s);
	tex_im->SetImageUnit(1.0, 1.0, 1.0);
	tex_im->SetImageDataPt(im->data);
}


int Initialize_geometry(Mesh * mesh)
{
	double mag;
	Cvec temp1, temp2, temp3, n;

	//First init_Mesh
	mesh->tri = nullptr;
	mesh->num = ICOS_NFACES;

	//malloc the memory space for mesh
	CMemManager<Tri>::hCreatePtr(mesh->tri, ICOS_NFACES);

	//Init the element in mesh
	for (int i = 0; i < ICOS_NFACES; i++)
	{
		Tri * const tri = mesh->tri + i;
		Cvec * const v = tri->v;
		int * const idx = tri->idx;

		int face_idx = i * 3;
		idx[0] = faces[face_idx];
		idx[1] = faces[face_idx + 1];
		idx[2] = faces[face_idx + 2];

		//Fill the Triangle vertex coordinate
		for (int j = 0; j < 3; j++)
		{
			const double mag_expected = sqrt(1 + gr * gr);

			int ele_idx = idx[j] * 3;
			v[j].x = vert[ele_idx];
			v[j].y = vert[ele_idx + 1];
			v[j].z = vert[ele_idx + 2];

			//Normalize to unit length
			mag = SIFT3D_CVEC_L2_NORM(v[j]);
			/*if (fabsf(mag - mag_expected) < 1E-10)
			{
			printf("No Problem\n");
			}*/
			SIFT3D_CVEC_SCALE(v[j], 1.0 / mag);
		}

		temp1.x = (v + 2)->x - (v + 1)->x;
		temp1.y = (v + 2)->y - (v + 1)->y;
		temp1.z = (v + 2)->z - (v + 1)->z;

		temp2.x = (v + 1)->x - v->x;
		temp2.y = (v + 1)->y - v->y;
		temp2.z = (v + 1)->z - v->z;

		SIFT3D_CVEC_CROSS(temp1, temp2, n);

		if (SIFT3D_CVEC_DOT(n, *v) < 0)
		{
			temp1 = *v;
			*v = *(v + 1);
			*(v + 1) = temp1;
		}

	}
	return 0;
}

void Initialize_Pyramid(std::vector<TexImage> &Pyramid, TexImage* Base_Im,
	int num_oct, int interval, int num_kp_levels, float sigma_default)
{
	int total_level = num_oct * interval;

	double sigma0 = sigma_default * pow(2.0, -1.0 / 3.0);

	//Initialize each level separately, for gss pyramid, be careful for the scale
	int o, s;
	int idx = 0;
	double scale_factor;
	int nx = Base_Im->GetDimX();
	int ny = Base_Im->GetDimY();
	int nz = Base_Im->GetDimZ();
	float ux = Base_Im->GetUnitX();
	float uy = Base_Im->GetUnitY();
	float uz = Base_Im->GetUnitZ();

	//double sigma0 = SIGMA_DEFAULT * pow(2.0, -1.0 / 3.0);

	//int nx_back = nx;
	//int ny_back = ny;
	//int nz_back = nz;

	TexImage level(nx, ny, nz);
	for (int i = 0; i < total_level; i++)
	{
		level.ReSetImageSize(nx, ny, nz);
		//level.MallocArrayMemory();
		level.SetImageUnit(ux, uy, uz);
		o = i / interval;
		s = i % interval;
		scale_factor = pow(2.0, o + (double)s / num_kp_levels);
		level.SetImageScale(scale_factor * sigma0);

		idx++;
		//Pyramid.push_back(level);
		Pyramid[i] = level;
		if (idx == interval)
		{
			idx = 0;
			//3 dimension
			nx /= 2;
			ny /= 2;
			nz /= 2;
			ux *= 2;
			uy *= 2;
			uz *= 2;
		}
		//level.SetImageDataPt(NULL);
		//level.~TexImage();
	}

	for (int i = 0; i < total_level; i++)
	{
		Pyramid[i].MallocArrayMemory();
	}

}


//!-------------------------------------------------------------------------------------------------------
//int Initialize_geometry(SIFT3D *sift3D)
//{
//	Mesh * const mesh = &sift3D->mesh;
//	double mag;
//
//	Cvec temp1, temp2, temp3, n;
//
//
//	//First init_Mesh
//	mesh->tri = NULL;
//	mesh->num = ICOS_NFACES;
//
//	//malloc the memory space for mesh
//	CMemManager<Tri>::hCreatePtr(mesh->tri,ICOS_NFACES);
//
//	//Init the element in mesh
//	for (int i = 0; i < ICOS_NFACES; i++)
//	{
//		Tri * const tri = mesh->tri + i;
//		Cvec * const v = tri->v;
//		int * const idx = tri->idx;
//
//		int face_idx = i * 3;
//		idx[0] = faces[face_idx];
//		idx[1] = faces[face_idx + 1];
//		idx[2] = faces[face_idx + 2];
//
//		//Fill the Triangle vertex coordinate
//		for (int j = 0; j < 3; j++)
//		{
//			const double mag_expected = sqrt(1 + gr * gr);
//
//			int ele_idx = idx[j] * 3;
//			v[j].x = vert[ele_idx];
//			v[j].y = vert[ele_idx + 1];
//			v[j].z = vert[ele_idx + 2];
//
//			//Normalize to unit length
//			mag = SIFT3D_CVEC_L2_NORM(v[j]);
//			/*if (fabsf(mag - mag_expected) < 1E-10)
//			{
//				printf("No Problem\n");
//			}*/
//			SIFT3D_CVEC_SCALE(v[j], 1.0 / mag);
//		}
//
//		temp1.x = (v + 2)->x - (v + 1)->x;
//		temp1.y = (v + 2)->y - (v + 1)->y;
//		temp1.z = (v + 2)->z - (v + 1)->z;
//
//		temp2.x = (v + 1)->x - v->x;
//		temp2.y = (v + 1)->y - v->y;
//		temp2.z = (v + 1)->z - v->z;
//
//		SIFT3D_CVEC_CROSS(temp1, temp2, n);
//
//		if (SIFT3D_CVEC_DOT(n, *v) < 0)
//		{
//			temp1 = *v;
//			*v = *(v + 1);
//			*(v + 1) = temp1;
//		}
//
//	}
//	return 0;
//}
//
//
////can be refactored to each pyramid initialization
//int Initialize_Pyramid(SIFT3D * sift3D)
//{
//	//levels should be 6 images in each Gaussian octave, 5 images in each DOG octave
//	int octave, interval_gaussian, interval_dog, gss_level, dog_level;
//
//	//Be careful here
//	double sigma0 = SIGMA_DEFAULT * pow(2.0, -1.0 / 3.0);// -1 level
//
//	const Image Image_ori = sift3D->Im_origin;
//	int dim[3];
//	double unit[3];
//	dim[0] = Image_ori.nx;
//	dim[1] = Image_ori.ny;
//	dim[2] = Image_ori.nz;
//	unit[0] = unit[1] = unit[2] = 1.0;
//
//	octave = sift3D->octave_num;
//	interval_gaussian = sift3D->interval + 3;
//	interval_dog = interval_gaussian - 1;
//
//	gss_level = octave * interval_gaussian;
//	dog_level = octave * interval_dog;
//
//	//Allocate level Image memory
//	CMemManager<Image>::hCreatePtr(sift3D->Gss_Pyramid, gss_level);
//	CMemManager<Image>::hCreatePtr(sift3D->DOG_Pyramid, dog_level);
//
//	//Initialize new levels
//	for (int i = 0; i < gss_level; i++)
//	{
//		Image * level = sift3D->Gss_Pyramid + i;
//		init_im(level);
//	}
//
//	for (int i = 0; i < dog_level; i++)
//	{
//		Image * level = sift3D->DOG_Pyramid + i;
//		init_im(level);
//	}
//
//	//Initialize each level separately, for gss pyramid, be careful for the scale
//	int o, s;
//	int idx = 0;
//	double scale_factor;
//	for (int i = 0; i < gss_level; i++)
//	{
//		Image * level = sift3D->Gss_Pyramid + i;
//		level->nx = dim[0];
//		level->ny = dim[1];
//		level->nz = dim[2];
//
//		level->ux = unit[0];
//		level->uy = unit[1];
//		level->uz = unit[2];
//
//		level->xs = 1;
//		level->ys = level->nx;
//		level->zs = level->nx * level->ny;
//
//		level->size = level->nx * level->ny * level->nz;
//		CMemManager<float>::hCreatePtr(level->data, level->size);
//
//		o = i / interval_gaussian;
//		s = i % interval_gaussian;
//		scale_factor = pow(2.0, o + (double)s / NUM_KP_LEVELS);
//		level->s = sigma0 * scale_factor;
//
//		idx++;
//		if (idx == interval_gaussian)
//		{
//			idx = 0;
//			//3 dimension
//			for (int j = 0; j < 3; j++)
//			{
//				dim[j] /= 2;
//				unit[j] *= 2;
//			}
//		}
//	}
//
//	dim[0] = Image_ori.nx;
//	dim[1] = Image_ori.ny;
//	dim[2] = Image_ori.nz;
//	unit[0] = unit[1] = unit[2] = 1.0;
//
//	//DOG levels
//	idx = 0;
//	for (int j = 0; j < dog_level; j++)
//	{
//		Image * level = sift3D->DOG_Pyramid + j;
//		level->nx = dim[0];
//		level->ny = dim[1];
//		level->nz = dim[2];
//
//		level->ux = unit[0];
//		level->uy = unit[1];
//		level->uz = unit[2];
//
//		level->xs = 1;
//		level->ys = level->nx;
//		level->zs = level->nx * level->ny;
//
//		level->size = level->nx * level->ny * level->nz;
//		CMemManager<float>::hCreatePtr(level->data, level->size);
//
//		o = j / interval_dog;
//		s = j % interval_dog;
//		level->s = sigma0 * pow(2.0, o + (double)s / NUM_KP_LEVELS);
//
//		idx++;
//		if (idx == interval_dog)
//		{
//			idx = 0;
//			//3 dimension
//			for (int k = 0; k < 3; k++)
//			{
//				dim[k] /= 2;
//				unit[k] *= 2;
//			}
//		}
//	}
//
//	return 0;
//}
//
//int Initialize_Parameter(SIFT3D * sift3D)
//{
//	//Calculate the number of octaves
//	const Image im = sift3D->Im_origin;
//
//	//The minimum size of a pyramid level is at least 8 pixels in any dimension
//	int last_octave = (int)log2((double)get_Min(im.nx, im.ny, im.nz)) - 3;
//	int num_octaves = last_octave + 1;// num = last - first + 1, first = 0 
//
//	sift3D->octave_num = num_octaves;
//	sift3D->interval = NUM_KP_LEVELS;
//
//	return 0;
//}


void Initialize_Keypoint(Keypoint& kp)
{
	kp.rx = kp.ry = kp.rz = -1.0;
	for (int i = 0; i < 9; i++)
	{
		kp.Rotation[i] = 0;
		kp.str_tensor[i] = 0;
	}

	
	//for (int j = 0; j < 4 * 4 * 4 * 12; j++)
	//{
	//	kp.desc[j] = 0;
	//}
}

int get_Min(int row, int col, int height)
{
	int res = row < col ? row : col;
	res = res < height ? res : height;

	return res;
}

void CopyVolumeImage(Image * src, Image * dst)
{
	int x = src->nx;
	int y = src->ny;
	int z = src->nz;

	dst->nx = x;
	dst->ny = y;
	dst->nz = z;

	//hCreatePtr(dst.data, x * y* z);
	for (int i = 0; i < x*y*z; i++)
	{
		dst->data[i] = src->data[i];
	}
}

//Initialize the values of the image struct, do not allocate memory
void init_im(Image *const im)
{
	im->data = NULL;
	im->ux = 1;
	im->uy = 1;
	im->uz = 1;

	im->nx = 0;
	im->ny = 0;
	im->nz = 0;

	im->xs = 0;
	im->ys = 0;
	im->zs = 0;

	im->size = 0;
	im->s = -1.0;
}

//Clean up memory for an Image
void im_free(Image * im)
{
	if (im->data != NULL)
		free(im->data);
}

void im_scale(Image * im)
{
	double max = im_max_abs(im);
	if (max == 0.0)
		return;

	for (int z = 0; z < im->nz; z++)
	{
		for (int y = 0; y < im->ny; y++)
		{
			for (int x = 0; x < im->nx; x++)
			{
				im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)] /= max;
			}
		}
	}

}

void data_scale(float * data, int nx, int ny, int nz, size_t xs, size_t ys, size_t zs)
{
	float max = 0.0;
	float tmp;

	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				tmp = data[ELT2(x, y, z, xs, ys, zs)];
				max = (fabs(tmp) > max) ? fabs(tmp) : max;
			}
		}
	}

	for (int z = 0; z < nz; z++)
	{
		for (int y = 0; y < ny; y++)
		{
			for (int x = 0; x < nx; x++)
			{
				data[ELT2(x, y, z, xs, ys, zs)] /= max;
			}
		}
	}

}

//absolute value
float im_max_abs(Image * im)
{
	float max = 0.0;
	float tmp;

	for (int z = 0; z < im->nz; z++)
	{
		for (int y = 0; y < im->ny; y++)
		{
			for (int x = 0; x < im->nx; x++)
			{
				tmp = im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)];
				max = (fabs(tmp) > max) ? fabs(tmp) : max;
			}
		}
	}
	return max;
}


float im_max_abs(TexImage * im)
{
	float max = 0.0;
	float tmp;

	for (int z = 0; z < im->GetDimZ(); z++)
	{
		for (int y = 0; y < im->GetDimY(); y++)
		{
			for (int x = 0; x < im->GetDimX(); x++)
			{
				tmp = im->GetImageDataWithIdx(x, y, z);
				//tmp = im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)];
				max = (fabs(tmp) > max) ? fabs(tmp) : max;
			}
		}
	}
	return max;
}

//New version for debug
void new_write_Im(TexImage* im)
{
	char file_name[50];
	int nx, ny, nz;
	errno_t err;

	FILE* out_im;
	sprintf_s(file_name, "D:\\debug\\im_win_new.txt");
	nx = im->GetDimX();
	ny = im->GetDimY();
	nz = im->GetDimZ();
	float tmp;

	printf("nx:%d ny:%d nz:%d\n", nx, ny, nz);
	err = fopen_s(&out_im, file_name, "w");
	if (err == 0)
	{
		for (int z = 0; z < nz; z++)
		{
			for (int y = 0; y < ny; y++)
			{
				for (int x = 0; x < nx; x++)
				{
					tmp = im->GetImageDataWithIdx(x, y, z);
					//tmp = im.data[ELT2((x), (y), (z), im.xs, im.ys, im.zs)];
					fprintf(out_im, "%.5lf, ", tmp);

					if (x == nx - 1)
						fprintf(out_im, "\n");
				}
			}
		}
	}
	fclose(out_im);
}

void new_write_Im(TexImage* im, int idx)
{
	char file_name[100];
	int nx, ny, nz;
	errno_t err;

	FILE* out_im;
	sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_gss_%d.txt", idx);
	nx = im->GetDimX();
	ny = im->GetDimY();
	nz = im->GetDimZ();
	double tmp;

	err = fopen_s(&out_im, file_name, "w");
	if (err == 0)
	{
		for (int z = 0; z < nz; z++)
		{
			for (int y = 0; y < ny; y++)
			{
				for (int x = 0; x < nx; x++)
				{
					tmp = im->GetImageDataWithIdx(x, y, z);
					//tmp = im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)];
					fprintf(out_im, "%.5lf, ", tmp);

					if (x == nx - 1)
						fprintf(out_im, "\n");
				}
			}
		}
	}
	fclose(out_im);
}

void new_write_GSS(std::vector<TexImage> &Gss_Pyramid, int octave_num, int interval)
{
	//Batch output gss level
	for (int i = 0; i < octave_num * interval; i++)
	{
		//Image* level = Gss_level + i;
		new_write_Im(&Gss_Pyramid[i], i);
	}
}

void new_write_DOG(std::vector<TexImage> &DoG_Pyramid, int octave_num, int interval)
{
	for (int i = 0; i < octave_num * interval; i++)
	{
		//Image* level = Gss_level + i;
		new_write_Im(&DoG_Pyramid[i], i);
	}
}


void new_write_Keypoint(std::vector<Keypoint> &kp)
{
	FILE* out_im;
	char file_name[30];
	sprintf_s(file_name, "D:\\debug2\\keypoint.txt");
	errno_t err;
	err = fopen_s(&out_im, file_name, "w");
	for (int i = 0; i < kp.size(); i++)
	{
		fprintf(out_im, "%.5lf,%.5lf,%.5lf,%.5lf\n",
			kp[i].x, kp[i].y, kp[i].z, kp[i].scale);

	}
	printf("size:%zd\n", kp.size());
	fclose(out_im);
}


















////for debug and output the intermediate variable
//void write_GSS(SIFT3D * sift3D)
//{
//	Image* Gss_level = sift3D->Gss_Pyramid;
//	int gss_length = sift3D->interval + 3;
//	int octave = sift3D->octave_num;
//
//	char file_name[20];
//	int nx, ny, nz;
//	errno_t err;
//	double tmp;
//
//	//Batch output gss level
//	for (int i = 0; i < octave * gss_length; i++)
//	{
//		Image* level = Gss_level + i;
//		write_Im(level, i);
//	}
//}
//
//
//void write_Im(Image * im, int idx)
//{
//	//Image im = sift3D->Im_origin;
//
//	char file_name[30];
//	int nx, ny, nz;
//	errno_t err;
//
//	FILE* out_im;
//	sprintf_s(file_name, "D:\\debug\\im_dog_win_%d.txt", idx);
//	nx = im->nx;
//	ny = im->ny;
//	nz = im->nz;
//	double tmp;
//
//	//printf("nx:%d ny:%d nz:%d\n", nx, ny, nz);
//	err = fopen_s(&out_im, file_name, "w");
//	if (err == 0)
//	{
//		for (int z = 0; z < nz; z++)
//		{
//			for (int y = 0; y < ny; y++)
//			{
//				for (int x = 0; x < nx; x++)
//				{
//					tmp = im->data[ELT2(x, y, z, im->xs, im->ys, im->zs)];
//					fprintf(out_im, "%.5lf, ", tmp);
//
//					if (x == nx - 1)
//						fprintf(out_im, "\n");
//				}
//			}
//		}
//	}
//	fclose(out_im);
//
//}
//
//
//void write_Im(SIFT3D * sift3D)
//{
//	Image im = sift3D->Im_origin;
//
//	char file_name[30];
//	int nx, ny, nz;
//	errno_t err;
//
//	FILE* out_im;
//	sprintf_s(file_name, "D:\\debug\\im_win.txt");
//	nx = im.nx;
//	ny = im.ny;
//	nz = im.nz;
//	double tmp;
//
//	printf("nx:%d ny:%d nz:%d\n", nx, ny, nz);
//	err = fopen_s(&out_im, file_name, "w");
//	if (err == 0)
//	{
//		for (int z = 0; z < nz; z++)
//		{
//			for (int y = 0; y < ny; y++)
//			{
//				for (int x = 0; x < nx; x++)
//				{
//					tmp = im.data[ELT2((x), (y), (z), im.xs, im.ys, im.zs)];
//					fprintf(out_im, "%.5lf, ", tmp);
//
//					if (x == nx - 1)
//						fprintf(out_im, "\n");
//				}
//			}
//		}
//	}
//	fclose(out_im);
//
//}
//
//
//void write_DOG(SIFT3D * sift3D)
//{
//	Image* DOG_level = sift3D->DOG_Pyramid;
//
//	int dog_length = sift3D->interval + 2;
//	int octave = sift3D->octave_num;
//
//	char file_name[20];
//	for (int i = 0; i < octave * dog_length; i++)
//	{
//		Image* level = DOG_level + i;
//		write_Im(level, i);
//	}
//}
//
//
//void write_Keypoint(SIFT3D * sift3D)
//{
//	vector<Keypoint> kp = sift3D->filter;
//	FILE* out_im;
//	char file_name[30];
//	sprintf_s(file_name, "D:\\debug\\keypoint_fix2.txt");
//	errno_t err;
//	err = fopen_s(&out_im, file_name, "w");
//	for (int i = 0; i < kp.size(); i++)
//	{
//		fprintf(out_im, "%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf,%.5lf\n",
//			kp[i].x, kp[i].y, kp[i].z, kp[i].scale, 
//			kp[i].Rotation[0], kp[i].Rotation[1], kp[i].Rotation[2], 
//			kp[i].Rotation[3], kp[i].Rotation[4], kp[i].Rotation[5],
//			kp[i].Rotation[6], kp[i].Rotation[7], kp[i].Rotation[8]);
//
//		//fprintf(out_im, "%.5lf,%.5lf,%.5lf,%.5lf,%d\n",
//		//	kp[i].x, kp[i].y, kp[i].z, kp[i].scale, kp[i].octave * 5 + kp[i].level - 1);// o = 3(start from 0), s = 1(start from 0) 
//	}
//	printf("size:%d\n", kp.size());
//	fclose(out_im);
//}
//
//
//void write_Desc(SIFT3D * sift3D)
//{
//	vector<Keypoint> kp = sift3D->filter;
//	FILE* out_desc;
//	char file_name[30];
//	sprintf_s(file_name, "D:\\debug\\description.txt");
//	errno_t err;
//	err = fopen_s(&out_desc, file_name, "w");
//	int num = 0;
//	for (int i = 0; i < kp.size(); i++)
//	{
//		Keypoint cur_kp = kp[i];
//		fprintf(out_desc, "%.5lf,%.5lf,%.5lf,%.5lf\n", cur_kp.rx, cur_kp.ry, cur_kp.rz, cur_kp.scale);
//		for (int j = 0; j < DESC_NUMEL; j++)
//		{
//			if (j % 12 == 0 && j > 0)
//			{
//				fprintf(out_desc, "\n");
//				//printf("\n");
//			}
//			fprintf(out_desc, "%.5lf,", *(cur_kp.desc + j));
//			//printf("%.5lf,", *(cur_kp.desc + j));
//			/*if (j < DESC_NUMEL - 1)
//				fprintf(out_desc, ",");*/
//			
//		}
//		fprintf(out_desc, "\n");
//	}
//
//	printf("size:%d\n", kp.size());
//	fclose(out_desc);
//}

void new_write_Desc(std::vector<Keypoint> &filter)
{
	FILE* out_desc;
	char file_name[100];
	sprintf_s(file_name, "D:\\debug3\\CPU\\description\\description_all.txt");
	errno_t err;
	err = fopen_s(&out_desc, file_name, "w");
	int num = 0;
	for (int i = 0; i < filter.size(); i++)
	{
		Keypoint cur_kp = filter[i];
		fprintf(out_desc, "%.5lf,%.5lf,%.5lf\n", cur_kp.rx, cur_kp.ry, cur_kp.rz);
		for (int j = 0; j < DESC_NUMEL; j++)
		{
			//if (j % 12 == 0 && j > 0)
			//{
			//	fprintf(out_desc, "\n");
			//	//printf("\n");
			//}
			fprintf(out_desc, "%.5lf,", *(cur_kp.desc + j));
			//printf("%.5lf,", *(cur_kp.desc + j));
			/*if (j < DESC_NUMEL - 1)
			fprintf(out_desc, ",");*/

		}
		fprintf(out_desc, "\n");
	}

	printf("size:%zd\n", filter.size());
	fclose(out_desc);
}


void write_sift_kp(std::vector<Cvec> &kp_coor, const char* file_name)
{
	FILE* out_desc;
	//char file_name[30];
	//sprintf_s(file_name, "D:\\debug2\\description.txt");
	errno_t err;
	err = fopen_s(&out_desc, file_name, "w");
	int num = 0;
	for (int i = 0; i < kp_coor.size(); i++)
	{
		Cvec cur_kp = kp_coor[i];
		fprintf(out_desc, "%.5lf,%.5lf,%.5lf\n", cur_kp.x, cur_kp.y, cur_kp.z);
	}

	printf("sift keypoint size:%zd\n", kp_coor.size());
	fclose(out_desc);
}

template<class Type>
Type stringToNum(string& str)
{
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}

Cvec explode(string& str, const char& ch)
{
	string next = "";
	vector<string> result;
	vector<float> num_res;
	for (string::const_iterator it = str.begin(); it != str.end(); it++)
	{
		if (*it == ch)
		{
			result.push_back(next);
			next.clear();
		}
		else
		{
			next += *it;
		}
	}
	if (!next.empty())
	{
		result.push_back(next);
	}

	for (int i = 0; i < result.size(); i++)
	{
		num_res.push_back(stringToNum<float>(result[i]));
	}

	Cvec cvec_res;
	cvec_res.x = num_res[0];
	cvec_res.y = num_res[1];
	cvec_res.z = num_res[2];

	return cvec_res;


}

void read_sift_kp(const char* file_name, std::vector<Cvec> &kp)
{
	ifstream input(file_name);

	string parsing;
	int i = 0;
	while (getline(input, parsing))
	{
		i++;
		//cout << "number idx: " << i << endl;
		Cvec tmp = explode(parsing, ',');
		kp.push_back(tmp);
	}
	cout << "File:"<< file_name <<", number of points:" << i << endl;
}


void new_write_level_Keypoint(std::vector<std::vector<Keypoint>>& level_extrema)
{
	printf("Writing Keypoint to txt now!\n");
	for (int i = 0; i < level_extrema.size(); i++)
	{
		FILE* out_im;
		char file_name[100];
		sprintf_s(file_name, "D:\\debug3\\CPU\\keypoint_level_%d.txt", i);
		errno_t err;
		err = fopen_s(&out_im, file_name, "w");

		//Write keypoint
		std::vector<Keypoint> level_kp = level_extrema[i];
		for (int i = 0; i < level_kp.size(); i++)
		{
			fprintf(out_im, "%.5lf,%.5lf,%.5lf,%.5lf\n",
				level_kp[i].x, level_kp[i].y, level_kp[i].z, level_kp[i].scale);

		}
		printf("level :%d  size is :%zd\n", i, level_kp.size());
		fclose(out_im);
	}
}


void write_Str_tensor(std::vector<Keypoint> extrema)
{
	char file_name[100];
	errno_t err;
	FILE* out_im;

	for (int i = 0; i < extrema.size(); i++)
	{
		Keypoint cur = extrema[i];
		int idx = cur.octave * NUM_KP_LEVELS + cur.level - 1;

		sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_str_tensor_level_%d.txt", idx);
		err = fopen_s(&out_im, file_name, "a");

		float* str_tensor_kp = cur.str_tensor;

		if (err == 0)
		{
			//write
			/*fprintf(out_im, "%.5lf, %.5lf, %.5lf, %.5lf\n", cur.x,
				cur.y, cur.z, cur.scale);

			fprintf(out_im, "str_tensor: %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf\n",
				str_tensor_kp[0], str_tensor_kp[1], str_tensor_kp[2],
				str_tensor_kp[3], str_tensor_kp[4], str_tensor_kp[5],
				str_tensor_kp[6], str_tensor_kp[7], str_tensor_kp[8]);*/

			fprintf(out_im, "%.5lf, %.5lf, %.5lf, %.5lf, %.5f, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf, %.5lf\n",
				cur.x, cur.y, cur.z, cur.scale,
				str_tensor_kp[0], str_tensor_kp[1], str_tensor_kp[2],
				str_tensor_kp[3], str_tensor_kp[4], str_tensor_kp[5],
				str_tensor_kp[6], str_tensor_kp[7], str_tensor_kp[8],
				cur.win.x, cur.win.y, cur.win.z);
		}

		fclose(out_im);
	}
}


void write_Eig(std::vector<Keypoint> extrema)
{
	char file_name[100];
	errno_t err;
	FILE* out_im;

	for (int i = 0; i < extrema.size(); i++)
	{
		Keypoint cur = extrema[i];
		int idx = cur.octave * NUM_KP_LEVELS + cur.level - 1;

		sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_eig_level_%d.txt", idx);
		err = fopen_s(&out_im, file_name, "a");

		float* str_eig = cur.eigvalue;
		float* str_eig_vec = cur.eigvector;

		if (err == 0)
		{
			//write
			fprintf(out_im, "%.6lf, %.6lf, %.6lf, %.6lf,\
					%.6lf, %.6lf, %.6lf,\
					%.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf\n",
				cur.x, cur.y, cur.z, cur.scale,
				str_eig[0], str_eig[1], str_eig[2],
				str_eig_vec[0], str_eig_vec[1], str_eig_vec[2],
				str_eig_vec[3], str_eig_vec[4], str_eig_vec[5],
				str_eig_vec[6], str_eig_vec[7], str_eig_vec[8]);
		}

		fclose(out_im);
	}
}


void write_Rot(std::vector<Keypoint> filter)
{
	char file_name[100];
	errno_t err;
	FILE* out_im;

	int num[15];
	memset(num, 0, sizeof(int) * 15);

	for (int i = 0; i < filter.size(); i++)
	{
		Keypoint cur = filter[i];
		int idx = cur.octave * NUM_KP_LEVELS + cur.level - 1;
		num[idx]++;

		sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_rot_level_%d.txt", idx);
		err = fopen_s(&out_im, file_name, "a");

		float* rot = cur.Rotation;

		if (err == 0)
		{
			//write
			fprintf(out_im, "%.6lf, %.6lf, %.6lf, %.6lf,\
					%.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf\n",
				cur.x, cur.y, cur.z, cur.scale,
				rot[0], rot[1], rot[2],
				rot[3], rot[4], rot[5],
				rot[6], rot[7], rot[8]);
		}

		fclose(out_im);
	}

	for (int i = 0; i < 15; i++)
	{
		printf("level: %d, num: %d\n", i, num[i]);
	}
}


void write_Desc(std::vector<Keypoint> filter)
{
	char file_name[100];
	errno_t err;
	FILE* out_im;

	int num[15];
	memset(num, 0, sizeof(int) * 15);

	for (int i = 0; i < filter.size(); i++)
	{
		Keypoint cur = filter[i];
		int idx = cur.octave * NUM_KP_LEVELS + cur.level - 1;
		num[idx]++;

		sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_desc_level_%d.txt", idx);
		err = fopen_s(&out_im, file_name, "a");

		float* desc = cur.desc;

		if (err == 0)
		{
			//write
			fprintf(out_im, "%.6lf, %.6lf, %.6lf, %.6lf\n",
				cur.x, cur.y, cur.z, cur.scale);

			for (int j = 0; j < 768; j++)
			{
				if (j % 12 == 0 && j > 0)
				{
					fprintf(out_im, "\n");
				}
				fprintf(out_im, "%.6lf,", desc[j]);

			}
			fprintf(out_im, "\n");
		}

		fclose(out_im);
	}

	for (int i = 0; i < 15; i++)
	{
		printf("level: %d, num: %d\n", i, num[i]);
	}
}


void write_debug_middle(float* loop_point, float* vrot, float* vbins, int* intersect_id, float* bary,
	float* dvbins, int* offset, float* desc_accum, int batchsize, Keypoint& kp, float* Rot)
{
	char file_name[100];

	//errno_t err;
	//FILE* out_im;
	sprintf_s(file_name, "D:\\debug3\\CPU\\cpu_desc_debug.txt");
	//err = fopen_s(&out_im, file_name, "a");

	std::ofstream out;
	out.setf(ios::fixed, ios::floatfield);
	out.precision(6);
	out.open(file_name, std::ios::app);

	out << kp.x << ", " << kp.y << ", " << kp.z << ", " << kp.scale << "\n";
	out << Rot[0] << ", " << Rot[1] << ", " << Rot[2] << ", " << Rot[3] 
		<< Rot[4] << ", " << Rot[5] << ", " << Rot[6] << ", " << Rot[7] 
		<< Rot[8] << "\n";

	for (int j = 0; j < batchsize; j++)
	{
		out << loop_point[j * 3 + 0] << ", " << loop_point[j * 3 + 1] << ", " << loop_point[j * 3 + 2]
			<< ", " << vrot[j * 3 + 0] << ", " << vrot[j * 3 + 1] << ", " << vrot[j * 3 + 2]
			<< ", " << vbins[j * 3 + 0] << ", " << vbins[j * 3 + 1] << ", " << vbins[j * 3 + 2]
			<< ", " << intersect_id[j]
			<< ", " << bary[j * 3 + 0] << ", " << bary[j * 3 + 1] << ", " << bary[j * 3 + 2]
			<< ", " << dvbins[j * 3 + 0] << ", " << dvbins[j * 3 + 1] << ", " << dvbins[j * 3 + 2] << "\n";
			//<< ", " << offset[j * 3 + 0] << ", " << offset[j * 3 + 1] << ", " << offset[j * 3 + 2]
			//<< ", " << desc_accum[j * 3 + 0] << ", " << desc_accum[j * 3 + 1] << ", " << desc_accum[j * 3 + 2] << "\n";
	
		/*for (int k = 0; k < 6; k++)
		{
			out <<  offset[j * 24 + k * 4 + 0] << ", " << desc_accum[j * 24 + k * 4 + 0] << ", " <<
					offset[j * 24 + k * 4 + 1] << ", " << desc_accum[j * 24 + k * 4 + 1] << ", " <<
					offset[j * 24 + k * 4 + 2] << ", " << desc_accum[j * 24 + k * 4 + 2] << ", " <<
					offset[j * 24 + k * 4 + 3] << ", " << desc_accum[j * 24 + k * 4 + 3] << "\n";

			if (k == 5)
				out << "\n";
		}*/
	
	}
	out << "\n\n\n\n\n";
	out.close();

	//if (err == 0)
	//{
		//write
	//printf("%.6lf, %.6lf, %.6lf, %.6lf\n", kp.x, kp.y, kp.z, kp.scale);
	//printf("%.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf\n",
	//	Rot[0], Rot[1], Rot[2], Rot[3], Rot[4], Rot[5], Rot[6], Rot[7], Rot[8]);
	//fprintf(out_im, "%.6lf, %.6lf, %.6lf, %.6lf\n", kp.x, kp.y, kp.z, kp.scale);
	//for (int j = 0; j < batchsize; j++)
	//{
	//	//fprintf(out_im, "%.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf\n",
	//	printf("%.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %d, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf, %.6lf\n",
	//		loop_point[j * 3 + 0], loop_point[j * 3 + 1], loop_point[j * 3 + 2],
	//		vrot[j * 3 + 0], vrot[j * 3 + 1], vrot[j * 3 + 2],
	//		vbins[j * 3 + 0], vbins[j * 3 + 1], vbins[j * 3 + 2],
	//		intersect_id[j],
	//		bary[j * 3 + 0], bary[j * 3 + 1], bary[j * 3 + 2],
	//		dvbins[j * 3 + 0], dvbins[j * 3 + 1], dvbins[j * 3 + 2],
	//		offset[j * 3 + 0], offset[j * 3 + 1], offset[j * 3 + 2],
	//		desc_accum[j * 3 + 0], desc_accum[j * 3 + 1], desc_accum[j * 3 + 2]);

	//}

	//printf("\n");
	//fprintf(out_im, "\n");


//}

//fclose(out_im);
}


void read_desc_all_level(char* file_name, std::vector<Keypoint>& kp)
{
	ifstream input(file_name);

	//input.open(file_name);
	std::stringstream ss;

	string parsing;
	int num_idx = 0;
	while (getline(input, parsing))
	{
		if (num_idx % 100 == 0)
			cout << num_idx << endl;

		ss.clear();
		ss.str(" ");
		ss << parsing;

		Keypoint new_kp;
		ss >> new_kp.x;
		ss >> new_kp.y;
		ss >> new_kp.z;

		for (int i = 0; i < 768; i++)
		{
			ss >> new_kp.desc[i];

		}

		kp.push_back(new_kp);

		num_idx++;
	}

	std::cout << "kp size :" << kp.size() << endl;
}
}
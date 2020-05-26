
#include "../../Include/Util/cTexImage.h"

#include <stdio.h>

#include "../../Include/Util/cMemManager.h"
#include "../../Include/cUtil.h"

TexImage::TexImage()
{
	_Data = nullptr;
	_numsize = 0;
	_nx = _ny = _nz = 0;
	_s = 0.0;
	_xs = _ys = _zs = 0;
	_ux = _uy = _uz = 0.0;
}

TexImage::TexImage(int width, int height, int depth)
{
	_Data = nullptr;
	_nx = width;
	_ny = height;
	_nz = depth;
	_numsize = width * height * depth * sizeof(float);

	_s = 0.0;
	_xs = 1;
	_ys = width;
	_zs = width * height;

	_ux = _uy = _uz = 1.0;

	//allocate the array memory
	//CMemManager<float>::hCreatePtr(_Data, _numsize);
}

void TexImage::ReSetImageSize(int width_new, int height_new, int depth_new)
{

	_Data = nullptr;
	_nx = width_new;
	_ny = height_new;
	_nz = depth_new;
	_numsize = width_new * height_new * depth_new * sizeof(float);

	_s = 0.0;
	_xs = 1;
	_ys = width_new;
	_zs = width_new * height_new;

	_ux = _uy = _uz = 1.0;
}

TexImage::~TexImage()
{
	if (_Data != nullptr)
	{
		CMemManager<float>::hDestroyPtr(_Data);
	}
}

void TexImage::MallocArrayMemory()
{
	if(_Data)
		CMemManager<float>::hDestroyPtr(_Data);

	CMemManager<float>::hCreatePtr(_Data, _nx * _ny * _nz);
}

void TexImage::SetImageSize(int width, int height, int depth)
{
	_nx = width;
	_ny = height;
	_nz = depth;
	_numsize = width * height * depth;

	_xs = 1;
	_ys = width;
	_zs = width * height;
}

void TexImage::SetImageScale(float scale)
{
	_s = scale;
}

void TexImage::SetImageUnit(float ux, float uy, float uz)
{
	_ux = ux;
	_uy = uy;
	_uz = uz;
}

void TexImage::SetImageDataPt(float* data)
{
	_Data = data;
}
//
//inline 
//void TexImage::SetImageDataWithIdx(float data, int coor_x, int coor_y, int coor_z)
//{
//	int idx = ELT2(coor_x, coor_y, coor_z, _xs, _ys, _zs);
//	_Data[idx] = data;
//}
//
//inline
//float TexImage::GetImageDataWithIdx(int coor_x, int coor_y, int coor_z)
//{
//	return _Data[ELT2(coor_x, coor_y, coor_z, _xs, _ys, _zs)];
//}
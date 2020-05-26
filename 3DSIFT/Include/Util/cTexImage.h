
#ifndef __CTEX_IMAGE_H__
#define __CTEX_IMAGE_H__

#define ELT2(x, y, z, x_stride, y_stride, z_stride) (long)((x) * x_stride + (y) * y_stride + (z) * z_stride)

//class for 3D volume data, using the texture memory
class TexImage
{
//protected:
public:
	//void* _cuData; //host or device data, binding data, for convience, set to device before, host data should be transferred first
	float* _Data; // how to access data in protected or private attribute
	//cudaArray* _cuData3D; // Used for destination Texture memory
	size_t _numsize;
	int _nx; // width
	int _ny; // height
	int _nz; // depth
	float _s; //scale
	size_t _xs;
	size_t _ys;
	size_t _zs;
	float _ux;
	float _uy;
	float _uz;


public:
	void SetImageSize(int width, int height, int depth);
	void ReSetImageSize(int width_new, int height_new, int depth_new);
	void SetImageScale(float scale);
	void SetImageUnit(float ux, float uy, float uz);
	void SetImageDataPt(float* data);
	void SetImageDataWithIdx(float data, int coor_x, int coor_y, int coor_z) {
		int idx = ELT2(coor_x, coor_y, coor_z, _xs, _ys, _zs);
		_Data[idx] = data;
	}
	void MallocArrayMemory();

public:
	TexImage();
	TexImage(int width, int height, int depth);
	virtual ~TexImage();

public:
	float GetScale() { return _s; }
	size_t GetXstride() { return _xs; }
	size_t GetYstride() { return _ys; }
	size_t GetZstride() { return _zs; }
	float GetUnitX() { return _ux; }
	float GetUnitY() { return _uy; }
	float GetUnitZ() { return _uz; }
	int GetDimX() { return _nx; }
	int GetDimY() { return _ny; }
	int GetDimZ() { return _nz; }
	float GetImageDataWithIdx(int coor_x, int coor_y, int coor_z) {
		return _Data[ELT2(coor_x, coor_y, coor_z, _xs, _ys, _zs)];
	}

};


#endif // !TEX_IAMGE_H

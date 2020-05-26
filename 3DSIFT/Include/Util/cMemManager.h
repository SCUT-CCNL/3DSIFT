
//!= Class responsible for the memory allocation and deallocation

#ifndef _CCMEMMANAGER_H_
#define _CCMEMMANAGER_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

template<class T>
class CMemManager
{
	//!= Host memory allocation methods.
public:
	// allocation
	static void hCreatePtr(T*& ptr, size_t size);
	static void hCreatePtr(T**& ptr, size_t row, size_t col);
	static void hCreatePtr(T***& ptr, size_t row, size_t col, size_t height);
	static void hCreatePtr(T****& ptr, size_t a, size_t b, size_t c, size_t d);
	static void hCreatePtr(T*****& ptr, size_t a, size_t b, size_t c, size_t d, size_t e);
	// delete
	static void hDestroyPtr(T*& ptr);
	static void hDestroyPtr(T**& ptr);
	static void hDestroyPtr(T***& ptr);
	static void hDestroyPtr(T****& ptr);
	static void hDestroyPtr(T*****& ptr);
};


template<class T>
void CMemManager<T>::hCreatePtr(T*& ptr, size_t size)
{
	ptr = (T*)malloc(size * sizeof(T));
	if (ptr != nullptr)
		memset(ptr, 0, size * sizeof(T));
	else
		fprintf(stderr, "Malloc memory fail\n");
}

template<class T>
void CMemManager<T>::hCreatePtr(T**& ptr, size_t row, size_t col)
{
	T * ptr1d = (T*)calloc(row*col, sizeof(T));
	ptr = (T**)malloc(row * sizeof(T*));

	for (size_t i = 0; i<row; i++)
	{
		ptr[i] = ptr1d + i*col;
	}
}

template<class T>
void CMemManager<T>::hCreatePtr(T***& ptr, size_t row, size_t col, size_t height)
{
	T *ptr1d = (T*)calloc(row*col*height, sizeof(T));
	T**ptr2d = (T**)malloc(row*col * sizeof(T*));
	ptr = (T***)malloc(row * sizeof(T**));

	for (size_t i = 0; i<row; i++)
	{
		for (size_t j = 0; j<col; j++)
		{
			ptr2d[i*col + j] = ptr1d + (i*col + j)*height;
		}
		ptr[i] = ptr2d + i*col;
	}
}

template<class T>
void CMemManager<T>::hCreatePtr(T****& ptr, size_t a, size_t b, size_t c, size_t d)
{
	T *ptr1d = (T*)calloc(a*b*c*d, sizeof(T));
	T**ptr2d = (T**)malloc(a*b*c * sizeof(T*));
	T***ptr3d = (T***)malloc(a*b * sizeof(T**));
	ptr = (T****)malloc(a * sizeof(T***));

	for (size_t i = 0; i<a; i++)
	{
		for (size_t j = 0; j<b; j++)
		{
			for (size_t k = 0; k<c; k++)
			{
				ptr2d[(i*b + j)*c + k] = ptr1d + ((i*b + j)*c + k)*d;
			}
			ptr3d[i*b + j] = ptr2d + (i*b + j)*c;
		}
		ptr[i] = ptr3d + i*b;
	}
}

template<class T>
void CMemManager<T>::hCreatePtr(T*****& ptr, size_t a, size_t b, size_t c, size_t d, size_t e)
{
	T *ptr1d = (T*)calloc(a*b*c*d*e, sizeof(T));
	T**ptr2d = (T**)malloc(a*b*c*d * sizeof(T*));
	T***ptr3d = (T***)malloc(a*b*c * sizeof(T**));
	T****ptr4d = (T****)malloc(a*b * sizeof(T***));
	ptr = (T*****)malloc(a * sizeof(T****));

	for (size_t i = 0; i<a; i++)
	{
		for (size_t j = 0; j<b; j++)
		{
			for (size_t k = 0; k<c; k++)
			{
				for (size_t l = 0; l < d; l++)
				{
					ptr2d[((i*b + j)*c + k)*d + l] = ptr1d + (((i*b + j)*c + k)*d + l)*e;
				}
				ptr3d[(i*b + j)*c + k] = ptr2d + ((i*b + j)*c + k)*d;
			}
			ptr4d[i*b + j] = ptr3d + (i*b + j)*c;
		}
		ptr[i] = ptr4d + i*b;
	}
}

//template<class T>
//void CMemManager<T>::hCreatePtr(T****& ptr, size_t a, size_t b, size_t c, size_t d, size_t e)
//{
//	T *ptr1d = (T*)calloc(a*b*c*d*e, sizeof(T));
//	T**ptr2d = (T**)malloc(a*b*c*d*sizeof(T*));
//	T***ptr3d = (T***)malloc(a*b*c*sizeof(T**));
//	T****ptr4d = (T****)malloc(a*b*sizeof(T***));
//	ptr = (T*****)malloc(a*sizeof(T****));
//
//	for (size_t i = 0; i<a; i++)
//	{
//		for (size_t j = 0; j<b; j++)
//		{
//			for (size_t k = 0; k<c; k++)
//			{
//				ptr2d[(i*b + j)*c + k] = ptr1d + ((i*b + j)*c + k)*d;
//			}
//			ptr3d[i*b + j] = ptr2d + (i*b + j)*c;
//		}
//		ptr[i] = ptr3d + i*b;
//	}
//}

template<class T>
void CMemManager<T>::hDestroyPtr(T*& ptr)
{
	if(ptr)
		free(ptr);
	ptr = nullptr;
}

template<class T>
void CMemManager<T>::hDestroyPtr(T**& ptr)
{
	free(ptr[0]);
	free(ptr);
	ptr = NULL;
}

template<class T>
void CMemManager<T>::hDestroyPtr(T***& ptr)
{
	free(ptr[0][0]);
	free(ptr[0]);
	free(ptr);
	ptr = NULL;
}

template<class T>
void CMemManager<T>::hDestroyPtr(T****& ptr)
{
	free(ptr[0][0][0]);
	free(ptr[0][0]);
	free(ptr[0]);
	free(ptr);
	ptr = NULL;
}

template<class T>
void CMemManager<T>::hDestroyPtr(T*****& ptr)
{
	free(ptr[0][0][0][0]);
	free(ptr[0][0][0]);
	free(ptr[0][0]);
	free(ptr[0]);
	free(ptr);
	ptr = NULL;
}

#endif // !_CMEMMANAGER_H_


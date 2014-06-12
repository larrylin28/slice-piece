#ifndef _CUDA_KERNEL_H
#define _CUDA_KERNEL_H

int freeSdf();
int getSdf(float* dists ,float* weights,float* n_x,float* n_y,float* n_z,int size);
int updateSdf(float* ptsx, float* ptsy, float* ptsz, float* norx, float* nory, float* norz,int size, float* mc,int width,int xres,int yres,float coeffx,float coeffy);
int initialSdf(float minx, float miny, float minz, int xlen,int ylen, int zlen, float devide, float epsi, float thre,float* dists, float* weights,float* nx,float* ny,float* nz,float* local);
int getAllDepths(float ix, float iy, float devide, int xl,int yl,float* mc, float* p_dists, float* p_nx, float* p_ny, float* p_nz);
int changeSdf(float minx, float miny, float minz, int xlen,int ylen, int zlen, float* local,float* olocalc);

struct kvertex
{
	float v[6];
	double dist;
	int tran;
	double upxyz[3];
	double upcount;
	double degree;
};

struct kmesh{
	int vertexes[3];
};

#endif
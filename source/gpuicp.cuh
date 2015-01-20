#ifndef _CUDA_GPUICP_H
#define _CUDA_GPUICP_H


int gpuICPIter(float* tr,int width, int height,int xres,int yres,float coeffx,float coeffy, float* ata, float* atb, float dthresh);
int gpuICPalloc(float* pto, float* ptn, float* noro, float* mc, float* tr,int width, int height);
void gpuICPfree();


#endif
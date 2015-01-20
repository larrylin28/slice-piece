#include "gpuicp.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>

#define ISZERO(x) x < 0.00000000001 && x > -0.00000000001

//dangerous code, be careful
#define GO_TO_ERROR(status, msg) if (status != cudaSuccess){fprintf(stderr, msg); goto Error;}


float *dev_pto = 0;
float *dev_ptn = 0;
float *dev_noro = 0;

float *dev_mc = 0;
float *dev_tr = 0;
float *dev_ata = 0;
float *dev_atb = 0;

__global__ void sum_linear(float* ata, float* atb, dim3 blockSize,int width, int height)
{
	for(int id = 0;id < width * height; id++)
	{
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				ata[i*6 + j] += ata[id*36 + i*6 + j];
			}
			atb[i] += atb[id*6 + i];
		}
	}

}

__global__ void sum_ICP(float* ata, float* atb, int iter, dim3 blockSize,int width, int height)
{
	int idx = threadIdx.x*blockSize.x+blockIdx.x;
	int idy = threadIdx.y*blockSize.y+blockIdx.y;

	int id = idy*width + idx;
	int idr = id + iter/2;

	if(id % iter == 0 && idr < width * height)
	{
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				ata[id*36 + i*6 + j] += ata[idr*36 + i*6 + j];
			}
			atb[id*6 + i] += atb[idr*6 + i];
		}
	}

}

//float* ptox, float* ptoy, float* ptoz, float* ptnx, float* ptny, float* ptnz, float* norox, float* noroy, float* noroz, float* nornx, float* norny, float* nornz, int size, float* mc, float* tr,int width,int xres,int yres,float coeffx,float coeffy
__global__ void cal_ICP(float* pto, float* ptn, float* noro, float* mc, float* tr, float* ata, float* atb,int nHalfXres,int nHalfYres, float fCoeffX, float fCoeffY,dim3 blockSize,dim3 threadSize,int width, int height, float dthresh)
{
	int idx = threadIdx.x*blockSize.x+blockIdx.x;
	int idy = threadIdx.y*blockSize.y+blockIdx.y;
	int tid = threadIdx.y*threadSize.x + threadIdx.x;

	__shared__ float shared_ata[2304];
	__shared__ float shared_atb[384];

	for(int i = 0; i < 6; i++)
	{
		for(int j = 0; j < 6; j++)
		{
			shared_ata[tid*36 + i*6 + j] = 0;
		}
		shared_atb[tid*6 + i] = 0;
	}

	float ptx = ptn[idy*width*4 + idx*4];
	float pty = ptn[idy*width*4 + idx*4 + 1];
	float ptz = ptn[idy*width*4 + idx*4 + 2];

	//if(ISZERO(ptx) && ISZERO(pty) && ISZERO(ptz)) return;

	float gptx = tr[0*4+0]*ptx+tr[0*4+1]*pty+tr[0*4+2]*ptz + tr[0*4+3]*1;
	float gpty = tr[1*4+0]*ptx+tr[1*4+1]*pty+tr[1*4+2]*ptz + tr[1*4+3]*1;
	float gptz = tr[2*4+0]*ptx+tr[2*4+1]*pty+tr[2*4+2]*ptz + tr[2*4+3]*1;

	float rx = mc[0*4+0]*gptx+mc[0*4+1]*gpty+mc[0*4+2]*gptz + mc[0*4+3]*1;
	float ry = mc[1*4+0]*gptx+mc[1*4+1]*gpty+mc[1*4+2]*gptz + mc[1*4+3]*1;
	float rz = mc[2*4+0]*gptx+mc[2*4+1]*gpty+mc[2*4+2]*gptz + mc[2*4+3]*1;

	int px = fCoeffX * (rx / 0.001f) / (rz / 0.001f) + nHalfXres;
    int py = nHalfYres - fCoeffY * (ry / 0.001f) / (rz / 0.001f);

	

	if(px >= 0 && px < 640 && py >= 0 && py < 480){
		float pox = pto[py*width*4 + px*4];
	    float poy = pto[py*width*4 + px*4 + 1];
		float poz = pto[py*width*4 + px*4 + 2];

		float nox = noro[py*width*4 + px*4];
	    float noy = noro[py*width*4 + px*4 + 1];
		float noz = noro[py*width*4 + px*4 + 2];

		float dis = (pox - gptx)*(pox - gptx) + (poy - gpty)*(poy - gpty) + (poz - gptz)*(poz - gptz);

        if(!((nox) != (nox)) && dis < dthresh)
		{
			float at[6];
			at[0] = gptz*noy - gpty*noz;
			at[1] = -gptz*nox + gptx*noz;
			at[2] = gpty*nox - gptx*noy;
			at[3] = nox;
			at[4] = noy;
			at[5] = noz;

			float b = nox*(pox - gptx) + noy*(poy - gpty) + noz*(poz - gptz);
			for(int i = 0; i < 6; i++)
			{
				for(int j = 0; j < 6; j++)
				{
					shared_ata[tid*36 + i*6 + j] = at[i] * at[j];
				}
				shared_atb[tid*6 + i] = at[i] * b;
			}
		}
	}

	
	__syncthreads();
	
	int iter = 8*8/2;
	while(iter != 0)
	{
		if(tid < iter)
		{
			for(int i = 0; i < 6; i++)
			{
				for(int j = 0; j < 6; j++)
				{
					shared_ata[tid*36 + i*6 + j] += shared_ata[(tid + iter)*36 + i*6 + j];
				}
				shared_atb[tid*6 + i] += shared_atb[(tid + iter)*6 + i];
			}
		}
		__syncthreads();
		iter /= 2;
	}
	
	__syncthreads();

	
	if(tid == 0)
	{
		int bid = blockIdx.y*blockSize.x + blockIdx.x;
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				ata[bid*36 + i*6 + j] = shared_ata[i*6 + j];
			}
			atb[bid*6 + i] = shared_atb[i];
		}
	}
	
}

int gpuICPalloc(float* pto, float* ptn, float* noro, float* mc, float* tr,int width, int height)
{
	dev_pto = 0;
	dev_ptn = 0;
    dev_noro = 0;

	dev_mc = 0;
	dev_tr = 0;
	dev_ata = 0;
	dev_atb = 0;

	cudaError_t cudaStatus = cudaSuccess;

	cudaStatus = cudaMalloc((void**)&dev_pto, width * height * 4 * sizeof(float));
	GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_ptn, width * height * 4 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_noro, width * height * 4 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_mc, 4*4 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_tr, 4*4 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_ata, (width * height / 8 / 8) * 36 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMalloc((void**)&dev_atb, (width * height / 8 / 8) * 6 * sizeof(float));
    GO_TO_ERROR(cudaStatus, "cudaMalloc failed!")

	cudaStatus = cudaMemcpy(dev_pto, pto, width * height * 4 * sizeof(float), cudaMemcpyHostToDevice);
	GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	cudaStatus = cudaMemcpy(dev_ptn, ptn, width * height * 4 * sizeof(float), cudaMemcpyHostToDevice);
    GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	cudaStatus = cudaMemcpy(dev_noro, noro, width * height * 4 * sizeof(float), cudaMemcpyHostToDevice);
    GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	cudaStatus = cudaMemcpy(dev_mc, mc, 4*4 * sizeof(float), cudaMemcpyHostToDevice);
    GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	cudaStatus = cudaMemcpy(dev_tr, tr, 4*4 * sizeof(float), cudaMemcpyHostToDevice);
    GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	return cudaStatus;

	Error:
    cudaFree(dev_pto);
	cudaFree(dev_ptn);
    cudaFree(dev_noro);
    cudaFree(dev_mc);
	cudaFree(dev_tr);
	cudaFree(dev_ata);
	cudaFree(dev_atb);

    return cudaStatus;
}

void gpuICPfree()
{
    cudaFree(dev_pto);
	cudaFree(dev_ptn);
    cudaFree(dev_noro);
    cudaFree(dev_mc);
	cudaFree(dev_tr);
	cudaFree(dev_ata);
	cudaFree(dev_atb);

	dev_pto = 0;
	dev_ptn = 0;
    dev_noro = 0;;
    dev_mc = 0;
	dev_tr = 0;
	dev_ata = 0;
	dev_atb = 0;

}

int gpuICPIter(float* tr,int width, int height,int xres,int yres,float coeffx,float coeffy, float* ata, float* atb, float dthresh)
{
	cudaError_t cudaStatus = cudaSuccess;

	cudaStatus = cudaMemcpy(dev_tr, tr, 4*4 * sizeof(float), cudaMemcpyHostToDevice);
    GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

    dim3 dimBlock(width/8,height/8);
	dim3 dimThread(8,8);
	cal_ICP<<<dimBlock, dimThread>>>(dev_pto, dev_ptn, dev_noro, dev_mc,dev_tr, dev_ata, dev_atb,xres,yres,coeffx,coeffy,dimBlock, dimThread, width, height, dthresh);

	cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }

	cudaStatus = cudaMemcpy(ata, dev_ata, (width * height / 8 / 8) * 36 * sizeof(float), cudaMemcpyDeviceToHost);
	GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	cudaStatus = cudaMemcpy(atb, dev_atb, (width * height / 8 / 8) * 6 * sizeof(float), cudaMemcpyDeviceToHost);
	GO_TO_ERROR(cudaStatus, "cudaMemcpy failed!")

	Error:

    return cudaStatus;

}

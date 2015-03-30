#include "icp.h"
#include "gpuicp.cuh"
#include <time.h>
#include <iostream>

int iter_max = 100;

Eigen::Matrix4f gpuICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh)
{
	clock_t start, finish, duration;

	start = clock(); 

	int width = old_cloud->width;
	int height = old_cloud->height;
	int size = width * height;
	float* pto = new float[size * 4];
	float* ptn = new float[size * 4];
	float* noro = new float[size * 4];
	for(int i = 0;i < size;i++)
	{
		pto[i*4] = old_cloud->points[i].x;
		pto[i*4 + 1] = old_cloud->points[i].y;
		pto[i*4 + 2] = old_cloud->points[i].z;
		ptn[i*4] = new_cloud->points[i].x;
		ptn[i*4 + 1] = new_cloud->points[i].y;
		ptn[i*4 + 2] = new_cloud->points[i].z;
		noro[i*4] = old_normals->points[i].normal_x;
		noro[i*4 + 1] = old_normals->points[i].normal_y;
		noro[i*4 + 2] = old_normals->points[i].normal_z;
	}

	float* mcf = new float[4*4];
	float* trf = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			mcf[i*4+j] = mc(i,j);
			trf[i*4+j] = tr(i,j);
		}
	}
	int status = 0;
	status = gpuICPalloc(pto, ptn, noro, mcf, trf, width, height);

	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"alloc:"<<duration<<std::endl;
	start = clock();
	
	int st = (width * height / 8 / 8);
	float* ata = new float[st * 36];
	float* atb = new float[st * 6];

	for(int k = 0; k < iter_max; k++)
	{
		
		int status = gpuICPIter(trf, width, height, xres, yres, coeffx, coeffy, ata, atb, dthresh);
		for(int iter = 1; iter < st; iter++)
		{
			for(int i = 0; i < 6; i++)
			{
				for(int j = 0; j < 6; j++)
				{
					ata[i*6 + j] += ata[iter*36 + i*6 + j];
				}
				atb[i] += atb[iter*6 + i];
			}
		}


		Eigen::MatrixXf A(6, 6);
		Eigen::VectorXf b(6);
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				A(i,j) = ata[i*6 + j];
			}
			b[i] = atb[i];
		}
		Eigen::VectorXf x = A.llt().solve(b);

		Eigen::Matrix4f tinc; 
	    tinc << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;
		tinc(0,1) = x[2];
		tinc(1,0) = -x[2];
		tinc(0,2) = -x[1];
		tinc(2,0) = x[1];
		tinc(1,2) = x[0];
		tinc(2,1) = -x[0];
		tinc(0,3) = x[3];
		tinc(1,3) = x[4];
		tinc(2,3) = x[5];

		tr = tinc*tr; 
		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				trf[i*4 + j] = tr(i,j);
			}
		}
	}
	delete[] ata;
	delete[] atb;

	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"iter:"<<duration<<std::endl;
	start = clock();
	
	gpuICPfree();
	delete[] pto;
	delete[] ptn;
	delete[] noro;
	delete[] mcf;
	delete[] trf;

	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"finish:"<<duration<<std::endl;
	return tr;
}

Eigen::Matrix4f cpuICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh)
{
	clock_t start, finish, duration;

	start = clock(); 

	int width = old_cloud->width;
	int height = old_cloud->height;
	int size = width * height;
	float* pto = new float[size * 3];
	float* ptn = new float[size * 3];
	float* noro = new float[size * 3];
	for(int i = 0;i < size;i++)
	{
		pto[i*3] = old_cloud->points[i].x;
		pto[i*3 + 1] = old_cloud->points[i].y;
		pto[i*3 + 2] = old_cloud->points[i].z;
		ptn[i*3] = new_cloud->points[i].x;
		ptn[i*3 + 1] = new_cloud->points[i].y;
		ptn[i*3 + 2] = new_cloud->points[i].z;
		noro[i*3] = old_normals->points[i].normal_x;
		noro[i*3 + 1] = old_normals->points[i].normal_y;
		noro[i*3 + 2] = old_normals->points[i].normal_z;
	}

	float* mcf = new float[4*4];
	float* trf = new float[4*4];

	
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			mcf[i*4+j] = mc(i,j);
			trf[i*4+j] = tr(i,j);
		}
	}

	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"alloc:"<<duration<<std::endl;
	start = clock();
	
	
	float* ata = new float[36];
	float* atb = new float[6];
	for(int k = 0; k < iter_max; k++)
	{

		for(int i = 0; i < 36; i++) ata[i] = 0;
		for(int i = 0; i < 6; i++) atb[i] = 0;
		
		//int status = cpuICPIter(trf, width, height, xres, yres, coeffx, coeffy, ata, atb);
		for(int idx = 0; idx < width; idx++)
		{
			for(int idy = 0; idy < height; idy++)
			{
				float ptx = ptn[idy*width*3 + idx*3];
				float pty = ptn[idy*width*3 + idx*3 + 1];
				float ptz = ptn[idy*width*3 + idx*3 + 2];

				float gptx = trf[0*4+0]*ptx+trf[0*4+1]*pty+trf[0*4+2]*ptz + trf[0*4+3]*1;
				float gpty = trf[1*4+0]*ptx+trf[1*4+1]*pty+trf[1*4+2]*ptz + trf[1*4+3]*1;
				float gptz = trf[2*4+0]*ptx+trf[2*4+1]*pty+trf[2*4+2]*ptz + trf[2*4+3]*1;

				float rx = mcf[0*4+0]*gptx+mcf[0*4+1]*gpty+mcf[0*4+2]*gptz + mcf[0*4+3]*1;
				float ry = mcf[1*4+0]*gptx+mcf[1*4+1]*gpty+mcf[1*4+2]*gptz + mcf[1*4+3]*1;
				float rz = mcf[2*4+0]*gptx+mcf[2*4+1]*gpty+mcf[2*4+2]*gptz + mcf[2*4+3]*1;

				int px = coeffx * (rx / 0.001f) / (rz / 0.001f) + xres;
				int py = yres - coeffy * (ry / 0.001f) / (rz / 0.001f);

				if(px >= 0 && px < 640 && py >= 0 && py < 480){
					float pox = pto[py*width*3 + px*3];
					float poy = pto[py*width*3 + px*3 + 1];
					float poz = pto[py*width*3 + px*3 + 2];

					float nox = noro[py*width*3 + px*3];
					float noy = noro[py*width*3 + px*3 + 1];
					float noz = noro[py*width*3 + px*3 + 2];

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
								ata[i*6 + j] += at[i] * at[j];
							}
							atb[i] += at[i] * b;
						}
					}
				}
			}
		}

		Eigen::MatrixXf A(6, 6);
		Eigen::VectorXf b(6);
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				A(i,j) = ata[i*6 + j];
			}
			b[i] = atb[i];
		}
		Eigen::VectorXf x = A.llt().solve(b);

		Eigen::Matrix4f tinc; 
	    tinc << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;
		tinc(0,1) = x[2];
		tinc(1,0) = -x[2];
		tinc(0,2) = -x[1];
		tinc(2,0) = x[1];
		tinc(1,2) = x[0];
		tinc(2,1) = -x[0];
		tinc(0,3) = x[3];
		tinc(1,3) = x[4];
		tinc(2,3) = x[5];

		tr = tinc*tr; 
		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				trf[i*4 + j] = tr(i,j);
			}
		}
	}
	delete[] ata;
	delete[] atb;
	
	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"iter:"<<duration<<std::endl;

	delete[] pto;
	delete[] ptn;
	delete[] noro;
	delete[] mcf;
	delete[] trf;

	return tr;
}


Eigen::Matrix4f cpuExtraICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh,
	                        int extraTag, int* oldTags, int* newTags, Eigen::Matrix4f extraTr)
{
	clock_t start, finish, duration;

	start = clock(); 

	int width = old_cloud->width;
	int height = old_cloud->height;
	int size = width * height;
	float* pto = new float[size * 3];
	float* ptn = new float[size * 3];
	float* noro = new float[size * 3];
	for(int i = 0;i < size;i++)
	{
		pto[i*3] = old_cloud->points[i].x;
		pto[i*3 + 1] = old_cloud->points[i].y;
		pto[i*3 + 2] = old_cloud->points[i].z;
		ptn[i*3] = new_cloud->points[i].x;
		ptn[i*3 + 1] = new_cloud->points[i].y;
		ptn[i*3 + 2] = new_cloud->points[i].z;
		noro[i*3] = old_normals->points[i].normal_x;
		noro[i*3 + 1] = old_normals->points[i].normal_y;
		noro[i*3 + 2] = old_normals->points[i].normal_z;
	}

	float* mcf = new float[4*4];
	float* trf = new float[4*4];

	float* etrf = NULL;
	float* emcf = NULL;
	Eigen::Matrix4f extraMc = extraTr.inverse();
	etrf = new float[4*4];
	emcf = new float[4*4];

	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			mcf[i*4+j] = mc(i,j);
			trf[i*4+j] = tr(i,j);
			etrf[i*4+j] = extraTr(i,j);
			emcf[i*4+j] = extraMc(i,j);
		}
	}

	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"alloc:"<<duration<<std::endl;
	start = clock();
	
	
	float* ata = new float[36];
	float* atb = new float[6];
	for(int k = 0; k < iter_max; k++)
	{

		for(int i = 0; i < 36; i++) ata[i] = 0;
		for(int i = 0; i < 6; i++) atb[i] = 0;
		
		//int status = cpuICPIter(trf, width, height, xres, yres, coeffx, coeffy, ata, atb);
		for(int idx = 0; idx < width; idx++)
		{
			for(int idy = 0; idy < height; idy++)
			{
				int newid = idy*width + idx;
				if(newTags[newid] != extraTag) continue;

				float ptx = ptn[idy*width*3 + idx*3];
				float pty = ptn[idy*width*3 + idx*3 + 1];
				float ptz = ptn[idy*width*3 + idx*3 + 2];

				float gptx = trf[0*4+0]*ptx+trf[0*4+1]*pty+trf[0*4+2]*ptz + trf[0*4+3]*1;
				float gpty = trf[1*4+0]*ptx+trf[1*4+1]*pty+trf[1*4+2]*ptz + trf[1*4+3]*1;
				float gptz = trf[2*4+0]*ptx+trf[2*4+1]*pty+trf[2*4+2]*ptz + trf[2*4+3]*1;

				float rx = mcf[0*4+0]*gptx+mcf[0*4+1]*gpty+mcf[0*4+2]*gptz + mcf[0*4+3]*1;
				float ry = mcf[1*4+0]*gptx+mcf[1*4+1]*gpty+mcf[1*4+2]*gptz + mcf[1*4+3]*1;
				float rz = mcf[2*4+0]*gptx+mcf[2*4+1]*gpty+mcf[2*4+2]*gptz + mcf[2*4+3]*1;

				int px = coeffx * (rx / 0.001f) / (rz / 0.001f) + xres;
				int py = yres - coeffy * (ry / 0.001f) / (rz / 0.001f);

				if(px >= 0 && px < 640 && py >= 0 && py < 480){
					float pox = pto[py*width*3 + px*3];
					float poy = pto[py*width*3 + px*3 + 1];
					float poz = pto[py*width*3 + px*3 + 2];

					float nox = noro[py*width*3 + px*3];
					float noy = noro[py*width*3 + px*3 + 1];
					float noz = noro[py*width*3 + px*3 + 2];

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
								ata[i*6 + j] += at[i] * at[j];
							}
							atb[i] += at[i] * b;
						}
					}
				}
			}
		}

		Eigen::MatrixXf A(6, 6);
		Eigen::VectorXf b(6);
		for(int i = 0; i < 6; i++)
		{
			for(int j = 0; j < 6; j++)
			{
				A(i,j) = ata[i*6 + j];
			}
			b[i] = atb[i];
		}
		Eigen::VectorXf x = A.llt().solve(b);

		Eigen::Matrix4f tinc; 
	    tinc << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;
		tinc(0,1) = x[2];
		tinc(1,0) = -x[2];
		tinc(0,2) = -x[1];
		tinc(2,0) = x[1];
		tinc(1,2) = x[0];
		tinc(2,1) = -x[0];
		tinc(0,3) = x[3];
		tinc(1,3) = x[4];
		tinc(2,3) = x[5];

		tr = tinc*tr; 
		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				trf[i*4 + j] = tr(i,j);
			}
		}
	}
	delete[] ata;
	delete[] atb;
	
	finish = clock(); 
	duration = (double)(finish - start) / CLOCKS_PER_SEC; 
	std::cout<<"iter:"<<duration<<std::endl;

	delete[] pto;
	delete[] ptn;
	delete[] noro;
	delete[] mcf;
	delete[] trf;

	return tr;
}

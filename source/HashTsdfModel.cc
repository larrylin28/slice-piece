#include "HashTsdfModel.h"
#include "hashSdf.cuh"

namespace Rgbd
{	

    HashTsdfModel::HashTsdfModel()
	{
		createHashSdf();
	}

	HashTsdfModel::~HashTsdfModel()
	{
	}

	void HashTsdfModel::dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader)
	{
	}

	void HashTsdfModel::rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz)
	{
	}

	void HashTsdfModel::freeData()
	{
	}
}
#include "SingleTsdfModel.h"

namespace Rgbd
{	
        SingleTsdfModel::SingleTsdfModel()
		{
		}

		SingleTsdfModel::~SingleTsdfModel()
		{
		}

		void SingleTsdfModel::dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY)

		{
		}

		void SingleTsdfModel::rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz)
		{
		}

		void SingleTsdfModel::freeData()
		{
		}
}
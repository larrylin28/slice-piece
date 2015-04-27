#include "Segmentation.h"
#include "SingleTsdfModel.h"

namespace Rgbd
{	
	
        SingleTsdfModel::SingleTsdfModel(Block& block)
		{
		}

		SingleTsdfModel::~SingleTsdfModel()
		{
		}

		void SingleTsdfModel::dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader)

		{
		}

		void SingleTsdfModel::dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader,
						int specTag, std::vector<int>& tags)
		{
		}

		void SingleTsdfModel::rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz)
		{
		}

		void SingleTsdfModel::getSample(float sdevide, float sminX, float sminY, float sminZ, int sXlen, int sYlen, int sZlen, float* absrate)
		{
		}

		void SingleTsdfModel::freeData()
		{
		}
}
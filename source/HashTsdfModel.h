#ifndef HASH_TSDF_MODEL_H
#define HASH_TSDF_MODEL_H

#include "SdfModel.h"

namespace Rgbd
{	

	class HashTsdfModel : public SdfModel
	{
	public:
		HashTsdfModel();
		~HashTsdfModel();

	public:
      
		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader);

		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz);

		virtual void freeData();

	private:
		
	};
}


#endif // HASH_TSDF_MODEL_H

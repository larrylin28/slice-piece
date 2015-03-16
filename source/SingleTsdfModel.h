#ifndef SINGLETSDF_MODEL_H
#define SINGLETSDF_MODEL_H

#include "SdfModel.h"
#include "Plane3D.h"

namespace Rgbd
{	
	

	class SingleTsdfModel : public SdfModel
	{
	public:
        SingleTsdfModel();
		~SingleTsdfModel();

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY);


		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz);

		virtual void freeData();

	};
}
#endif // SINGLETSDF_MODEL_H

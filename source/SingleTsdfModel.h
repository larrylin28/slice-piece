#ifndef SINGLETSDF_MODEL_H
#define SINGLETSDF_MODEL_H

#include "SdfModel.h"
#include "Plane3D.h"
#include "Segmentation.h"

namespace Rgbd
{	

	class SingleTsdfModel : public SdfModel
	{
	public:
        SingleTsdfModel(Block& block);
		~SingleTsdfModel();

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader);

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader,
						int specTag, std::vector<int>& tags);


		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz);

		virtual void getSample(float sdevide, float sminX, float sminY, float sminZ, int sXlen, int sYlen, int sZlen, float* absrate);

		virtual void freeData();

	};
}
#endif // SINGLETSDF_MODEL_H

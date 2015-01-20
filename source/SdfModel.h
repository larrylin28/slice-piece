#ifndef SDF_MODEL_H
#define SDF_MODEL_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace Rgbd
{
	typedef pcl::PointXYZRGB PointT;
    typedef pcl::PointCloud<PointT> PointCloud;
	typedef pcl::PointCloud<pcl::Normal> PointNormal;
	

	class SdfModel
	{
		public:
		SdfModel()
		{
		}

		virtual ~SdfModel()
		{
		};

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY) = 0;

		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz) = 0;

		virtual void freeData() = 0;

	};
}

#endif // SDF_MODEL_H

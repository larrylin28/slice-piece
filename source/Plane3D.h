#ifndef PLANE3D_H
#define PLANE3D_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace Rgbd
{
	typedef pcl::PointXYZRGB PointT;
    typedef pcl::PointCloud<PointT> PointCloud;
	typedef pcl::PointCloud<pcl::Normal> PointNormal;

	

	struct Plane3D{
		double a;
		double b;
		double c;
		double d;
		int inliers;

		int id; //∂ÓÕ‚ ≈≈–Ú”√
	};

	struct coffCal
	{
		double Ex;
		double Ey;
		double Ez;
		double Ex2;
		double Ey2;
		double Ez2;
		double Exy;
		double Eyz;
		double Exz;

		coffCal():Ex(0),Ey(0),Ez(0),Ex2(0),Ey2(0),Ez2(0),Exy(0),Eyz(0),Exz(0)
		{
		}

		void addPoint(PointT pp)
		{
			Ex += pp.x;
			Ey += pp.y;
			Ez += pp.z;
			Ex2 += pp.x*pp.x;
			Ey2 += pp.y*pp.y;
			Ez2 += pp.z*pp.z;
			Exy += pp.x*pp.y;
			Eyz += pp.y*pp.z;
			Exz += pp.x*pp.z;
		}

		Plane3D getCoff();
	};

	double p_distance(Plane3D& a, Plane3D& b);
	double p_distance(Plane3D& p, double a, double b, double c, double d);

	double p_angle(Plane3D& a, Plane3D& b);
	double p_angle(Plane3D& p, double a, double b, double c, double d);

	Plane3D coff_cloud(pcl::PointCloud<pcl::Normal>::Ptr normals, PointCloud::Ptr cloud);

	void dbplane_segment(PointCloud::Ptr cloud, PointNormal::Ptr normals,
		                 int minInlier, double sthre, double athre, double dthre,
						 std::vector<Plane3D>& coffs, std::vector<int>& seg);
}

#endif // PLANE3D_H

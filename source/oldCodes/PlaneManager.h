#ifndef PLANE_MANAGER_H_
#define PLANE_MANAGER_H_

#include <vector>
#include <pcl/common/eigen.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include "KeyPointManager.h"

#endif

using namespace std;


struct Plane3D{
	double a;
	double b;
	double c;
	double d;
	int inliers;
};

struct plane_info{
	double a;
	double b;
	double c;
	double d;
	int inliers;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr bounds;
	Point3D* boundPoints;
	bool skeletonGenerated;
	int clusteredID;
};

struct plane_match{
	plane_info source;
	plane_info target;
};

class PlaneManager{
   public:
     vector<plane_info> cloud_planes;
	 vector<plane_info>* clutered_planes;
	 double treshold;
	 int minInliers;
	 Eigen::Vector3f x_coord; 
	 Eigen::Vector3f y_coord; 
	 Eigen::Vector3f z_coord; 
	 int coordbuildstep;

	 PlaneManager();
	 int addPlane(double a,double b,double c,double d,int inliers,pcl::PointCloud<pcl::PointXYZRGB>::VectorType contour);
	 int addPlane(double a,double b,double c,double d,int inliers);
	 void clusterPlane();
	 void generateSkeleton();
	 void printPlaneInfo();
	 void PlaneManager::getCoordinary();
	 Point3D getCrossPoint(Plane3D p1,Plane3D p2,Plane3D p3);
	 double det(double a1,double b1,double c1,double a2,double b2,double c2,double a3,double b3,double c3);
	 
};
#ifndef KEYPOINT_MANAGER_H_
#define KEYPOINT_MANAGER_H_

#include <vector>
#include <pcl/common/eigen.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

#endif

using namespace std;

struct KeyPointID{
	int frameID;
	int keyID;
};

struct Point3D{
	double x;
	double y;
	double z;
};

class virtualKeypoint{
public:
	vector<Point3D> keypoints;
	Point3D currentMean;

	void addPoint(Point3D p);
	void addPoint(double x,double y,double z);
};

class KeyPointManager{
public:
	int virtualKeypointsCount;
	vector<virtualKeypoint> virtualKeypoints;

	void addvirtualPoint(double x,double y,double z);
};
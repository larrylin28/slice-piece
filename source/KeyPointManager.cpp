#include "KeyPointManager.h"
#include <iostream>

using namespace std;

void virtualKeypoint::addPoint(double x,double y,double z){
	Point3D p = *(new Point3D());
	p.x = x;
	p.y = y;
	p.z = z;
	addPoint(p);
}

void virtualKeypoint::addPoint(Point3D p){
	int size = keypoints.size();
	keypoints.push_back(p);
	if(size <= 0){
		currentMean = p;
	}else{
		currentMean.x = (currentMean.x*size + p.x)/(size+1);
		currentMean.y = (currentMean.y*size + p.y)/(size+1);
		currentMean.z = (currentMean.z*size + p.z)/(size+1);
	}
}

void KeyPointManager::addvirtualPoint(double x,double y,double z){
	virtualKeypoint vk = *(new virtualKeypoint());
	vk.addPoint(x,y,z);
	virtualKeypoints.push_back(vk);
	virtualKeypointsCount++;
}
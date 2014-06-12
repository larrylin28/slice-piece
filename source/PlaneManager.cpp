#include "PlaneManager.h"
#include <iostream>

using namespace std;

//自定义比较函数  
bool planeInfoCompare( const plane_info &v1, const plane_info &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{  
	return v1.inliers > v2.inliers;  
}

PlaneManager::PlaneManager(){
	coordbuildstep = 0;
	treshold = 1;
	minInliers = 100000;
	clutered_planes = NULL;
}

void PlaneManager::clusterPlane(){
	//确定K值和初始点
	double tr = 0.5;
	if(cloud_planes.size() == 0) return;
	int K = 1;
	clutered_planes = new vector<plane_info>();
    plane_info* m = new plane_info();
	m->a = cloud_planes.at(0).a;
	m->b = cloud_planes.at(0).b;
	m->c = cloud_planes.at(0).c;
	m->d = cloud_planes.at(0).d;
	m->inliers = 0;
	m->bounds = (pcl::PointCloud<pcl::PointXYZRGB>::Ptr)(new pcl::PointCloud<pcl::PointXYZRGB>);
	clutered_planes->push_back(*m);
	for(int i = 1;i < cloud_planes.size();i++){
		bool addK = true;
		for(int j = 0;j < K;j++){
			//double dis = pcl::euclideanDistance(p, mean->at(j));
			double dis = sqrt((clutered_planes->at(j).a - cloud_planes.at(i).a)*(clutered_planes->at(j).a - cloud_planes.at(i).a) + (clutered_planes->at(j).b - cloud_planes.at(i).b)*(clutered_planes->at(j).b - cloud_planes.at(i).b) + (clutered_planes->at(j).c - cloud_planes.at(i).c)*(clutered_planes->at(j).c - cloud_planes.at(i).c) + (clutered_planes->at(j).d - cloud_planes.at(i).d)*(clutered_planes->at(j).d - cloud_planes.at(i).d));
			if(dis < tr){
				addK = false;
				break;
			}
		}
		if(addK){
			K++;
			plane_info* m = new plane_info();
			m->a = cloud_planes.at(i).a;
			m->b = cloud_planes.at(i).b;
			m->c = cloud_planes.at(i).c;
			m->d = cloud_planes.at(i).d;
			m->inliers = 0;
			m->bounds = (pcl::PointCloud<pcl::PointXYZRGB>::Ptr)(new pcl::PointCloud<pcl::PointXYZRGB>);
			clutered_planes->push_back(*m);
		}
	}
	//进行聚类
	for(int loop = 0;loop < 1000;loop++){
		int* count = new int[K];
		double* amean = new double[K];
		double* bmean = new double[K];
		double* cmean = new double[K];
		double* dmean = new double[K];
		for(int i = 0;i < K;i++){
			count[i] = 0;
			amean[i] = 0.0;
			bmean[i] = 0.0;
			cmean[i] = 0.0;
			dmean[i] = 0.0;
		}
		for(int i = 0;i < cloud_planes.size();i++){
			for(int j = 0;j < K;j++){
				double dis = sqrt((clutered_planes->at(j).a - cloud_planes.at(i).a)*(clutered_planes->at(j).a - cloud_planes.at(i).a) + (clutered_planes->at(j).b - cloud_planes.at(i).b)*(clutered_planes->at(j).b - cloud_planes.at(i).b) + (clutered_planes->at(j).c - cloud_planes.at(i).c)*(clutered_planes->at(j).c - cloud_planes.at(i).c) + (clutered_planes->at(j).d - cloud_planes.at(i).d)*(clutered_planes->at(j).d - cloud_planes.at(i).d));
				if(dis < tr){
					count[j] += cloud_planes.at(i).inliers;
					amean[j] += cloud_planes.at(i).a * cloud_planes.at(i).inliers;
					bmean[j] += cloud_planes.at(i).b * cloud_planes.at(i).inliers;
					cmean[j] += cloud_planes.at(i).c * cloud_planes.at(i).inliers;
					dmean[j] += cloud_planes.at(i).d * cloud_planes.at(i).inliers;
					break;
				}
			}
		}
		bool iterend = true;
		double iterbound = 0.01;
		for(int i = 0;i < K;i++){
			amean[i] /= count[i];
			bmean[i] /= count[i];
			cmean[i] /= count[i];
			dmean[i] /= count[i];
			if(amean[i] - clutered_planes->at(i).a > iterbound || bmean[i] - clutered_planes->at(i).b > iterbound && cmean[i] - clutered_planes->at(i).c > iterbound && dmean[i] - clutered_planes->at(i).d > iterbound){
				iterend = false;
			}
			clutered_planes->at(i).a = amean[i];
			clutered_planes->at(i).b = bmean[i];
			clutered_planes->at(i).c = cmean[i];
			clutered_planes->at(i).d = dmean[i];
		}
		if(iterend){ break;}
	}
	for(int i = 0;i < cloud_planes.size();i++){
		for(int j = 0;j < K;j++){
			double dis = sqrt((clutered_planes->at(j).a - cloud_planes.at(i).a)*(clutered_planes->at(j).a - cloud_planes.at(i).a) + (clutered_planes->at(j).b - cloud_planes.at(i).b)*(clutered_planes->at(j).b - cloud_planes.at(i).b) + (clutered_planes->at(j).c - cloud_planes.at(i).c)*(clutered_planes->at(j).c - cloud_planes.at(i).c) + (clutered_planes->at(j).d - cloud_planes.at(i).d)*(clutered_planes->at(j).d - cloud_planes.at(i).d));
			if(dis < tr){
				clutered_planes->at(j).inliers += cloud_planes.at(i).inliers;
				cloud_planes.at(i).clusteredID = j;
				for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = cloud_planes.at(i).bounds->points.begin();it != cloud_planes.at(i).bounds->points.end();it++){
			        clutered_planes->at(j).bounds->points.push_back(*it);
		        } 
				break;
			}
		}
	}
}
int PlaneManager::addPlane(double a,double b,double c,double d,int inliers,pcl::PointCloud<pcl::PointXYZRGB>::VectorType contour){
	plane_info plane = *(new plane_info());
	plane.a = a;
	plane.b = b;
	plane.c = c;
	plane.d = d;
	plane.inliers = inliers;
	plane.bounds = (pcl::PointCloud<pcl::PointXYZRGB>::Ptr)(new pcl::PointCloud<pcl::PointXYZRGB>);
	plane.bounds->points = contour;
	plane.clusteredID = -1;
	cloud_planes.push_back(plane);
	return cloud_planes.size()-1;
}
int PlaneManager::addPlane(double a,double b,double c,double d,int inliers){
	plane_info plane = *(new plane_info());
	plane.a = a;
	plane.b = b;
	plane.c = c;
	plane.d = d;
	plane.inliers = inliers;
	plane.bounds = (pcl::PointCloud<pcl::PointXYZRGB>::Ptr)(new pcl::PointCloud<pcl::PointXYZRGB>);
	plane.clusteredID = -1;
	cloud_planes.push_back(plane);
	return cloud_planes.size()-1;
}
/*
void PlaneManager::addPlane(double a,double b,double c,double d,int inliers,pcl::PointCloud<pcl::PointXYZRGB>::VectorType contour){
	bool needAdd = true;
	for(int i = 0;i < cloud_planes.size();i++) {
		plane_info p = cloud_planes.at(i);
	   double dis1 = sqrt((p.a - a)*(p.a - a) + (p.b - b)*(p.b - b) + (p.c - c)*(p.c - c) + (p.d - d)*(p.d - d));
	   double dis2 = sqrt((p.a + a)*(p.a + a) + (p.b + b)*(p.b + b) + (p.c + c)*(p.c + c) + (p.d + d)*(p.d + d));
	   if(dis1 < treshold || dis2 < treshold){
		   cloud_planes.at(i).a = (p.a*p.inliers + a*inliers)/(p.inliers + inliers);
		   cloud_planes.at(i).b = (p.b*p.inliers + b*inliers)/(p.inliers + inliers);
		   cloud_planes.at(i).c = (p.c*p.inliers + c*inliers)/(p.inliers + inliers);
		   cloud_planes.at(i).d = (p.d*p.inliers + d*inliers)/(p.inliers + inliers);
		   cloud_planes.at(i).inliers = p.inliers + inliers;
		   for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = contour.begin();it != contour.end();it++){
			   cloud_planes.at(i).bounds->points.push_back(*it);
		   }
		   needAdd = false;
		   break;
	   }
    }
	if(needAdd){
		   plane_info plane = *(new plane_info());
		   plane.a = a;
		   plane.b = b;
		   plane.c = c;
		   plane.d = d;
		   plane.inliers = inliers;
		   plane.bounds = (pcl::PointCloud<pcl::PointXYZRGB>::Ptr)(new pcl::PointCloud<pcl::PointXYZRGB>);
		   plane.bounds->points = contour;
		   cloud_planes.push_back(plane);
	}
}
*/
void PlaneManager::printPlaneInfo()
{
	for(vector<plane_info>::const_iterator iter = clutered_planes->begin(); iter != clutered_planes->end(); ++iter) {
		 plane_info p = *iter;
		 cout<<"("<<p.a<<","<<p.b<<","<<p.c<<","<<p.d<<") inliers:"<<p.inliers;
	}
	cout<<"x_coord:("<<x_coord[0]<<","<<x_coord[1]<<","<<x_coord[2]<<")";
	cout<<"y_coord:("<<y_coord[0]<<","<<y_coord[1]<<","<<y_coord[2]<<")";
	cout<<"z_coord:("<<z_coord[0]<<","<<z_coord[1]<<","<<z_coord[2]<<")";
}
void PlaneManager::getCoordinary()
{
	if(clutered_planes == NULL) clusterPlane();
	std::sort(clutered_planes->begin(),clutered_planes->end(),planeInfoCompare);  
	for(vector<plane_info>::const_iterator iter = clutered_planes->begin(); iter != clutered_planes->end(); ++iter) {
		 plane_info p = *iter;
		 if(coordbuildstep == 0){
			 x_coord[0] = p.a;
			 x_coord[1] = p.b;
			 x_coord[2] = p.c;
			 double norm = sqrt(x_coord[0]*x_coord[0] + x_coord[1]*x_coord[1] + x_coord[2]*x_coord[2]);
			 x_coord[0] = x_coord[0]/norm;
			 x_coord[1] = x_coord[1]/norm;
			 x_coord[2] = x_coord[2]/norm;
			 coordbuildstep = 1;
		 }else if(coordbuildstep == 1){
			 double angle =  (x_coord[0] * p.a + x_coord[1] * p.b + x_coord[2] * p.c)/sqrt(x_coord[0] * x_coord[0] + x_coord[1] * x_coord[1] + x_coord[2] * x_coord[2])/sqrt(p.a * p.a + p.b * p.b + p.c * p.c);
		     double delt = (x_coord[0] * p.a + x_coord[1] * p.b + x_coord[2] * p.c)/sqrt(x_coord[0] * x_coord[0] + x_coord[1] * x_coord[1] + x_coord[2] * x_coord[2]);
			 if(angle < 0.1 && angle > -0.1){
				 y_coord[0] = p.a;
				 y_coord[1] = p.b;
				 y_coord[2] = p.c;
				 double norm = sqrt(y_coord[0]*y_coord[0] + y_coord[1]*y_coord[1] + y_coord[2]*y_coord[2]);
				 y_coord[0] = y_coord[0]/norm;
				 y_coord[1] = y_coord[1]/norm;
				 y_coord[2] = y_coord[2]/norm;
				 coordbuildstep = 2;
				 break;
			 }
		 }
	}
	//求z方向
	//叉乘x方向和y方向
	z_coord[0] = x_coord[1]*y_coord[2] - x_coord[2]*y_coord[1];
    z_coord[1] = x_coord[2]*y_coord[0] - x_coord[0]*y_coord[2];
	z_coord[2] = x_coord[0]*y_coord[1] - x_coord[1]*y_coord[0];
	//使z标准化
	double znorm = sqrt(z_coord[0]*z_coord[0] + z_coord[1]*z_coord[1] + z_coord[2]*z_coord[2]);
	z_coord[0] = z_coord[0]/znorm;
	z_coord[1] = z_coord[1]/znorm;
	z_coord[2] = z_coord[2]/znorm;

	//修正x方向
	double delt = (x_coord[0] * z_coord[0]+ x_coord[1] * z_coord[1] + x_coord[2] * z_coord[2])/sqrt(z_coord[0] * z_coord[0] + z_coord[1] * z_coord[1] + z_coord[2] * z_coord[2]);
	x_coord[0] = x_coord[0] - delt*z_coord[0];
	x_coord[1] = x_coord[1] - delt*z_coord[1];
	x_coord[2] = x_coord[2] - delt*z_coord[2];
	//使x标准化
	double xnorm = sqrt(x_coord[0]*x_coord[0] + x_coord[1]*x_coord[1] + x_coord[2]*x_coord[2]);
	x_coord[0] = x_coord[0]/xnorm;
	x_coord[1] = x_coord[1]/xnorm;
	x_coord[2] = x_coord[2]/xnorm;

	//修正y方向
	y_coord[0] = x_coord[1]*z_coord[2] - x_coord[2]*z_coord[1];
    y_coord[1] = x_coord[2]*z_coord[0] - x_coord[0]*z_coord[2];
	y_coord[2] = x_coord[0]*z_coord[1] - x_coord[1]*z_coord[0];
	//使y标准化
	double ynorm = sqrt(y_coord[0]*y_coord[0] + y_coord[1]*y_coord[1] + y_coord[2]*y_coord[2]);
	y_coord[0] = y_coord[0]/ynorm;
	y_coord[1] = y_coord[1]/ynorm;
	y_coord[2] = y_coord[2]/ynorm;

	coordbuildstep = 2;
}
/*
void PlaneManager::getCoordinary()
{
	
	std::sort(cloud_planes.begin(),cloud_planes.end(),planeInfoCompare);  
	for(vector<plane_info>::const_iterator iter = cloud_planes.begin(); iter != cloud_planes.end(); ++iter) {
		 plane_info p = *iter;
		 if(coordbuildstep == 0){
			 x_coord[0] = p.a;
			 x_coord[1] = p.b;
			 x_coord[2] = p.c;
			 double norm = sqrt(x_coord[0]*x_coord[0] + x_coord[1]*x_coord[1] + x_coord[2]*x_coord[2]);
			 x_coord[0] = x_coord[0]/norm;
			 x_coord[1] = x_coord[1]/norm;
			 x_coord[2] = x_coord[2]/norm;
			 coordbuildstep = 1;
		 }else if(coordbuildstep == 1){
			 double angle =  (x_coord[0] * p.a + x_coord[1] * p.b + x_coord[2] * p.c)/sqrt(x_coord[0] * x_coord[0] + x_coord[1] * x_coord[1] + x_coord[2] * x_coord[2])/sqrt(p.a * p.a + p.b * p.b + p.c * p.c);
		     double delt = (x_coord[0] * p.a + x_coord[1] * p.b + x_coord[2] * p.c)/sqrt(x_coord[0] * x_coord[0] + x_coord[1] * x_coord[1] + x_coord[2] * x_coord[2]);
			 if(angle < 0.1 && angle > -0.1){
				 y_coord[0] = p.a;
				 y_coord[1] = p.b;
				 y_coord[2] = p.c;
				 double norm = sqrt(y_coord[0]*y_coord[0] + y_coord[1]*y_coord[1] + y_coord[2]*y_coord[2]);
				 y_coord[0] = y_coord[0]/norm;
				 y_coord[1] = y_coord[1]/norm;
				 y_coord[2] = y_coord[2]/norm;
				 coordbuildstep = 2;
				 break;
			 }
		 }
	}
	//求z方向
	//叉乘x方向和y方向
	z_coord[0] = x_coord[1]*y_coord[2] - x_coord[2]*y_coord[1];
    z_coord[1] = x_coord[2]*y_coord[0] - x_coord[0]*y_coord[2];
	z_coord[2] = x_coord[0]*y_coord[1] - x_coord[1]*y_coord[0];
	//使z标准化
	double znorm = sqrt(z_coord[0]*z_coord[0] + z_coord[1]*z_coord[1] + z_coord[2]*z_coord[2]);
	z_coord[0] = z_coord[0]/znorm;
	z_coord[1] = z_coord[1]/znorm;
	z_coord[2] = z_coord[2]/znorm;

	//修正x方向
	double delt = (x_coord[0] * z_coord[0]+ x_coord[1] * z_coord[1] + x_coord[2] * z_coord[2])/sqrt(z_coord[0] * z_coord[0] + z_coord[1] * z_coord[1] + z_coord[2] * z_coord[2]);
	x_coord[0] = x_coord[0] - delt*z_coord[0];
	x_coord[1] = x_coord[1] - delt*z_coord[1];
	x_coord[2] = x_coord[2] - delt*z_coord[2];
	//使x标准化
	double xnorm = sqrt(x_coord[0]*x_coord[0] + x_coord[1]*x_coord[1] + x_coord[2]*x_coord[2]);
	x_coord[0] = x_coord[0]/xnorm;
	x_coord[1] = x_coord[1]/xnorm;
	x_coord[2] = x_coord[2]/xnorm;

	//修正y方向
	y_coord[0] = x_coord[1]*z_coord[2] - x_coord[2]*z_coord[1];
    y_coord[1] = x_coord[2]*z_coord[0] - x_coord[0]*z_coord[2];
	y_coord[2] = x_coord[0]*z_coord[1] - x_coord[1]*z_coord[0];
	//使y标准化
	double ynorm = sqrt(y_coord[0]*y_coord[0] + y_coord[1]*y_coord[1] + y_coord[2]*y_coord[2]);
	y_coord[0] = y_coord[0]/ynorm;
	y_coord[1] = y_coord[1]/ynorm;
	y_coord[2] = y_coord[2]/ynorm;

	coordbuildstep = 2;
}
*/

Point3D PlaneManager::getCrossPoint(Plane3D p1,Plane3D p2,Plane3D p3){
	double d = det(p1.a,p1.b,p1.c,p2.a,p2.b,p2.c,p3.a,p3.b,p3.c);
	double dx = det(p1.d,p1.b,p1.c,p2.d,p2.b,p2.c,p3.d,p3.b,p3.c);
	double dy = det(p1.a,p1.d,p1.c,p2.a,p2.d,p2.c,p3.a,p3.d,p3.c);
	double dz = det(p1.a,p1.b,p1.d,p2.a,p2.b,p2.d,p3.a,p3.b,p3.d);
	Point3D ret = *(new Point3D());
	ret.x = dx/d;
	ret.y = dy/d;
	ret.z = dz/d;
	return ret;
}

double PlaneManager::det(double a1,double b1,double c1,double a2,double b2,double c2,double a3,double b3,double c3){
	return a1 * b2 * c3 +c1 * a2 * b3 + b1 * c2 * a3 - c1 * b2 *a3 - b1 * a2 * c3 - a1 * c2 * b3;
}

void PlaneManager::generateSkeleton(){
	if(coordbuildstep < 2) getCoordinary();
	for(int i = 0;i < clutered_planes->size();i++) {
		 plane_info p = clutered_planes->at(i);
		 clutered_planes->at(i).skeletonGenerated = false;
		 if(p.inliers > minInliers){
			 //确定扫描线
			 double xVec = p.a * x_coord[0] + p.b * x_coord[1] + p.c * x_coord[2];
			 double yVec = p.a * y_coord[0] + p.b * y_coord[1] + p.c * y_coord[2];
			 double zVec = p.a * z_coord[0] + p.b * z_coord[1] + p.c * z_coord[2];
			 double x_scan[3];
			 double y_scan[3];
			 double normal[3];
			 if(zVec < 0.2 && zVec > -0.2){
				 y_scan[0] = z_coord[0];
				 y_scan[1] = z_coord[1];
				 y_scan[2] = z_coord[2];
				 if(xVec < 0.2 && xVec > -0.2){
					 x_scan[0] = x_coord[0];
					 x_scan[1] = x_coord[1];
					 x_scan[2] = x_coord[2];

					 normal[0] = y_coord[0];
					 normal[1] = y_coord[1];
					 normal[2] = y_coord[2];

				 }else if(yVec < 0.2 && yVec > -0.2){
					 x_scan[0] = y_coord[0];
					 x_scan[1] = y_coord[1];
					 x_scan[2] = y_coord[2];

					 normal[0] = x_coord[0];
					 normal[1] = x_coord[1];
					 normal[2] = x_coord[2];
				 }else{
					 continue;
				 }
			 }else if(xVec < 0.2 && xVec > -0.2 && yVec < 0.2 && yVec > -0.2){
				 x_scan[0] = x_coord[0];
				 x_scan[1] = x_coord[1];
				 x_scan[2] = x_coord[2];
				 y_scan[0] = y_coord[0];
				 y_scan[1] = y_coord[1];
				 y_scan[2] = y_coord[2];

				normal[0] = z_coord[0];
				normal[1] = z_coord[1];
				normal[2] = z_coord[2];

			 }else{
				 continue;
			 }

			 //扫描确定边界 粗略版 确定矩形边界
			 double xMax = 0,xMin = 0,yMax = 0,yMin = 0,normD = 0;
			 bool initial = false;
			 for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = p.bounds->begin();it != p.bounds->end();it++){
			   pcl::PointXYZRGB point = *it;
			   double d1 = (point.x * x_scan[0] + point.y * x_scan[1] + point.z * x_scan[2]);
			   double d2 = (point.x * y_scan[0] + point.y * y_scan[1] + point.z * y_scan[2]);

			   double dn = (point.x * normal[0] + point.y * normal[1] + point.z * normal[2]);
			   normD += dn;

			   if(!initial){
				   xMax = xMin = d1;
				   yMax = yMin = d2;
				   initial = true;
			   }else{
				   if(d1 > xMax) xMax = d1;
				   if(d1 < xMin) xMin = d1;
				   if(d2 > yMax) yMax = d2;
				   if(d2 < yMin) yMin = d2;
			   }

		   }
		   //确定矩形边界4点：
			 Plane3D pp = *(new Plane3D);
			 pp.a = normal[0];
			 pp.b = normal[1];
			 pp.c = normal[2];
			 pp.d = normD/p.bounds->size();


			 Plane3D px1 = *(new Plane3D);
			 px1.a = x_scan[0];
			 px1.b = x_scan[1];
			 px1.c = x_scan[2];
			 px1.d = xMax;

			 Plane3D px2 = *(new Plane3D);
			 px2.a = x_scan[0];
			 px2.b = x_scan[1];
			 px2.c = x_scan[2];
			 px2.d = xMin;

			 Plane3D py1 = *(new Plane3D);
			 py1.a = y_scan[0];
			 py1.b = y_scan[1];
			 py1.c = y_scan[2];
			 py1.d = yMax;

			 Plane3D py2 = *(new Plane3D);
			 py2.a = y_scan[0];
			 py2.b = y_scan[1];
			 py2.c = y_scan[2];
			 py2.d = yMin;

			 Point3D cross0 = getCrossPoint(pp,px1,py1);
			 Point3D cross1 = getCrossPoint(pp,py1,px2);
             Point3D cross2 = getCrossPoint(pp,px2,py2);
			 Point3D cross3 = getCrossPoint(pp,py2,px1);

			 clutered_planes->at(i).boundPoints = new Point3D[4];
			 clutered_planes->at(i).boundPoints[0].x = cross0.x;
			 clutered_planes->at(i).boundPoints[0].y = cross0.y;
			 clutered_planes->at(i).boundPoints[0].z = cross0.z;
			 clutered_planes->at(i).boundPoints[1].x = cross1.x;
			 clutered_planes->at(i).boundPoints[1].y = cross1.y;
			 clutered_planes->at(i).boundPoints[1].z = cross1.z;
			 clutered_planes->at(i).boundPoints[2].x = cross2.x;
			 clutered_planes->at(i).boundPoints[2].y = cross2.y;
			 clutered_planes->at(i).boundPoints[2].z = cross2.z;
			 clutered_planes->at(i).boundPoints[3].x = cross3.x;
			 clutered_planes->at(i).boundPoints[3].y = cross3.y;
			 clutered_planes->at(i).boundPoints[3].z = cross3.z;
			 
			 clutered_planes->at(i).a = pp.a;
			 clutered_planes->at(i).b = pp.b;
			 clutered_planes->at(i).c = pp.c;
			 clutered_planes->at(i).d = pp.d;
			 clutered_planes->at(i).skeletonGenerated = true;
		 }
	}	
}
/*
void PlaneManager::generateSkeleton(){
	if(coordbuildstep < 2) getCoordinary();
	for(int i = 0;i < cloud_planes.size();i++) {
		 plane_info p = cloud_planes.at(i);
		 cloud_planes.at(i).skeletonGenerated = false;
		 if(p.inliers > minInliers){
			 //确定扫描线
			 double xVec = p.a * x_coord[0] + p.b * x_coord[1] + p.c * x_coord[2];
			 double yVec = p.a * y_coord[0] + p.b * y_coord[1] + p.c * y_coord[2];
			 double zVec = p.a * z_coord[0] + p.b * z_coord[1] + p.c * z_coord[2];
			 double x_scan[3];
			 double y_scan[3];
			 if(xVec < 0.1 && xVec > -0.1){
				 x_scan[0] = x_coord[0];
				 x_scan[1] = x_coord[1];
				 x_scan[2] = x_coord[2];
				 if(yVec < 0.1 && yVec > -0.1){
					 y_scan[0] = y_coord[0];
					 y_scan[1] = y_coord[1];
					 y_scan[2] = y_coord[2];
				 }else{
					 y_scan[0] = z_coord[0];
					 y_scan[1] = z_coord[1];
					 y_scan[2] = z_coord[2];
				 }
			 }else{
				 x_scan[0] = y_coord[0];
				 x_scan[1] = y_coord[1];
				 x_scan[2] = y_coord[2];
				 y_scan[0] = z_coord[0];
				 y_scan[1] = z_coord[1];
				 y_scan[2] = z_coord[2];
			 }

			 //扫描确定边界 粗略版 确定矩形边界
			 double xMax = 0,xMin = 0,yMax = 0,yMin = 0;
			 bool initial = false;
			 for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = p.bounds->begin();it != p.bounds->end();it++){
			   pcl::PointXYZRGB point = *it;
			   double d1 = -(point.x * x_scan[0] + point.y * x_scan[1] + point.z * x_scan[2]);
			   double d2 = -(point.x * y_scan[0] + point.y * y_scan[1] + point.z * y_scan[2]);
			   if(!initial){
				   xMax = xMin = d1;
				   yMax = yMin = d2;
				   initial = true;
			   }else{
				   if(d1 > xMax) xMax = d1;
				   if(d1 < xMin) xMin = d1;
				   if(d2 > yMax) yMax = d2;
				   if(d2 < yMin) yMin = d2;
			   }

		   }
		   //确定矩形边界4点：
			 Plane3D pp = *(new Plane3D);
			 pp.a = p.a;
			 pp.b = p.b;
			 pp.c = p.c;
			 pp.d = -p.d;

			 Plane3D px1 = *(new Plane3D);
			 px1.a = x_scan[0];
			 px1.b = x_scan[1];
			 px1.c = x_scan[2];
			 px1.d = xMax;

			 Plane3D px2 = *(new Plane3D);
			 px2.a = x_scan[0];
			 px2.b = x_scan[1];
			 px2.c = x_scan[2];
			 px2.d = xMin;

			 Plane3D py1 = *(new Plane3D);
			 py1.a = y_scan[0];
			 py1.b = y_scan[1];
			 py1.c = y_scan[2];
			 py1.d = yMax;

			 Plane3D py2 = *(new Plane3D);
			 py2.a = y_scan[0];
			 py2.b = y_scan[1];
			 py2.c = y_scan[2];
			 py2.d = yMin;

			 Point3D cross0 = getCrossPoint(pp,px1,py1);
			 Point3D cross1 = getCrossPoint(pp,py1,px2);
             Point3D cross2 = getCrossPoint(pp,px2,py2);
			 Point3D cross3 = getCrossPoint(pp,py2,px1);

			 cloud_planes.at(i).boundPoints = new Point3D[4];
			 cloud_planes.at(i).boundPoints[0].x = cross0.x;
			 cloud_planes.at(i).boundPoints[0].y = cross0.y;
			 cloud_planes.at(i).boundPoints[0].z = cross0.z;
			 cloud_planes.at(i).boundPoints[1].x = cross1.x;
			 cloud_planes.at(i).boundPoints[1].y = cross1.y;
			 cloud_planes.at(i).boundPoints[1].z = cross1.z;
			 cloud_planes.at(i).boundPoints[2].x = cross2.x;
			 cloud_planes.at(i).boundPoints[2].y = cross2.y;
			 cloud_planes.at(i).boundPoints[2].z = cross2.z;
			 cloud_planes.at(i).boundPoints[3].x = cross3.x;
			 cloud_planes.at(i).boundPoints[3].y = cross3.y;
			 cloud_planes.at(i).boundPoints[3].z = cross3.z;

			 cloud_planes.at(i).skeletonGenerated = true;
		 }
	}
}
*/







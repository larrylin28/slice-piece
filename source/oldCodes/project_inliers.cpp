#include <iostream>
#include <string>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/filters/project_inliers.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/range_image_visualizer.h>
#include <pcl/point_cloud.h>
#include <XnCppWrapper.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/filters/passthrough.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/organized_multi_plane_segmentation.h>
#include <pcl/features/integral_image_normal.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d.h>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/legacy/legacy.hpp>
#include <pcl/registration/gicp.h>
#include "slam3D.h"
//#include "KeyPointManager.h"
#include "PlaneManager.h"
#include <time.h>

#include "glwindows.h"
#include "NeHeGL.h"
#include "sdf.h"
#include "icp.h"

using namespace std;


#define GETCOLORR(i) (i%4)*1.0/3
#define GETCOLORG(i) ((i+1)%5)*1.0/4
#define GETCOLORB(i) ((1000-i)%4)*1.0/3
typedef pcl::PointXYZRGB PointT;
typedef pcl::PointXYZRGBNormal PointN;
typedef pcl::PointCloud<PointT> PointCloud;

vector<IplImage*> depths,images;
vector<XnUInt64> timestamps;
vector<XnInt32> frameIDs;
double XtoZ, YtoZ;
vector<const XnDepthPixel *> depth_maps;
vector<xn::DepthMetaData *> depth_maps1;
vector<const XnRGB24Pixel *> image_maps;
vector<xn::ImageMetaData *> image_maps1;

//SLAM用的存储结构
vector<vector<cv::KeyPoint>> keypoints; //存储每帧包含的SURF特征点
//vector<vector<vector<KeyPointID>>> keyMatches; //存储每帧包含的特征点所匹配的特征点
vector<vector<int>> foundPlanes;  //存储每帧提取出来的平面对应的ID

vector<cv::Mat> descriptors;
vector<PointCloud::Ptr> pointClouds;
vector<PointCloud::Ptr> transpointClouds;
vector<pcl::PointCloud<pcl::Normal>::Ptr> transpointNormals;
vector<Eigen::Matrix4f> trans;

#define isnan(x) ((x) != (x)) 

struct frameMean{
	int frame;
	double m[3];
};

void put_Point(float* m,PointT pt,pcl::Normal pn){
	m[0] = pt.x;
	m[1] = pt.y;
	m[2] = pt.z;
	m[3] = pn.normal_x;
	m[4] = pn.normal_y;
	m[5] = pn.normal_z;
}

bool lessmark0(const frameMean& Item1, const frameMean& Item2)  
 {  
     return Item1.m[0] > Item2.m[0];  
 }

bool lessmark1(const frameMean& Item1, const frameMean& Item2)  
 {  
     return Item1.m[1] > Item2.m[1];  
 }

bool lessmark2(const frameMean& Item1, const frameMean& Item2)  
 {  
     return Item1.m[2] > Item2.m[2];  
 }


void CheckOpenNIError( XnStatus result, string status )  
{   
	if( result != XN_STATUS_OK )   
		cerr << status << " Error: " << xnGetStatusString( result ) << endl;  
} 



double* CrossProduct(double a[3], double b[3])
{
    double* c = new double[3];

    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];

    return c;
}

double DotProduct(double a[3], double b[3])
{
    double result;
    result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    return result;
}

double Normalize(double v[3])
{
    double result;

    result = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    return result;
}

double* Rotate(double** rotateMatrix, double u[3]){
	double* ret = new double[3];
    ret[0]=rotateMatrix[0][0]*u[0]+rotateMatrix[0][1]*u[1]+rotateMatrix[0][2]*u[2];
    ret[1]=rotateMatrix[1][0]*u[0]+rotateMatrix[1][1]*u[1]+rotateMatrix[1][2]*u[2];
    ret[2]=rotateMatrix[2][0]*u[0]+rotateMatrix[2][1]*u[1]+rotateMatrix[2][2]*u[2];

	return ret;
}

double** RotationMatrix(double angle, double u[3])
{
    double norm = Normalize(u);
    double** rotatinMatrix = new double*[3];
	for(int id = 0; id < 3;id++){
		rotatinMatrix[id] = new double[3];
	}
    
    u[0] = u[0] / norm;
    u[1] = u[1] / norm;
    u[2] = u[2] / norm;

    rotatinMatrix[0][0] = cos(angle) + u[0] * u[0] * (1 - cos(angle));
    rotatinMatrix[0][1] = u[0] * u[1] * (1 - cos(angle) - u[2] * sin(angle));
    rotatinMatrix[0][2] = u[1] * sin(angle) + u[0] * u[2] * (1 - cos(angle));

    rotatinMatrix[1][0] = u[2] * sin(angle) + u[0] * u[1] * (1 - cos(angle));
    rotatinMatrix[1][1] = cos(angle) + u[1] * u[1] * (1 - cos(angle));
    rotatinMatrix[1][2] = -u[0] * sin(angle) + u[1] * u[2] * (1 - cos(angle));
      
    rotatinMatrix[2][0] = -u[1] * sin(angle) + u[0] * u[2] * (1 - cos(angle));
    rotatinMatrix[2][1] = u[0] * sin(angle) + u[1] * u[2] * (1 - cos(angle));
    rotatinMatrix[2][2] = cos(angle) + u[2] * u[2] * (1 - cos(angle));

    return rotatinMatrix;
}

//求旋转矩阵
double** Calculation(double vectorBefore[3], double vectorAfter[3])
{
    double* rotationAxis;
    double rotationAngle;
    rotationAxis = CrossProduct(vectorBefore, vectorAfter);
    rotationAngle = acos(DotProduct(vectorBefore, vectorAfter) / Normalize(vectorBefore) / Normalize(vectorAfter));
    return RotationMatrix(rotationAngle, rotationAxis);
}

Eigen::Matrix4f CalTransport(double vectorBefore[4],double vectorAfter[4]){
	double vb[3];
	double va[3];
	for(int i = 0;i < 3;i++){
		vb[i] = vectorBefore[i];
		va[i] = vectorAfter[i];
	}
	double** rotateMatrix = Calculation(vb, va);

	Eigen::Matrix4f p_transformation_matrix;
	
	p_transformation_matrix(0,0) = rotateMatrix[0][0];
	p_transformation_matrix(0,1) = rotateMatrix[0][1];
	p_transformation_matrix(0,2) = rotateMatrix[0][2];
	p_transformation_matrix(0,3) = 0;

	p_transformation_matrix(1,0) = rotateMatrix[1][0];
	p_transformation_matrix(1,1) = rotateMatrix[1][1];
	p_transformation_matrix(1,2) = rotateMatrix[1][2];
	p_transformation_matrix(1,3) = 0;

	p_transformation_matrix(2,0) = rotateMatrix[2][0];
	p_transformation_matrix(2,1) = rotateMatrix[2][1];
	p_transformation_matrix(2,2) = rotateMatrix[2][2];
	p_transformation_matrix(2,3) = 0;

	p_transformation_matrix(3,0) = 0;
	p_transformation_matrix(3,1) = 0;
	p_transformation_matrix(3,2) = 0;
	p_transformation_matrix(3,3) = 1;

	return p_transformation_matrix;

}

double** MutiMatrix(double** a, double** b, int m,int n,int l){

	double** ret= new double*[n];
	for(int id = 0; id < n;id++){
		ret[id] = new double[l];
	}
	

	for(int i = 0;i < n;i++){
		for(int j = 0;j < l;j++){
			 ret[i][j] = 0;
			 for(int k = 0;k < m;k++){
				 ret[i][j] += a[i][k] * b[k][j];
			 }
		}
	}
	return ret;
}

double** MutiMatrix3(double** a, double** b){
    return MutiMatrix(a,b,3,3,3);
}

pcl::PointCloud<pcl::PointXYZRGB>::Ptr convertToPointCloud(xn::DepthGenerator& rDepthGen,const XnDepthPixel *dm,const XnRGB24Pixel *im,XnUInt64 timestamp, XnInt32 frameID)		
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);

	// Not supported in file yet:
	cloud->header.stamp = timestamp;
	cloud->header.seq = frameID;
	// End not supported in file yet
	cloud->height = 480;
	cloud->width = 640;
	cloud->is_dense = false;

	cloud->points.resize(cloud->height * cloud->width);

	register int centerX = (cloud->width >> 1);
	int centerY = (cloud->height >> 1);

	register const XnDepthPixel* depth_map = dm;
	register const XnRGB24Pixel* rgb_map = im;

  unsigned int uPointNum = cloud->width * cloud->height;


  XnPoint3D* pDepthPointSet = new XnPoint3D[ uPointNum ];
  unsigned int i, j, idxShift, idx;
  for( j = 0; j < cloud->height; ++j )
  {
    idxShift = j * cloud->width;
    for( i = 0; i < cloud->width; ++i )
    {
      idx = idxShift + i;
      pDepthPointSet[idx].X = i;
      pDepthPointSet[idx].Y = j;
      pDepthPointSet[idx].Z = depth_map[idx];
    }
  }

  XnPoint3D* p3DPointSet = new XnPoint3D[ uPointNum ];
  rDepthGen.ConvertProjectiveToRealWorld( uPointNum, pDepthPointSet, p3DPointSet );
  delete[] pDepthPointSet;
  register int depth_idx = 0;
  for (int v = -centerY; v < centerY; ++v)
	{
		for (register int u = -centerX; u < centerX; ++u, ++depth_idx)
		{
			pcl::PointXYZRGB& pt = cloud->points[depth_idx];

			pt.z = p3DPointSet[depth_idx].Z * 0.001f;
			pt.x = p3DPointSet[depth_idx].X * 0.001f;
			pt.y = p3DPointSet[depth_idx].Y * 0.001f;
			pt.r = (float) rgb_map[depth_idx].nRed;
			pt.g = (float) rgb_map[depth_idx].nGreen;
			pt.b = (float) rgb_map[depth_idx].nBlue;
		}
	}
	return cloud;
}

vector<plane_info>* KMeansCluster(vector<plane_info>* planes,double treshold,double iterbound)
{
   //确定K值和初始点
	//double treshold = 0.05;
	if(planes->size() == 0) return new vector<plane_info>();
	int K = 1;
	vector<plane_info>* mean = new vector<plane_info>();
    plane_info* m = new plane_info();
	m->a = planes->at(0).a;
	m->b = planes->at(0).b;
	m->c = planes->at(0).c;
	m->d = planes->at(0).d;
	mean->push_back(*m);
	for(int i = 1;i < planes->size();i++){
		bool addK = true;
		for(int j = 0;j < K;j++){
			//double dis = pcl::euclideanDistance(p, mean->at(j));
			double dis = sqrt((mean->at(j).a - planes->at(i).a)*(mean->at(j).a - planes->at(i).a) + (mean->at(j).b - planes->at(i).b)*(mean->at(j).b - planes->at(i).b) + (mean->at(j).c - planes->at(i).c)*(mean->at(j).c - planes->at(i).c) + (mean->at(j).d - planes->at(i).d)*(mean->at(j).d - planes->at(i).d));
			if(dis < treshold){
				addK = false;
				break;
			}
		}
		if(addK){
			K++;
			plane_info* m = new plane_info();
			m->a = planes->at(i).a;
			m->b = planes->at(i).b;
			m->c = planes->at(i).c;
			m->d = planes->at(i).d;
			mean->push_back(*m);
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
		for(int i = 0;i < planes->size();i++){
			for(int j = 0;j < K;j++){
				double dis = sqrt((mean->at(j).a - planes->at(i).a)*(mean->at(j).a - planes->at(i).a) + (mean->at(j).b - planes->at(i).b)*(mean->at(j).b - planes->at(i).b) + (mean->at(j).c - planes->at(i).c)*(mean->at(j).c - planes->at(i).c) + (mean->at(j).d - planes->at(i).d)*(mean->at(j).d - planes->at(i).d));
				if(dis < treshold){
					count[j] += planes->at(i).inliers;
					amean[j] += planes->at(i).a * planes->at(i).inliers;
					bmean[j] += planes->at(i).b * planes->at(i).inliers;
					cmean[j] += planes->at(i).c * planes->at(i).inliers;
					dmean[j] += planes->at(i).d * planes->at(i).inliers;
					break;
				}
			}
		}
		bool iterend = true;
		for(int i = 0;i < K;i++){
			amean[i] /= count[i];
			bmean[i] /= count[i];
			cmean[i] /= count[i];
			dmean[i] /= count[i];
			if(amean[i] - mean->at(i).a > iterbound || bmean[i] - mean->at(i).b > iterbound && cmean[i] - mean->at(i).c > iterbound && dmean[i] - mean->at(i).d > iterbound){
				iterend = false;
			}
			mean->at(i).a = amean[i];
			mean->at(i).b = bmean[i];
			mean->at(i).c = cmean[i];
			mean->at(i).d = dmean[i];
		}
		if(iterend){ break;}
	}
	return mean;
}



double* ransacmethod(double ddx[],double ddy[],double ddz[],int length,double theshold){
	double* pnow = new double[3];
	//double distance;
	//double mindistance = 10000;
	int maxinliercount = -1;
	bool* mininlier;
	int minid = -1;
	//找最优点及其内点
	for(int i = 0;i < length;i++){
		//distance = 0;
		int inliercount = 0;
		bool* inlier = new bool[length];
		for(int j = 0;j < length;j++){
			//计算distance
			double dis = sqrt((ddx[i]-ddx[j])*(ddx[i]-ddx[j]) + (ddy[i]-ddy[j])*(ddy[i]-ddy[j]) + (ddz[i]-ddz[j])*(ddz[i]-ddz[j]));
			//distance += dis;
			if(dis < theshold){
			    inlier[j] = true;
				inliercount++;
			}else{
				inlier[j] = false;
			}
		}
		if(inliercount > maxinliercount){
			maxinliercount = inliercount;
			//mindistance = distance;
			minid = i;
			mininlier = inlier;
		}
	}
	//由内点估算最优模型
	pnow[0] = 0;
	pnow[1] = 0;
	pnow[2] = 0;
	int count = 0;
	for(int i = 0;i < length;i++){
		if(mininlier[i]){
			pnow[0] += ddx[i];
			pnow[1] += ddy[i];
			pnow[2] += ddz[i];
			count++;
		}
	}
	pnow[0] /= count;
	pnow[1] /= count;
	pnow[2] /= count;
	return pnow;
}

/*
int findMatchKeypoint(int frameID,int keyID,int planeID,PlaneManager p_manager){
	int size = keyMatches[frameID][keyID].size();
	//if(planeID >= p_manager.virtual_planes->size()) return -1;
	int pksize = p_manager.virtual_planes->at(planeID).virtualKeypointsCount;
	if(size <= 0 || pksize  <= 0) return -1;
	int* count = new int[pksize];
	for(int i = 0;i < pksize;i++){
		count[i] = 0;
	}
	for(int i = 0;i < size;i++){
		KeyPointID id = keyMatches[frameID][keyID][i];
		if(id.frameID >= frameID) break;
		virtualKeyPointID pmatch = planeMatches[id.frameID][id.keyID];
		if(pmatch.planeID == planeID){
			count[pmatch.keyID]++;
		}
	}
	int maxcount = 0;
	int maxID = -1;
	for(int i = 0;i < pksize;i++){
		if(count[i] > maxcount){
			maxID = i;
			maxcount = count[i];
		}
	}
	if(maxcount > 0) return maxID;
	else return -1;
}


bool do_match_with_plane(vector<point2PlaneMatch> matches,Eigen::Matrix4f &transformation_matrix,int id_s,PlaneManager p_manager){
	pcl::CorrespondencesPtr correspondences (new pcl::Correspondences);
	correspondences->resize(matches.size());

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
	keypoints_source->width    = 1*matches.size();
	keypoints_source->height   = 1;
	keypoints_source->is_dense = false;
	keypoints_source->points.resize (keypoints_source->width * keypoints_source->height);

	keypoints_target->width    = 1*matches.size();
	keypoints_target->height   = 1;
	keypoints_target->is_dense = false;
	keypoints_target->points.resize (keypoints_target->width * keypoints_target->height);

	for(int i = 0;i < matches.size();i++){
		double x1 = keypoints[id_s][matches[i].keypoint].pt.x;
		double y1 = keypoints[id_s][matches[i].keypoint].pt.y;
		PointT searchPoint_s = pointClouds[id_s]->points[(int)(y1+0.5)*pointClouds[id_s]->width+(int)(x1+0.5)];
		pcl::PointXYZRGB& pt_s = keypoints_source->points[i];
		pt_s.x = searchPoint_s.x;
		pt_s.y = searchPoint_s.y;
		pt_s.z = searchPoint_s.z;

		pcl::PointXYZRGB& pt_t = keypoints_target->points[i];
		virtualPlane plane = p_manager.virtual_planes->at(matches[i].plane);
		pt_t.x = plane.virtualKeypoints[matches[i].virtualPoint].currentMean.x;
		pt_t.y = plane.virtualKeypoints[matches[i].virtualPoint].currentMean.y;
		pt_t.z = plane.virtualKeypoints[matches[i].virtualPoint].currentMean.z;

		(*correspondences)[i].index_query = i;
		(*correspondences)[i].index_match = i;
		(*correspondences)[i].distance = 
		pcl::euclideanDistance(keypoints_source->at(i), keypoints_target->at(i));
	}

	pcl::CorrespondencesPtr ransac_correspondences (new pcl::Correspondences);
	double inlier_threshold = 0.005;
	int iteration_time = 1000;
	pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZRGB> rejector;
	rejector.setInputCloud(keypoints_source);
	rejector.setTargetCloud(keypoints_target);
	rejector.setInputCorrespondences(correspondences);
	rejector.setInlierThreshold(inlier_threshold);
	rejector.setMaxIterations(iteration_time);
	rejector.getCorrespondences(*ransac_correspondences);

	Eigen::Matrix4f trans = rejector.getBestTransformation();
	if(ransac_correspondences->size() < 3) return false;
	else{
		transformation_matrix << trans;
		return true;
	}
	return false;
}
*/
bool try_do_match(int id_s,int id_t,pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source,pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target,double &variance,vector<PointCloud::Ptr>* sourceClouds, vector<PointCloud::Ptr>* targetClouds,int id_t_incloud){

	cv::BruteForceMatcher<cv::L2<float>> matcher;
	vector<cv::DMatch> matches;
	matcher.match(descriptors[id_s],descriptors[id_t],matches);
	keypoints_source->width    = 1*matches.size();
	keypoints_source->height   = 1;
	keypoints_source->is_dense = false;
	keypoints_source->points.resize (keypoints_source->width * keypoints_source->height);

	keypoints_target->width    = 1*matches.size();
	keypoints_target->height   = 1;
	keypoints_target->is_dense = false;
	keypoints_target->points.resize (keypoints_target->width * keypoints_target->height);

	
	pcl::CorrespondencesPtr correspondences (new pcl::Correspondences);
	correspondences->resize(matches.size());

	for(int i = 0;i < matches.size();i++){
		
		double x1 = keypoints[id_s][matches[i].queryIdx].pt.x;
		double y1 = keypoints[id_s][matches[i].queryIdx].pt.y;
		double x2 = keypoints[id_t][matches[i].trainIdx].pt.x;
		double y2 = keypoints[id_t][matches[i].trainIdx].pt.y;

		PointT searchPoint_s = sourceClouds->at(id_s)->points[(int)(y1+0.5)* sourceClouds->at(id_s)->width+(int)(x1+0.5)];
		pcl::PointXYZRGB& pt_s = keypoints_source->points[i];
		pt_s.x = searchPoint_s.x;
		pt_s.y = searchPoint_s.y;
		pt_s.z = searchPoint_s.z;
		pt_s.r = searchPoint_s.r;
		pt_s.g = searchPoint_s.g;
		pt_s.b = searchPoint_s.b;


		PointT searchPoint_t = targetClouds->at(id_t_incloud)->points[(int)(y2+0.5)*targetClouds->at(id_t_incloud)->width+(int)(x2+0.5)];
		pcl::PointXYZRGB& pt_t = keypoints_target->points[i];
		pt_t.x = searchPoint_t.x;
		pt_t.y = searchPoint_t.y;
		pt_t.z = searchPoint_t.z;
		pt_t.r = searchPoint_t.r;
		pt_t.g = searchPoint_t.g;
		pt_t.b = searchPoint_t.b;

		(*correspondences)[i].index_query = i;
		(*correspondences)[i].index_match = i;
		(*correspondences)[i].distance = 
		pcl::euclideanDistance(keypoints_source->at(i), keypoints_target->at(i));
	}

	double inlier_threshold = 0.005;
	int iteration_time = 1000;
	pcl::CorrespondencesPtr ransac_correspondences (new pcl::Correspondences);
	pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZRGB> rejector;
	rejector.setInputCloud(keypoints_source);
	rejector.setTargetCloud(keypoints_target);
	rejector.setInputCorrespondences(correspondences);
	rejector.setInlierThreshold(inlier_threshold);
	rejector.setMaxIterations(iteration_time);
	rejector.getCorrespondences(*ransac_correspondences);

	Eigen::Matrix4f trans = rejector.getBestTransformation();

	if(ransac_correspondences->size() < 3) return false;
	else{
		
		variance = 0;
		for(int i = 0;i < ransac_correspondences->size();i++){
			double dist = ransac_correspondences->at(i).distance;
			double weight =  ransac_correspondences->at(i).weight;
			variance += dist*dist*weight;
		}
		return true;
	}
	
}

bool real_do_match(pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source,pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target,Eigen::Matrix4f &p_transformation_matrix){

	int size = keypoints_source->size();
	pcl::CorrespondencesPtr correspondences (new pcl::Correspondences);
	correspondences->resize(size);
	for(int i = 0;i < size;i++){
		(*correspondences)[i].index_query = i;
		(*correspondences)[i].index_match = i;
		(*correspondences)[i].distance = 
		pcl::euclideanDistance(keypoints_source->at(i), keypoints_target->at(i));
	}

	double inlier_threshold = 0.005;
	int iteration_time = 1000;
	pcl::CorrespondencesPtr ransac_correspondences (new pcl::Correspondences);
	pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZRGB> rejector;
	rejector.setInputCloud(keypoints_source);
	rejector.setTargetCloud(keypoints_target);
	rejector.setInputCorrespondences(correspondences);
	rejector.setInlierThreshold(inlier_threshold);
	rejector.setMaxIterations(iteration_time);
	rejector.getCorrespondences(*ransac_correspondences);

	Eigen::Matrix4f trans = rejector.getBestTransformation();

	if(ransac_correspondences->size() < 3) return false;
	else{
		p_transformation_matrix << trans;
		return true;
	}
	
}

bool do_match(int id_s,int id_t,pcl::CorrespondencesPtr &ransac_correspondences,Eigen::Matrix4f &p_transformation_matrix,double &variance,bool id_s_transfered){

	cv::BruteForceMatcher<cv::L2<float>> matcher;
	vector<cv::DMatch> matches;
	matcher.match(descriptors[id_s],descriptors[id_t],matches);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
	keypoints_source->width    = 1*matches.size();
	keypoints_source->height   = 1;
	keypoints_source->is_dense = false;
	keypoints_source->points.resize (keypoints_source->width * keypoints_source->height);

	keypoints_target->width    = 1*matches.size();
	keypoints_target->height   = 1;
	keypoints_target->is_dense = false;
	keypoints_target->points.resize (keypoints_target->width * keypoints_target->height);

	
	pcl::CorrespondencesPtr correspondences (new pcl::Correspondences);
	correspondences->resize(matches.size());

	for(int i = 0;i < matches.size();i++){
		
		double x1 = keypoints[id_s][matches[i].queryIdx].pt.x;
		double y1 = keypoints[id_s][matches[i].queryIdx].pt.y;
		double x2 = keypoints[id_t][matches[i].trainIdx].pt.x;
		double y2 = keypoints[id_t][matches[i].trainIdx].pt.y;

		PointT searchPoint_s;
		if(id_s_transfered){
			searchPoint_s = transpointClouds[id_s]->points[(int)(y1+0.5)*pointClouds[id_s]->width+(int)(x1+0.5)];
		}else{
			searchPoint_s = pointClouds[id_s]->points[(int)(y1+0.5)*pointClouds[id_s]->width+(int)(x1+0.5)];
		}
		pcl::PointXYZRGB& pt_s = keypoints_source->points[i];
		pt_s.x = searchPoint_s.x;
		pt_s.y = searchPoint_s.y;
		pt_s.z = searchPoint_s.z;
		pt_s.r = searchPoint_s.r;
		pt_s.g = searchPoint_s.g;
		pt_s.b = searchPoint_s.b;


		PointT searchPoint_t = transpointClouds[id_t]->points[(int)(y2+0.5)*pointClouds[id_t]->width+(int)(x2+0.5)];
		pcl::PointXYZRGB& pt_t = keypoints_target->points[i];
		pt_t.x = searchPoint_t.x;
		pt_t.y = searchPoint_t.y;
		pt_t.z = searchPoint_t.z;
		pt_t.r = searchPoint_t.r;
		pt_t.g = searchPoint_t.g;
		pt_t.b = searchPoint_t.b;

		(*correspondences)[i].index_query = i;
		(*correspondences)[i].index_match = i;
		(*correspondences)[i].distance = 
		pcl::euclideanDistance(keypoints_source->at(i), keypoints_target->at(i));
	}

	double inlier_threshold = 0.005;
	int iteration_time = 1000;
	pcl::registration::CorrespondenceRejectorSampleConsensus<pcl::PointXYZRGB> rejector;
	rejector.setInputCloud(keypoints_source);
	rejector.setTargetCloud(keypoints_target);
	rejector.setInputCorrespondences(correspondences);
	rejector.setInlierThreshold(inlier_threshold);
	rejector.setMaxIterations(iteration_time);
	rejector.getCorrespondences(*ransac_correspondences);

	Eigen::Matrix4f trans = rejector.getBestTransformation();

	if(ransac_correspondences->size() < 3) return false;
	else{
		
		variance = 0;
		for(int i = 0;i < ransac_correspondences->size();i++){
			double dist = ransac_correspondences->at(i).distance;
			double weight =  ransac_correspondences->at(i).weight;
			variance += dist*dist*weight;
		}
		p_transformation_matrix << trans;
		return true;
	}
	
}

bool More_Plane3D(const Plane3D& a,const Plane3D& b)
{
	return a.inliers > b.inliers;
}

void dbplane_segment(pcl::PointCloud<pcl::Normal>::Ptr normals,PointCloud::Ptr cloud, double thre)
{
	int minInlier = 1000;
	vector<int> seg;
	vector<Plane3D> coffs;
	int csize = 0;
	for(int i = 0;i < cloud->size();i++)
	{
		seg.push_back(-1);
	}
	for(int i = 0;i < cloud->size();i++)
	{
		if(seg[i] == -1 && !isnan(normals->at(i).normal_x))
		{
			int id = csize;
			int count = 1;
			seg[i] = id;
			queue<int> search;
			search.push(i);

			//最小二乘拟合平面参数
			double Ex = 0,Ey = 0,Ez = 0,Ex2 = 0,Ey2 = 0,Ez2 = 0,Exy = 0,Eyz = 0,Exz = 0;

			while(!search.empty())
			{
				int p = search.front();
				search.pop();
				int y = p / cloud->width;
				int x = p - y*cloud->width;
				pcl::Normal np = normals->at(p);
				pcl::PointXYZRGB pp = cloud->at(p);

			   //参数更新
			   Ex += pp.x;
			   Ey += pp.y;
			   Ez += pp.z;
			   Ex2 += pp.x*pp.x;
			   Ey2 += pp.y*pp.y;
			   Ez2 += pp.z*pp.z;
			   Exy += pp.x*pp.y;
			   Eyz += pp.y*pp.z;
			   Exz += pp.x*pp.z;

				for(int jx = x-1;jx <= x+1;jx++)
				{
					if(jx < 0) continue;
					if(jx >= cloud->width) continue;
					for(int jy = y-1;jy <= y+1;jy++)
					{
						if(jy < 0) continue;
						if(jy >= cloud->height) continue;
						int nid = jy*cloud->width+jx;
						if(seg[nid] == -1 && !isnan(normals->at(nid).normal_x))
						{
							pcl::Normal nn = normals->at(nid);
							pcl::PointXYZRGB pn = cloud->at(nid);
							double disn = sqrt((np.normal_x-nn.normal_x)*(np.normal_x-nn.normal_x)+(np.normal_y-nn.normal_y)*(np.normal_y-nn.normal_y)+(np.normal_z-nn.normal_z)*(np.normal_z-nn.normal_z));
							double disp = sqrt((pp.x-pn.x)*(pp.x-pn.x)+(pp.y-pn.y)*(pp.y-pn.y)+(pp.z-pn.z)*(pp.z-pn.z));
							//double dis = sqrt(disn+disp);
							if(disn < thre && disp < thre)
							{
								search.push(nid);
								seg[nid] = id;
								count++;
							}
						}
					}
				}
			}
			double fem = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
			Plane3D coff;
			if(fem - 0 < 1e-10 && fem - 0 > -1e-10){
			  coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
			  coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
			  coff.c = 1;
			  coff.d = 0;
			}else{
			  coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
			  coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
			  coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
			  coff.d = 1;
			}
			double u = sqrt(coff.a*coff.a + coff.b*coff.b + coff.c*coff.c);
			coff.a /= u;
			coff.b /= u;
			coff.c /= u;
			coff.d /= u;

			if(count > minInlier)
			{
				//寻找平面的更多内点
				for(int j = i+1;j < cloud->size();j++)
				{
					if(seg[j] == -1 && !isnan(normals->at(j).normal_x))
		            {
						pcl::PointXYZRGB pj = cloud->at(j);
						if(coff.a * pj.x + coff.b * pj.y + coff.c * pj.z + coff.d < thre)
						{
							seg[j] = id;
							
							//参数更新
							Ex += pj.x;
							Ey += pj.y;
							Ez += pj.z;
							Ex2 += pj.x*pj.x;
							Ey2 += pj.y*pj.y;
							Ez2 += pj.z*pj.z;
							Exy += pj.x*pj.y;
							Eyz += pj.y*pj.z;
							Exz += pj.x*pj.z;
							count++;
						}
					}
				}
				//重新计算参数
				double fem2 = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
				if(fem2 - 0 < 1e-10 && fem2 - 0 > -1e-10){
				  coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
				  coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
				  coff.c = 1;
				  coff.d = 0;
				}else{
				  coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
				  coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
				  coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
				  coff.d = 1;
				}
				u = sqrt(coff.a*coff.a + coff.b*coff.b + coff.c*coff.c);
				coff.a /= u;
				coff.b /= u;
				coff.c /= u;
				coff.d /= u;
			}

			coff.inliers = count;
			coffs.push_back(coff);
			csize++;
		}
	}
	sort(coffs.begin(),coffs.end(),More_Plane3D);
	for(int i =0;i < 3;i++) cout<<"("<<coffs[i].a<<","<<coffs[i].b<<","<<coffs[i].c<<","<<coffs[i].d<<") "<<"inliers:"<<coffs[i].inliers<<endl;
	cout<<endl;
}

Plane3D coff_cloud(pcl::PointCloud<pcl::Normal>::Ptr normals,PointCloud::Ptr cloud)
{
	//最小二乘拟合平面参数
	double Ex = 0,Ey = 0,Ez = 0,Ex2 = 0,Ey2 = 0,Ez2 = 0,Exy = 0,Eyz = 0,Exz = 0;
	for(int i = 0;i < cloud->size();i++)
	{
		if(!isnan(normals->at(i).normal_x))
		{
			
			   pcl::PointXYZRGB pp = cloud->at(i);
			   //参数更新
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
	}
	double fem = -Ex2*Ey2*Ez2 + Ex2*Eyz*Eyz + Ez2*Exy*Exy + Ey2*Exz*Exz - 2*Exy*Eyz*Exz;
	Plane3D coff;
	if(fem - 0 < 1e-10 && fem - 0 > -1e-10){
		coff.a = (Exz*Ey2 - Exy*Eyz)/(Exy*Exy - Ex2*Ey2);
		coff.b = -(Exy*Exz - Ex2*Eyz)/(Exy*Exy - Ex2*Ey2);
		coff.c = 1;
		coff.d = 0;
	}else{
		coff.a = (Ex*Eyz*Eyz - Exz*Ey*Eyz - Ex*Ey2*Ez2 + Exy*Ey*Ez2 + Exz*Ey2*Ez - Exy*Eyz*Ez)/fem;
		coff.b = (Exz*Exz*Ey - Ex*Exz*Eyz + Ex*Exy*Ez2 - Exy*Exz*Ez - Ex2*Ey*Ez2 + Ex2*Eyz*Ez)/fem;
		coff.c = (Exy*Exy*Ez + Ex*Exz*Ey2 - Ex*Exy*Eyz - Exy*Exz*Ey + Ex2*Ey*Eyz - Ex2*Ey2*Ez)/fem;
		coff.d = 1;
	}
	double u = sqrt(coff.a*coff.a + coff.b*coff.b + coff.c*coff.c);
	coff.a /= u;
	coff.b /= u;
	coff.c /= u;
	coff.d /= u;

	//cout<<"("<<coff.a<<","<<coff.b<<","<<coff.c<<","<<coff.d<<") "<<"inliers:"<<coff.inliers<<endl;
	return coff;
}

void do_SLAM(string sFilename,int argc,char** argv){
	/*
	vector<vector<cv::KeyPoint>> keypoints;
	vector<cv::Mat> descriptors;
	vector<PointCloud::Ptr> pointClouds;
	vector<PointCloud::Ptr> transpointClouds;
	vector<Eigen::Matrix4f*> trans;
	*/


	PlaneManager p_manager = *(new PlaneManager());
	

	XnStatus result = XN_STATUS_OK;    
	xn::DepthMetaData depthMD;  
	xn::ImageMetaData imageMD;  

	int x=1;

	//OpenCV  
	/*
	IplImage*  imgDepth16u=cvCreateImage(cvSize(640,480),IPL_DEPTH_16U,1); 
	IplImage* imgRGB8u=cvCreateImage(cvSize(640,480),IPL_DEPTH_8U,3);
	IplImage*  depthShow=cvCreateImage(cvSize(640,480),IPL_DEPTH_8U,1);
	IplImage* imageShow=cvCreateImage(cvSize(640,480),IPL_DEPTH_8U,3);
	*/
	cv::Mat  imgDepth16u( 480,640,CV_16UC1);
    cv::Mat  imgRGB8u( 480,640,CV_8UC3);
    cv::Mat  depthShow( 480,640,CV_8UC1);
    cv::Mat  imageShow( 480,640,CV_8UC3);
	cvNamedWindow("depth",1);  
	cvNamedWindow("image",1); 
	char key=0;  

	// context   
	xn::Context context;   
	result = context.Init();   
	CheckOpenNIError( result, "initialize context" );    

	result = context.OpenFileRecording( sFilename.c_str() );
	CheckOpenNIError(result, "Open record");

	// creategenerator    
	xn::DepthGenerator depthGenerator;    
	result = depthGenerator.Create( context );   
	CheckOpenNIError( result, "Create depth generator" );    
	xn::ImageGenerator imageGenerator;  
	result = imageGenerator.Create( context );
	CheckOpenNIError( result, "Create image generator" );  
 
	//map mode    
	XnMapOutputMode mapMode;   
	mapMode.nXRes = 640;    
	mapMode.nYRes = 480;   
	mapMode.nFPS = 30;   
	result = depthGenerator.SetMapOutputMode( mapMode );    
	result = imageGenerator.SetMapOutputMode( mapMode );    
  
	// correct view port    
	depthGenerator.GetAlternativeViewPointCap().SetViewPoint( imageGenerator );   
 
	//read data  
	result = context.StartGeneratingAll();      
	result = context.WaitNoneUpdateAll();    

	XnFieldOfView fov;
	depthGenerator.GetFieldOfView(fov);

	//Set up Globals
	XtoZ = atan(fov.fHFOV/2)*2;
	YtoZ = atan(fov.fVFOV/2)*2;
	Eigen::Matrix4f zero;  
	Eigen::Matrix4f lzero; 
	Eigen::Matrix4f dzero; 
	zero << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;

	lzero<< 0, 1, 0, 0,  
			1, 0, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;

	dzero<< 0, 0, -1, 0,  
			0, 1, 0, 0,  
			-1, 0, 0, 0,  
			0, 0, 0, 1;

	//Sdf::sdf = new Sdf(zero, 1, 1.2, 4.0, -2.0, -3.5, -1.0, 0.01, 0.1, 0.001); //xy
	//Sdf::sdf = new Sdf(zero, 4.0, 1, 1.2, -1.0, -2.0, -3.5, 0.01, 0.1, 0.001); //xz

	Sdf::sdf = NULL;//new Sdf(zero, 1.2, 1, 4.0, -3.5, -2.0, -1.0, 0.01, 0.1, 0.001);
	//Sdf::tsdf = new TsdfNode();
	//Sdf::tsdf->Initial(-3.0,-2.0,0.0,1.2,1.0,4.0,0.05);
	//Sdf::tsdf->split();
	clock_t start, finish;

	float dthresh = 0.5;

	int k = 0,d = 0;
	while( (key!=27) && !(result = context.WaitAndUpdateAll( ))  && d < 440)
	{    
		//get meta data  
		depthGenerator.GetMetaData(depthMD);   
		imageGenerator.GetMetaData(imageMD);
		d++;
		//if( d < 150) continue;
		if( d % 10 != 1) continue;
 
		//OpenCV output  
		memcpy(imgDepth16u.data,depthMD.Data(),640*480*2);  
		imgDepth16u.convertTo(depthShow,CV_8U,255/2096.0);
		//cvConvertScale(imgDepth16u,depthShow,255/4096.0,0);  
		memcpy(imgRGB8u.data,imageMD.Data(),640*480*3); 
		cvtColor(imgRGB8u,imageShow,CV_RGB2BGR);  
		imshow("depth", depthShow);  
		imshow("image",imageShow); 

		
		cv::SurfFeatureDetector detector;
		vector<cv::KeyPoint> keypoint;
		detector.detect(imageShow,keypoint);
		cv::SurfDescriptorExtractor extractor;
		cv::Mat descriptor;
		extractor.compute(imageShow,keypoint,descriptor);
		keypoints.push_back(keypoint);
		descriptors.push_back(descriptor);

		
		XnUInt64 timestamp = depthMD.Timestamp();
        XnInt32 frameID = depthMD.FrameID();
		PointCloud::Ptr cloud = convertToPointCloud(depthGenerator,depthMD.Data(),imageMD.RGB24Data(),timestamp,frameID);
		pointClouds.push_back(cloud);

		foundPlanes.push_back(*(new vector<int>()));


		//front-end 用RANSAC方法做第一次粗匹配
		if(k == 0){
			transpointClouds.push_back(cloud);
			Eigen::Matrix4f tran = Eigen::Matrix4f::Identity();
			trans.push_back(tran);
		}
		else{
			XnFieldOfView fov;
			depthGenerator.GetFieldOfView(fov);
			XnMapOutputMode outputMode;
			depthGenerator.GetMapOutputMode(outputMode);

			double fXToZ = tan(fov.fHFOV/2)*2;	
			double fYToZ = tan(fov.fVFOV/2)*2;

			XnUInt32 nHalfXres = outputMode.nXRes / 2;	
			XnUInt32 nHalfYres = outputMode.nYRes / 2;

			XnDouble fCoeffX = outputMode.nXRes / fXToZ;
			XnDouble fCoeffY = outputMode.nYRes / fYToZ;

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
	        pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
			double variance = 0;
			pcl::CorrespondencesPtr ransac_correspondences (new pcl::Correspondences);
			Eigen::Matrix4f init;
			do_match(k,k-1,ransac_correspondences,init,variance,false);

			//Eigen::Matrix4f init2 = trans[k - 1];
			Eigen::Matrix4f tran = cpuICP(transpointNormals[k - 1], transpointClouds[k - 1], cloud, trans[k - 1].inverse(), init, nHalfXres, nHalfYres, fCoeffX, fCoeffY,dthresh);

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*(pointClouds.at(k)), *transCloud, tran);
			transpointClouds.push_back(transCloud);
			trans.push_back(tran);
		}

		pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
		pcl::IntegralImageNormalEstimation<PointT, pcl::Normal> ne;
	    ne.setNormalEstimationMethod (ne.AVERAGE_3D_GRADIENT);
	    ne.setMaxDepthChangeFactor(0.03f);
	    ne.setNormalSmoothingSize(20.0f);
	    ne.setInputCloud(transpointClouds[k]);
	    ne.compute(*normals);
		transpointNormals.push_back(normals);

		/*
		if(k > 0){
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
	        pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
			double variance = 0;
			try_do_match(k,k-1,keypoints_source,keypoints_target,variance,&pointClouds,&transpointClouds,k-1);
			
			vector<int> CandidateFrame;
			for(int i = 0, j = k-2;i < 3 && j >= 0;i++,j--){
				CandidateFrame.push_back(j);
			}
			//for(int i = (k/10);i >= 0;i--){
			//	if(k - i*10 > 3){
			//		CandidateFrame.push_back(i*10);
			//	}
			//}
			int count = 0;
			for(int i = 0;i < CandidateFrame.size();i++){
				pcl::PointCloud<pcl::PointXYZRGB>::Ptr add_source(new pcl::PointCloud<pcl::PointXYZRGB>);
	            pcl::PointCloud<pcl::PointXYZRGB>::Ptr add_target(new pcl::PointCloud<pcl::PointXYZRGB>);
				bool match_flag = try_do_match(k,CandidateFrame[i],add_source,add_target,variance,&pointClouds,&transpointClouds,CandidateFrame[i]);
				if(match_flag && variance < 0.1 && variance > 0.00001){
					for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = add_source->begin();it != add_source->end();it++){
						keypoints_source->points.push_back(*it);
		            } 
					for(pcl::PointCloud<pcl::PointXYZRGB>::VectorType::iterator it = add_target->begin();it != add_target->end();it++){
						keypoints_target->points.push_back(*it);
		            } 
				}
			}
			
			Eigen::Matrix4f tran;
			bool match_flag = real_do_match(keypoints_source,keypoints_target,tran);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*(pointClouds.at(k)), *transCloud, tran);
			transpointClouds.push_back(transCloud);
			trans.push_back(tran);
		}
		*/
		k++;

		key=cvWaitKey(10); 
	}
	
   //manager.initializeOptimization();
   boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1 (new pcl::visualization::PCLVisualizer ("Viewer_1"));

	vis1->initCameraParameters ();
	vis1->setBackgroundColor (0.1, 0.1, 0.1);
	vis1->addText("粗配准", 10, 10, "粗配准v1");

	for(int i = 0;i < pointClouds.size();i++){

		//pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
		//double* est = new double[7];
		//manager.optimizer->vertex(i)->getEstimateData(est);
		//Eigen::Matrix4f m = *(manager.getEuleMatrix(est));
		//pcl::transformPointCloud(*(transpointClouds.at(i)), *transCloud, *(manager.getEuleMatrix(est)));

		pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb_cloud(transpointClouds[i]);
		stringstream ss;
		ss<<"cloud"<<i;
		vis1->addPointCloud<pcl::PointXYZRGB> (transpointClouds[i], rgb_cloud, ss.str());
		vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
	}


	while (!vis1->wasStopped())
	{
		//在此处可以添加其他处理
		vis1->spinOnce (100);
		//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}

	//GL_Mesh::gl_vertexes = new std::vector<GLVertex>();
	//GL_Mesh::gl_meshes = new std::map<std::tuple<int,int,int>,GLMesh>();
	//GL_Mesh::gl_edges = new std::map<std::pair<int,int>,int>();
	//GL_Mesh::gl_trans = new std::vector<Eigen::Matrix4f>();
	int lp =0;
	Plane3D current_coff;
	for(vector<PointCloud::Ptr>::iterator it = transpointClouds.begin();it != transpointClouds.end();it++,lp++){
		//求平面
		pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
		pcl::IntegralImageNormalEstimation<PointT, pcl::Normal> ne;
	    ne.setNormalEstimationMethod (ne.AVERAGE_3D_GRADIENT);
	    ne.setMaxDepthChangeFactor(0.03f);
	    ne.setNormalSmoothingSize(20.0f);
	    ne.setInputCloud(*it);
	    ne.compute(*normals);

		Plane3D coff = coff_cloud(normals,*it);
		if(!Sdf::sdf)
		{
			current_coff = coff;
			//dbplane_segment(normals,*it,0.02);
			double before[4];
			double after[4];

			before[0] = coff.a;
			before[1] = coff.b;
			before[2] = coff.c;
			before[3] = coff.d;

			after[0] = 0;
			after[1] = 0;
			after[2] = 1;
			after[3] = 0;

			Eigen::Matrix4f local = CalTransport(before,after);
			//Eigen::Matrix4f local = zero;
			Eigen::Matrix4f local_c = local.inverse();

			bool flag = false;
			double maxX,maxY,maxZ,minX,minY,minZ;
			PointCloud::Ptr cloud = *it;
			for(int i = 0;i < cloud->size();i++)
			{
				pcl::PointXYZRGB pp = cloud->at(i);
				Eigen::Vector4f v;
				v[0] = pp.x ;
				v[1] = pp.y;
				v[2] = pp.z;
				v[3] = 1;
				Eigen::Vector4f t = local_c * v;
				if(!flag)
				{
					maxX = t[0];
					maxY = t[1];
					maxZ = t[2];
					minX = t[0];
					minY = t[1];
					minZ = t[2];
					flag = true;
				}else{
					if(t[0] > maxX) maxX = t[0];
					if(t[1] > maxY) maxY = t[1];
					if(t[2] > maxZ) maxZ = t[2];
					if(t[0] < minX) minX = t[0];
					if(t[1] < minY) minY = t[1];
					if(t[2] < minZ) minZ = t[2];
				}
			}
		
		
			//cout<<maxX-minX<<","<<maxY-minY<<","<<maxZ-minZ<<endl;
			double dev = 0.01;
			//Sdf::sdf = new Sdf(zero, 1.2, 1, 4.0, -3.5, -2.5, -3.0, 0.01, 0.1, 0.001);
			Sdf::sdf = new Sdf(local, maxX, maxY, maxZ, minX, minY, minZ, dev, dev*10, dev/10);
			//Sdf::sdf = new Sdf(local, 1.2, 1, 1.0, -1.5, -1.0, 0.0, dev, dev*10, dev/10);
		}else{
			double dis = (coff.a*current_coff.a)+(coff.b*current_coff.b)+(coff.c*current_coff.c);
			//cout<<"cur_coff"<<current_coff.a<<","<<dis<<endl;
			if(dis < 0.7) current_coff = coff;
			
			
	    	current_coff = coff;
			//dbplane_segment(normals,*it,0.02);
			double before[4];
			double after[4];

			before[0] = coff.a;
			before[1] = coff.b;
			before[2] = coff.c;
			before[3] = coff.d;

			after[0] = 0;
			after[1] = 0;
			after[2] = 1;
			after[3] = 0;

			//Eigen::Matrix4f local = CalTransport(before,after);
			Eigen::Matrix4f local = zero;
			Eigen::Matrix4f local_c = local.inverse();

			bool flag = false;
			double maxX,maxY,maxZ,minX,minY,minZ;
			PointCloud::Ptr cloud = *it;
			for(int i = 0;i < cloud->size();i++)
			{
				pcl::PointXYZRGB pp = cloud->at(i);
				Eigen::Vector4f v;
				v[0] = pp.x ;
				v[1] = pp.y;
				v[2] = pp.z;
				v[3] = 1;
				Eigen::Vector4f t = local_c * v;
				if(!flag)
				{
					maxX = t[0];
					maxY = t[1];
					maxZ = t[2];
					minX = t[0];
					minY = t[1];
					minZ = t[2];
					flag = true;
				}else{
					if(t[0] > maxX) maxX = t[0];
					if(t[1] > maxY) maxY = t[1];
					if(t[2] > maxZ) maxZ = t[2];
					if(t[0] < minX) minX = t[0];
					if(t[1] < minY) minY = t[1];
					if(t[2] < minZ) minZ = t[2];
				}
			}
		
		
			//cout<<maxX-minX<<","<<maxY-minY<<","<<maxZ-minZ<<endl;
			double dx = (maxX - minX);
			double dy = (maxY - minY);
			double dz = (maxZ - minZ);
			cout<<"size:"<<dx<<","<<dy<<","<<dz<<endl;
			//Sdf::sdf->changeLocal(local, 1.2, 1, 4.0, -3.5, -2.0, -1.0);
			//start = clock(); 
			Sdf::sdf->changeLocal(local, maxX, maxY, maxZ, minX, minY, minZ);
			//finish = clock(); 
			//double duration = (double)(finish - start) / CLOCKS_PER_SEC; 
			//cout<<duration<<" ";
		}

		DepthCalculator dg;
		dg.cloud = (*it);
		dg.depthGenerator = &depthGenerator;
		dg.normals = normals;
		dg.tran = trans[lp];

		//Sdf::tsdf->Traversal(&dg);
		//GL_Mesh::gl_meshes = NULL;
		//GLWindows();
		start = clock();
		Sdf::sdf->Traversal2(&dg);
		finish = clock(); 
		double duration = (double)(finish - start) / CLOCKS_PER_SEC; 
		cout<<duration<<" ";
		//if(lp >= 0)
		//{
			GL_Mesh::gl_trans->push_back(trans[lp]);
			start = clock(); 
			Sdf::sdf->castNextMeshes(GL_Mesh::gl_trans->size()-1);
			finish = clock(); 
			duration = (double)(finish - start) / CLOCKS_PER_SEC; 
			cout<<duration<<endl;
			//GL_Mesh::gl_trans->push_back(dzero * trans[lp]);
		    //Sdf::sdf->castNextMeshes(GL_Mesh::gl_trans->size()-1);
			//GL_Mesh::gl_trans->push_back(lzero * trans[lp]);
		    //Sdf::sdf->castNextMeshes(GL_Mesh::gl_trans->size()-1);
		//}


		
		/*
		pcl::OrganizedMultiPlaneSegmentation<PointT,pcl::Normal,pcl::Label>  mps;
		mps.setMinInliers(1000);
		mps.setAngularThreshold(0.017453*2.0);
		mps.setDistanceThreshold(0.02);
		mps.setInputNormals(normals);
		mps.setInputCloud(*it);
		std::vector<pcl::PlanarRegion<PointT>, Eigen::aligned_allocator<pcl::PlanarRegion<PointT> > > regions;
		mps.segmentAndRefine(regions);

		for(int j=0;j < regions.size();j++){
		   Eigen::Vector4f cof = regions[j].getCoefficients();
		   int addplane = p_manager.addPlane(cof[0],cof[1],cof[2],cof[3],regions[j].getCount(),regions[j].getContour());
		   foundPlanes[lp].push_back(addplane);
	    }
		*/
	}
	//Sdf::sdf->getData();
	//GLWindows();
	
	//GL_Mesh::gl_trans->push_back(trans[0]);
	//GL_Mesh::gl_trans->push_back(trans[5]);
	//GL_Mesh::gl_trans->push_back(trans[10]);
	//GL_Mesh::gl_trans->push_back(trans[15]);
	//GL_Mesh::gl_trans->push_back(trans[20]);
	
	//for(int i = 0;i < GL_Mesh::gl_trans->size();i++)
	//{
	//	Sdf::sdf->castNextMeshes(i);
	//}
	//GL_Mesh::gl_trans->push_back(trans[15]);
	//GL_Mesh::gl_trans->push_back(trans[20]);
	//GL_Mesh::gl_trans->push_back(trans[25]);
	//GL_Mesh::gl_trans->push_back(trans[30]);
	//Sdf::sdf->castIsoMeshes();

	//GL_Mesh::gl_vertex_map = new std::map<std::tuple<int,int,int>,GLMergeVertex>();
	//Sdf::sdf->castMeshes(trans[0]);
	//Sdf::sdf->castMeshes(trans[5]);
	//Sdf::sdf->castMeshes(trans[10]);
	Sdf::sdf->outputObj();
	//GLWindows();

	Sdf::sdf->freeData();

	/*
	p_manager.clusterPlane();
	p_manager.getCoordinary();
	p_manager.printPlaneInfo();
	p_manager.generateSkeleton();

	for(int k =0;k < p_manager.clutered_planes->size();k++){
		if(p_manager.clutered_planes->at(k).skeletonGenerated && p_manager.clutered_planes->at(k).boundPoints != NULL){
			stringstream ss;
			ss<<"cloud"<<k;
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr planeBounds (new pcl::PointCloud<pcl::PointXYZRGB>);
			planeBounds->width    = 4;
			planeBounds->height   = 1;
			planeBounds->is_dense = false;
			planeBounds->points.resize (planeBounds->width * planeBounds->height);
			for(int l = 0;l < 4;l++){
				pcl::PointXYZRGB& pt_s = planeBounds->points[l];
				pt_s.x = p_manager.clutered_planes->at(k).boundPoints[l].x;
				pt_s.y = p_manager.clutered_planes->at(k).boundPoints[l].y;
				pt_s.z = p_manager.clutered_planes->at(k).boundPoints[l].z;
				pt_s.r = 1;
				pt_s.g = 0;
				pt_s.b = 0;
			     
			}
			vis1->addPolygon<pcl::PointXYZRGB>(planeBounds,GETCOLORR(k),GETCOLORG(k),GETCOLORB(k),ss.str());
		}
    }
	vis1->spinOnce (100);
	while (!vis1->wasStopped())
	{
		//在此处可以添加其他处理
		vis1->spinOnce (100);
		//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}

	int framesize = transpointClouds.size();
	int planesize = p_manager.clutered_planes->size();
	int count = 0;
	vector<frameMean>** planeSets = new vector<frameMean>*[planesize];
	for(int i = 0;i < planesize;i++){
		planeSets[i] = new vector<frameMean>();
	}

	lp =0;
	double tr = 0.03;
	int mininlier = 20;
	//while(true){
	//cin>>tr>>mininlier;
    int* countP = new int[planesize];
	double* sx = new double[planesize];
	double* sy = new double[planesize];
	double* sz = new double[planesize];
	for(vector<PointCloud::Ptr>::iterator it = transpointClouds.begin();it != transpointClouds.end();it++,lp++){
		 for(int i = 0;i < planesize;i++){
			countP[i] = 0;
			sx[i] = 0.0;
			sy[i] = 0.0;
			sz[i] = 0.0;
	     }
		 PointCloud::Ptr cloud = *it;
		 for(PointCloud::iterator it2 = cloud->begin();it2 != cloud->end();it2++){
			 PointT p = *it2;
			 int cc = 0;
			 for(vector<plane_info>::iterator it3 = p_manager.clutered_planes->begin();it3 != p_manager.clutered_planes->end();it3++,cc++){
				 plane_info plane = *it3;
				 double dis = p.x *  plane.a + p.y *  plane.b + p.z *  plane.c - plane.d;
				 if(dis < tr && dis > -tr){
					 sx[cc] += p.x;
					 sy[cc] += p.y;
					 sz[cc] += p.z;
					 countP[cc]++;
				 }
	         }
		 }
		  for(int i = 0;i < planesize;i++){
			  if(countP[i] > cloud->size() / mininlier){
				  frameMean* fa = new frameMean();
			      fa->frame = lp;
				  fa->m[0] = sx[i]/countP[i];
				  fa->m[1] = sy[i]/countP[i];
				  fa->m[2] = sz[i]/countP[i];
				  planeSets[i]->push_back(*fa);
			  }
	     }
	}
	
	for(int i = 0;i < planesize;i++){
		if(planeSets[i]->size() > 0){
			double minx,miny,minz,maxx,maxy,maxz;
			vector<frameMean>::iterator it = planeSets[i]->begin();
			minx = maxx = (*it).m[0];
			miny = maxy = (*it).m[1];
			minz = maxz = (*it).m[2];
			it++;
			for(;it != planeSets[i]->end();it++){
				if((*it).m[0] < minx) minx = (*it).m[0];
				if((*it).m[0] > maxx) maxx = (*it).m[0];
				if((*it).m[1] < miny) miny = (*it).m[1];
				if((*it).m[1] > maxy) maxy = (*it).m[1];
				if((*it).m[2] < minz) minz = (*it).m[2];
				if((*it).m[2] > maxz) maxz = (*it).m[2];
			}
			double xsize = maxx-minx;
			double ysize = maxy-miny;
			double zsize = maxz-minz;
			if(xsize >= ysize && xsize >= zsize){
				sort(planeSets[i]->begin(), planeSets[i]->end(), lessmark0);
			}else if(ysize >= xsize && ysize >= zsize){
				sort(planeSets[i]->begin(), planeSets[i]->end(), lessmark1);
			}else{
				sort(planeSets[i]->begin(), planeSets[i]->end(), lessmark2);
			}
		}
	}

	for(int i = 0;i < planesize;i++){
		cout<<"plane"<<i<<":";
		for(vector<frameMean>::iterator it = planeSets[i]->begin();it != planeSets[i]->end();it++){
			cout<<(*it).frame<<",";
		}
		string s = p_manager.clutered_planes->at(i).skeletonGenerated?"true":"false";
		cout<<s<<endl;
	}

	vis1->spinOnce (100);
	while (!vis1->wasStopped())
	{
		//在此处可以添加其他处理
		vis1->spinOnce (100);
		//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}


	int planewall[] = {2,3,4,5};
	int wallsize = 4;
	vector<Eigen::Matrix4f*>** Walltrans = new vector<Eigen::Matrix4f*>*[wallsize];
	for(int j = 0;j < wallsize;j++){
		Walltrans[j] = new vector<Eigen::Matrix4f*>();
		vector<PointCloud::Ptr>* newtranspointClouds = new vector<PointCloud::Ptr>();
		int pssize = planeSets[planewall[j]]->size();
		int* ps = new int[pssize];
		ps[0] = planeSets[planewall[j]]->at(0).frame;
		newtranspointClouds->push_back(transpointClouds[ps[0]]);
		Walltrans[j]->push_back(&zero);
		for(int i =1;i < pssize;i++){
			ps[i] = planeSets[planewall[j]]->at(i).frame;
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
			double variance = 0;
			bool match_flag = try_do_match(ps[i],ps[i-1],keypoints_source,keypoints_target,variance,&transpointClouds,newtranspointClouds,i-1);
			Eigen::Matrix4f* tran = new Eigen::Matrix4f();
			match_flag = real_do_match(keypoints_source,keypoints_target,*tran);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*(transpointClouds.at(ps[i])), *transCloud, *tran);
			newtranspointClouds->push_back(transCloud);
			Walltrans[j]->push_back(tran);
		}
	}

   boost::shared_ptr<pcl::visualization::PCLVisualizer> vis2 (new pcl::visualization::PCLVisualizer ("Viewer_2"));

	vis2->initCameraParameters ();
	vis2->setBackgroundColor (0.1, 0.1, 0.1);
	vis2->addText("粗配准", 10, 10, "粗配准v1");

	for(int i = 0;i <  planeSets[planewall[0]]->size();i++){
		int id = planeSets[planewall[0]]->at(i).frame;

		pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
		Eigen::Matrix4f* m = Walltrans[0]->at(i);
		pcl::transformPointCloud(*(transpointClouds.at(id)), *transCloud, *m);

		pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb_cloud(transCloud);
		stringstream ss;
		ss<<"cloud"<<i;
		vis2->addPointCloud<pcl::PointXYZRGB> (transCloud, rgb_cloud, ss.str());
		vis2->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
	}

   	while (!vis2->wasStopped())
	{
		//在此处可以添加其他处理
		vis2->spinOnce (100);
		//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
	}
	*/
}

int
 main2 (int argc, char** argv)
{
  //do_SLAM("F:/1.ONI",argc,argv);	
  do_SLAM("D:/PCL/record.ONI",argc,argv);
  //do_SLAM("D:/PCL/DATA/16.ONI",argc,argv);
  return 0;
}
 


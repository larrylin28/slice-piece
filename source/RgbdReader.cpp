#include "RGBDReader.h"
#include "SdfModel.h"
#include "Segmentation.h"
#include "icp.h"

//pcl
#include <pcl/registration/correspondence_estimation.h>
#include <pcl/registration/correspondence_rejection_sample_consensus.h>
#include <pcl/registration/ia_ransac.h>
#include <pcl/features/integral_image_normal.h>
#include <pcl/visualization/pcl_visualizer.h>

//openNI
#include <XnCppWrapper.h>

//openCV
#include <opencv/highgui.h>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/legacy/legacy.hpp>

#include <iostream>

#define GETCOLORR(i) (i%4)*1.0/3
#define GETCOLORG(i) ((i+1)%5)*1.0/4
#define GETCOLORB(i) (((1000-i)+(i/4))%4)*1.0/3

using namespace std;

namespace Rgbd
{
	cv::Mat  imgDepth16u( 480,640,CV_16UC1);
	cv::Mat  imgRGB8u( 480,640,CV_8UC3);
	cv::Mat  depthShow( 480,640,CV_8UC1);
	cv::Mat  imageShow( 480,640,CV_8UC3);

	RgbdDatas::RgbdDatas():context(new xn::Context),
		                  depthMD(new xn::DepthMetaData), 
						  imageMD(new xn::ImageMetaData),
						  depthGenerator(new xn::DepthGenerator),
						  imageGenerator(new xn::ImageGenerator)
	{
	}
	
	RgbdDatas::~RgbdDatas()
	{
	}

	bool RgbdDatas::get_match(int id_s,int id_t,Eigen::Matrix4f &p_transformation_matrix,double &variance,bool id_s_transfered){
		pcl::CorrespondencesPtr ransac_correspondences (new pcl::Correspondences);
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

	//convert rgbd data to point cloud
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr RGBD2PointCloud(xn::DepthGenerator& rDepthGen,const XnDepthPixel *dm,const XnRGB24Pixel *im,XnUInt64 timestamp, XnInt32 frameID)		
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

	inline void CheckError_( XnStatus result, string status )  
	{   
		if( result != XN_STATUS_OK )   
			cerr << status << " Error: " << xnGetStatusString( result ) << endl;  
	} 

	RgbdReader::RgbdReader():datas(new RgbdDatas())
	{
	}

	RgbdReader::~RgbdReader()
	{
	}

	void RgbdReader::Read(string sFilename, int startFrame, int skip, int frameCount, double dthresh)
	{
		initReading(sFilename);
		int readed_frames = 0;
		int real_frames = 0;
		char key = 0;  
		XnStatus result = XN_STATUS_OK;

		XnFieldOfView fov;
		datas->depthGenerator->GetFieldOfView(fov);
		XnMapOutputMode outputMode;
		datas->depthGenerator->GetMapOutputMode(outputMode);

		double fXToZ = tan(fov.fHFOV/2)*2;	
		double fYToZ = tan(fov.fVFOV/2)*2;

		nHalfXres = outputMode.nXRes / 2;	
		nHalfYres = outputMode.nYRes / 2;

		fCoeffX = outputMode.nXRes / fXToZ;
		fCoeffY = outputMode.nYRes / fYToZ;

		while( (key!=27) && !(result = datas->context->WaitAndUpdateAll( ))  && readed_frames < frameCount)
		{  
			//get meta data  
			datas->depthGenerator->GetMetaData(*(datas->depthMD));   
			datas->imageGenerator->GetMetaData(*(datas->imageMD));

			//skip some frames
			real_frames++;
			if( real_frames <= startFrame) continue;
		    if( real_frames % skip != 1) continue;
			cout<<"record_frame:"<<real_frames<<endl;
			
			//tranform to point cloud
			readFrame();
			//do slam for point cloud
			doSlam(readed_frames, dthresh);

			readed_frames++;
		    key=cvWaitKey(10);
		}
	}

	void RgbdReader::ShowPointCoud()
	{
		boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1 (new pcl::visualization::PCLVisualizer ("Viewer_1"));

		vis1->initCameraParameters ();
		vis1->setBackgroundColor (0.1, 0.1, 0.1);
		vis1->addText("point clouds", 10, 10, "point clouds");

		for(int i = 0;i < datas->transpointClouds.size();i++){

			//pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			//double* est = new double[7];
			//manager.optimizer->vertex(i)->getEstimateData(est);
			//Eigen::Matrix4f m = *(manager.getEuleMatrix(est));
			//pcl::transformPointCloud(*(transpointClouds.at(i)), *transCloud, *(manager.getEuleMatrix(est)));

			pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb_cloud(datas->transpointClouds[i]);
			stringstream ss;
			ss<<"cloud"<<i;
			vis1->addPointCloud<pcl::PointXYZRGB> (datas->transpointClouds[i], rgb_cloud, ss.str());
			vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
		}


		while (!vis1->wasStopped())
		{
			//在此处可以添加其他处理
			vis1->spinOnce (100);
			//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
		}
	}


	bool RgbdReader::initReading(string sFilename)
	{
		XnStatus result = XN_STATUS_OK;      

	
		cvNamedWindow("depth",1);  
		cvNamedWindow("image",1); 
		char key=0;  

		// context    
		result = datas->context->Init();   
		CheckError_( result, "initialize context" );    

		result = datas->context->OpenFileRecording( sFilename.c_str() );
		CheckError_(result, "Open record");
  
		result = datas->depthGenerator->Create( *(datas->context) );   
		CheckError_( result, "Create depth generator" );    
		xn::ImageGenerator imageGenerator;  
		result = datas->imageGenerator->Create( *(datas->context) );
		CheckError_( result, "Create image generator" );  
 
		//map mode    
		XnMapOutputMode mapMode;   
		mapMode.nXRes = 640;    
		mapMode.nYRes = 480;   
		mapMode.nFPS = 30;   
		result = datas->depthGenerator->SetMapOutputMode( mapMode );    
		result = datas->imageGenerator->SetMapOutputMode( mapMode );    
  
		// correct view port    
		datas->depthGenerator->GetAlternativeViewPointCap().SetViewPoint( (*datas->imageGenerator) );   
 
		//read data  
		result = datas->context->StartGeneratingAll();      
		result = datas->context->WaitNoneUpdateAll();    

		XnFieldOfView fov;
		datas->depthGenerator->GetFieldOfView(fov);

		//Set up parameters
		datas->XtoZ = atan(fov.fHFOV/2)*2;
		datas->YtoZ = atan(fov.fVFOV/2)*2;

		return (result == XN_STATUS_OK);   
	}


	void RgbdReader::readFrame()
	{
		//OpenCV output  
		memcpy(imgDepth16u.data,datas->depthMD->Data(),640*480*2);  
		imgDepth16u.convertTo(depthShow,CV_8U,255/2096.0);
		//cvConvertScale(imgDepth16u,depthShow,255/4096.0,0);  
		memcpy(imgRGB8u.data,datas->imageMD->Data(),640*480*3); 
		cvtColor(imgRGB8u,imageShow,CV_RGB2BGR);  
		imshow("depth", depthShow);  
		imshow("image",imageShow); 

		
		cv::SurfFeatureDetector detector;
		vector<cv::KeyPoint> keypoint;
		detector.detect(imageShow,keypoint);
		cv::SurfDescriptorExtractor extractor;
		cv::Mat descriptor;
		extractor.compute(imageShow,keypoint,descriptor);
		datas->keypoints.push_back(keypoint);
		datas->descriptors.push_back(descriptor);

		
		XnUInt64 timestamp = datas->depthMD->Timestamp();
		XnInt32 frameID = datas->depthMD->FrameID();
		PointCloud::Ptr cloud = RGBD2PointCloud(*datas->depthGenerator,datas->depthMD->Data(),datas->imageMD->RGB24Data(),timestamp,frameID);
		datas->pointClouds.push_back(cloud);

	}

	int RgbdReader::getFrameCount()
	{
		return datas->transpointClouds.size();
	}

	Eigen::Matrix4f RgbdReader::getTrans(int frameId)
	{
		return datas->trans[frameId];
	}

	void RgbdReader::putDataToSdf(SdfModel* sdf, int frameId)
	{
		sdf->dataFusion(datas->transpointClouds[frameId], datas->transpointNormals[frameId], datas->trans[frameId], nHalfXres, nHalfYres, fCoeffX , fCoeffY, this);
	}

	void RgbdReader::initaltags(int frameId, std::vector<int>& oldTags, std::vector<int>& newTags, double dthresh)
	{
		initalTags(datas->transpointClouds[frameId - 1], datas->transpointClouds[frameId], datas->trans[frameId - 1].inverse(), datas->trans[frameId],  nHalfXres, nHalfYres, fCoeffX, fCoeffY, dthresh, &(oldTags[0]), &(newTags[0]));
	}

	Eigen::Matrix4f RgbdReader::getExtraSlam(int frameId, double dthresh, int extraTag, vector<int>& oldTags, vector<int>& newTags, Eigen::Matrix4f& extraTr)
	{
		Eigen::Matrix4f init = Eigen::Matrix4f::Identity();
		Eigen::Matrix4f tran = cpuExtraICP(datas->transpointNormals[frameId - 1], datas->transpointClouds[frameId - 1], datas->transpointClouds[frameId], datas->trans[frameId - 1].inverse(), init, nHalfXres, nHalfYres, fCoeffX, fCoeffY,dthresh, extraTag, &(oldTags[0]), &(newTags[0]), extraTr);

		Eigen::Matrix3f rot   = tran.block<3, 3> (0, 0);
        Eigen::Vector3f trans = tran.block<3, 1> (0, 3);
		Eigen::Vector3f line;
		line << 1, 1, 1;
		Eigen::Vector3f rline = rot* line + trans - line;
		double dis = sqrt(rline[0] * rline[0] + rline[1] * rline[1] + rline[2] * rline[2]);
		cout<<dis;
		if(dis > 0.1) return init; 
		return tran;
	}

	void RgbdReader::doSlam(int frameID, double dthresh)
	{ 
		
		if(frameID == 0){
			datas->transpointClouds.push_back(datas->pointClouds[frameID]);
			Eigen::Matrix4f tran = Eigen::Matrix4f::Identity();
			datas->trans.push_back(tran);
		}
		else{

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_source(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_target(new pcl::PointCloud<pcl::PointXYZRGB>);
			double variance = 0;
			Eigen::Matrix4f init;
			datas->get_match(frameID,frameID-1,init,variance,false);

			//Eigen::Matrix4f init2 = datas->trans[frameID - 1];
			Eigen::Matrix4f tran = cpuICP(datas->transpointNormals[frameID - 1], datas->transpointClouds[frameID - 1], datas->pointClouds[frameID], datas->trans[frameID - 1].inverse(), init, nHalfXres, nHalfYres, fCoeffX, fCoeffY,dthresh);

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*(datas->pointClouds[frameID]), *transCloud, tran);
			datas->transpointClouds.push_back(transCloud);
			datas->trans.push_back(tran);
		}

		pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
		pcl::IntegralImageNormalEstimation<PointT, pcl::Normal> ne;
		ne.setNormalEstimationMethod (ne.AVERAGE_3D_GRADIENT);
		ne.setMaxDepthChangeFactor(0.03f);
		ne.setNormalSmoothingSize(20.0f);
		ne.setInputCloud(datas->transpointClouds[frameID]);
		ne.compute(*normals);
		datas->transpointNormals.push_back(normals);
		
	}

	std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> RgbdReader::getPointCloud(int frameId)
	{
		std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> ret(datas->transpointClouds[frameId], datas->transpointNormals[frameId]);
		return ret;
	}

	std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> RgbdReader::subPointCloud(Eigen::Matrix4f transform, std::vector<int>& Tags, int blockId, int frameId)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
		pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);

		cloud->height = 480;
		cloud->width = 640;
		cloud->is_dense = false;

		cloud->points.resize(cloud->height * cloud->width);

		normals->height = 480;
		normals->width = 640;
		normals->is_dense = false;

		normals->points.resize(normals->height * normals->width);

		Eigen::Matrix3f rot   = transform.block<3, 3> (0, 0);
        Eigen::Vector3f trans = transform.block<3, 1> (0, 3);

		for(int i = 0; i < cloud->size(); i++)
		{
			if(Tags[i] == blockId)
			{
				cloud->points[i].getVector3fMap () = rot * datas->transpointClouds[frameId]->points[i].getVector3fMap () + trans;
				cloud->points[i].r = datas->transpointClouds[frameId]->points[i].r;
				cloud->points[i].g = datas->transpointClouds[frameId]->points[i].g;
				cloud->points[i].b = datas->transpointClouds[frameId]->points[i].b;
				normals->points[i].getNormalVector3fMap() = rot * datas->transpointNormals[frameId]->points[i].getNormalVector3fMap();
			}
			else
			{
                cloud->points[i].x = 0;
				cloud->points[i].y = 0;
				cloud->points[i].z = 0;
				normals->points[i].normal_x = 0xFFFFFFFF;
				normals->points[i].normal_y = 0xFFFFFFFF;
				normals->points[i].normal_z = 0xFFFFFFFF;

			}
		}
		std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> ret(cloud, normals);
		return ret;
	}
}
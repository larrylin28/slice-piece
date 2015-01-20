#ifndef RGBDREADER_H_
#define RGBDREADER_H_

//pcl;
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>


//opencv;
#include <opencv/cv.h>


namespace xn
{
	class Context;
	class DepthMetaData;  
	class ImageMetaData;  
	class DepthGenerator;        
	class ImageGenerator; 
}


namespace Rgbd
{
	typedef pcl::PointXYZRGB PointT;
    typedef pcl::PointXYZRGBNormal PointN;
    typedef pcl::PointCloud<PointT> PointCloud;
	typedef boost::shared_ptr<xn::Context> ContextPtr;
	typedef boost::shared_ptr<xn::DepthMetaData> DepthMetaPtr;  
	typedef boost::shared_ptr<xn::ImageMetaData> ImageMetaPtr;  
	typedef boost::shared_ptr<xn::DepthGenerator> DepthGeneratorPtr;        
	typedef boost::shared_ptr<xn::ImageGenerator> ImageGeneratorPtr;

	class SdfModel;
	class Segmentation;

	class RgbdDatas
	{
	public:
		friend class RgbdReader;
	
		RgbdDatas();
		~RgbdDatas();

	private:
		//parameters;
		double XtoZ;
		double YtoZ;

		//OpenNI datas;
		ContextPtr context;
		DepthMetaPtr depthMD;  
		ImageMetaPtr imageMD;  
		DepthGeneratorPtr depthGenerator;        
		ImageGeneratorPtr imageGenerator;  

		//point cloud datas;
		std::vector<cv::Mat> descriptors;
		std::vector<std::vector<cv::KeyPoint>> keypoints;
		std::vector<PointCloud::Ptr> pointClouds;
		std::vector<PointCloud::Ptr> transpointClouds;
		std::vector<pcl::PointCloud<pcl::Normal>::Ptr> transpointNormals;
		std::vector<Eigen::Matrix4f> trans;

		//get ransac transformation matrix
		bool get_match(int id_s,int id_t, Eigen::Matrix4f &p_transformation_matrix,double &variance,bool id_s_transfered);
	};
	typedef boost::shared_ptr<RgbdDatas> DatasPtr;

	class RgbdReader
	{
	public:
		RgbdReader();
		~RgbdReader();
		void Read(std::string sFilename, int startFrame, int skip, int frameCount, double dthresh);
		void ShowPointCoud();
		void putDataToSdf(SdfModel* sdf, int frameId);
		void showSegmentation(int frameId, Segmentation* seg);
		int getFrameCount();
		Eigen::Matrix4f getTrans(int frameId);

	private:
		bool initReading(std::string sFilename);
		void readFrame();
		void doSlam(int frameID, double dthresh);

		DatasPtr datas;
		//parameters
		int nHalfXres;	
		int nHalfYres;

		double fCoeffX;
		double fCoeffY;
	};
}
#endif 

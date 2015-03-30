#ifndef _ICP_H
#define _ICP_H

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

Eigen::Matrix4f gpuICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh);
Eigen::Matrix4f cpuICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh);

Eigen::Matrix4f cpuExtraICP(pcl::PointCloud<pcl::Normal>::Ptr old_normals, pcl::PointCloud<pcl::PointXYZRGB>::Ptr old_cloud, pcl::PointCloud<pcl::PointXYZRGB>::Ptr new_cloud, Eigen::Matrix4f mc, Eigen::Matrix4f tr,int xres,int yres,float coeffx,float coeffy, float dthresh, int extraTag = -1, int* oldTags = NULL, int* newTags = NULL, Eigen::Matrix4f extraTr = Eigen::Matrix4f::Identity());

#endif
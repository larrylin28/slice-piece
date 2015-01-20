#ifndef GRAPH_MANAGER_H_
#define GRAPH_MANAGER_H_


#include "g2o/core/sparse_optimizer.h"
#include "g2o/types/slam3d/se3quat.h"
#endif

class GraphManager{
public:
	//相机姿态Graph
	//相机姿态顶点VertexSE3，每个顶点对应一帧，包含参数为:(x, y, z, qx, qy, qz, qw)^T 其中x,y,z为平移矩阵的tx,ty,tz，qx,qy,qz,qw为旋转矩阵的四元数 
	//相机姿态Graph的边EdgeSE3，每条边代表两帧之间的旋转关系；
	//SparseOptimizer优化器
	g2o::SparseOptimizer *optimizer;

	GraphManager();
	~GraphManager();
	//当前解决方案每个顶点初始都是原始坐标，每条边的初始是两帧之间通过RANSAC得到的旋转矩阵
	void addVertex(int id);
	void addEdge(int from,int to,Eigen::Matrix4f* trans,double variance,int inliers);
	void initializeOptimization();

	//数学公式：
	g2o::SE3Quat getEstimate(const Eigen::Matrix4d& m);
	Eigen::Matrix4f* getEuleMatrix(double* est);
	void getJacobian();
};
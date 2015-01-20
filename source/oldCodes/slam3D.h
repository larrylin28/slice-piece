#ifndef GRAPH_MANAGER_H_
#define GRAPH_MANAGER_H_


#include "g2o/core/sparse_optimizer.h"
#include "g2o/types/slam3d/se3quat.h"
#endif

class GraphManager{
public:
	//�����̬Graph
	//�����̬����VertexSE3��ÿ�������Ӧһ֡����������Ϊ:(x, y, z, qx, qy, qz, qw)^T ����x,y,zΪƽ�ƾ����tx,ty,tz��qx,qy,qz,qwΪ��ת�������Ԫ�� 
	//�����̬Graph�ı�EdgeSE3��ÿ���ߴ�����֮֡�����ת��ϵ��
	//SparseOptimizer�Ż���
	g2o::SparseOptimizer *optimizer;

	GraphManager();
	~GraphManager();
	//��ǰ�������ÿ�������ʼ����ԭʼ���꣬ÿ���ߵĳ�ʼ����֮֡��ͨ��RANSAC�õ�����ת����
	void addVertex(int id);
	void addEdge(int from,int to,Eigen::Matrix4f* trans,double variance,int inliers);
	void initializeOptimization();

	//��ѧ��ʽ��
	g2o::SE3Quat getEstimate(const Eigen::Matrix4d& m);
	Eigen::Matrix4f* getEuleMatrix(double* est);
	void getJacobian();
};
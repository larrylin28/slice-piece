#ifndef SEGMENTATION_H
#define SEGMENTATION_H

//pcl
#include "Plane3D.h"
#include "MinCut.h"

namespace Rgbd
{
	class Block
	{
	public:

		//3d space feature
		Eigen::Matrix4f local;
		Plane3D coff;
		double bound[6]; //minx, miny, minz, maxx, maxy, maxz

		//color space feature
		double ELab[3];
		double VLab[3];

		int inliers;
		double energy;

	public:
		Block();
		void clear();
		void addInlier(PointT& p, pcl::Normal& n);
		void updateCore();
		double getDistance(PointT& p);
		double getDistance(Block& b);

	private:
		void getProjectCoord(PointT& p, double coord[3]);
		void getProjectCoord(int x, int y, int z, double coord[3]);

	private:
		//for core calculation
		double EL;
		double Ea;
		double Eb;
		double VL;
		double Va;
		double Vb;
	};

	class BlockGraph
	{
	public:
		BlockGraph(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<Plane3D>& coffs);
		~BlockGraph();
		int graphCut(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<int>& newTag, double disthre);
	private:
		double initialEdge(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, Graph& graph, double disthre);
		void updateDistanceCache(PointT& p, pcl::Normal& n, int fromId, int toId);
		double blockSpaceDistance(int a, int b);
		double blockSpaceDistanceOneDirection(int a, int b);
		double blockColorDistance(int a, int b);
	    double blockEdgeWeight(int a, int b);
	private:
		bool connectedWithKeyEdge(CutTreeNode* root); //judge if is a connected graph
		bool connectedWithKeyEdgeD(CutTreeNode* root); //judge if is a connected graph (Direct Graph)
		double toleration(CutTreeNode* root);
		bool KeyEdgeJudgement(int block_a, int block_b);
		void splitNode(CutTreeNode* root, double maxCutValue);
		std::vector<Block> blocks;
		double** spaceDistanceCache;
		double maxDis;
		double minDis;
		bool** isKeyEdge;
	};

	class Segmentation
	{
	public:
		//parameters
		//int inth; //inlierTreshold
		//int areath; //areaTreshold
	public:
		Segmentation(){}
		~Segmentation(){}

		int nextSegment(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags);
		void calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			                  std::vector<int>& tags, int size, std::vector<int>& newTags, std::vector<Plane3D>& newCoffs);
	};
}

#endif // SEGMENTATION_H

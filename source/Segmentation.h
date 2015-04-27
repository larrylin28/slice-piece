#ifndef SEGMENTATION_H
#define SEGMENTATION_H

//pcl
#include "Plane3D.h"
#include "MinCut.h"

namespace Rgbd
{
	struct Bound
	{
		double bound[6];

		inline double& minX(){ return bound[0];}
		inline double& minY(){ return bound[1];}
		inline double& minZ(){ return bound[2];}
		inline double& maxX(){ return bound[3];}
		inline double& maxY(){ return bound[4];}
		inline double& maxZ(){ return bound[5];}

	    inline Bound()
		{
			bound[0] = 1;
			bound[3] = 0; 
		}

		inline double volume()
		{
			return (maxX() - minX()) * (maxY() - minY()) * (maxZ() - minZ());
		}
	};

	class Block
	{
	public:

		//3d space feature
		Eigen::Matrix4f local;
		//Plane3D coff;
		coffCal coff;
		Bound bound;
		//double bound[6]; //minx, miny, minz, maxx, maxy, maxz

		int belongId; //valid only if it be merged by others;

		//color space feature
		double ELab[3];
		double VLab[3];

		int inliers;
		double energy;

	public:
		Block();
		void clear();
		void addInlier(PointT& p, pcl::Normal& n);
		void combine(Block& b);
		void updateBound(Bound& b, PointT& p, pcl::Normal& n);
		void updateBound(Bound& b, double x, double y, double z);
		void updateCore();
		double getDistance(PointT& p);
		double getDistance(Block& b);
		double intersection(Block& b);
		double CombinerRate(Block& b);
		bool canCombiner(Block& b);
		double insideRate(PointT& p);
		bool isInBlock(PointT& p);
		bool isInBlock(Block& b);
		void output();
		void getProjectCoord(PointT& p, double coord[3]);
		void getProjectCoord(double x, double y, double z, double coord[3]);
		void getProjectCoord(Eigen::Matrix4f loc, double x, double y, double z, double coord[3]);

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
		BlockGraph(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<coffCal>& coffs, std::vector<Block>& v_blocks);
		~BlockGraph();
		int graphCut(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<int>& newTag, double disthre);
		void getBlocks(std::vector<Block>& b);
	private:
		double initialEdge(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, Graph& graph, double disthre);
		void updateDistanceCache(PointT& p, pcl::Normal& n, int fromId, int toId);
		double blockSpaceDistance(int a, int b);
		double blockSpaceDistanceOneDirection(int a, int b);
		double blockColorDistance(int a, int b);
	    double blockEdgeWeight(int a, int b, bool neighbour);
		double blockEdgeWeightOneDirection(int a, int b, bool neighbour);
		double cutThreshold(CutTreeNode* root); //calculate a threshold for a cut node
		bool connectedWithKeyEdge(CutTreeNode* root); //judge if is a connected graph
		bool connectedWithKeyEdgeD(CutTreeNode* root); //judge if is a connected graph (Direct Graph)
		bool KeyEdgeJudgement(int block_a, int block_b);
		void splitNode(CutTreeNode* root, double maxCutValue);
		int setNewTags(CutTreeNode* root, std::vector<int>& newTags); //return new tag count;
	private:
		std::vector<Block> blocks;
		std::vector<Block> virtual_blocks; //block from frames before
		double** spaceDistanceCache;

		double maxDis;
		double minDis;
		std::vector<double> aver_dists;
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

		int nextSegment(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<int>& initTags,
			            std::vector<Block>& blocks, std::vector<Block>& seged_blocks);
		void calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			                  std::vector<int>& tags, int size, std::vector<int>& newTags, 
							  std::vector<int>& foundTags, std::vector<coffCal>& newCoffs);
		void calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			                  std::vector<int>& tags, int size, std::vector<int>& newTags, 
							  std::vector<coffCal>& oldCoffs, std::vector<coffCal>& newCoffs);
		int calSegedTags(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& initTags, 
			             std::vector<int>& tags, std::vector<int>& newTags, std::vector<int>& foundTags,
			             std::vector<Plane3D>& planes, std::vector<Block>& seged_blocks);


	public:
		//std::vector<Block> seged_blocks;
	};
}

#endif // SEGMENTATION_H

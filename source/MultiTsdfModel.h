#ifndef MULTITSDF_MODEL_H
#define MULTITSDF_MODEL_H

#include "Segmentation.h"
#include "SingleTsdfModel.h"

namespace Rgbd
{	
	

	class MultiTsdfModel : public SdfModel
	{
	public:
        MultiTsdfModel();
		~MultiTsdfModel();

		void showSegmentation(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<Block>& blocks, int tagcount);

		int addNewBlock(PointCloud::Ptr cloud, PointNormal::Ptr normals, Block& b, int tag, std::vector<int>& tags);

		void blockCluster(std::vector<int>& blocktags);

		void buildBlocks(PointCloud::Ptr cloud, PointNormal::Ptr normals, Eigen::Matrix4f tran);

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY);


		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz);

		virtual void freeData();

	private:
		Segmentation seg;
		std::vector<Block> seged_blocks;
		std::vector<SingleTsdfModel> nodes;

		//double fractionThre;
		//std::vector<Block> fractions;
	};
}
#endif // MULTITSDF_MODEL_H

#ifndef MULTITSDF_MODEL_H
#define MULTITSDF_MODEL_H

#include "SingleTsdfModel.h"

namespace Rgbd
{	
	

	class MultiTsdfModel : public SdfModel
	{
	public:
        MultiTsdfModel();
		~MultiTsdfModel();

		void showNewRegistration(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, int extraTag, int frameId);

		void showSegmentation(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, std::vector<Block>& blocks, int tagcount);

		int addNewBlock(PointCloud::Ptr cloud, PointNormal::Ptr normals, Block& b, int tag, std::vector<int>& tags);

		void blockCluster(std::vector<int>& blocktags);

		void buildBlocks(PointCloud::Ptr cloud, PointNormal::Ptr normals, Eigen::Matrix4f tran);

		void doNewRegistration(PointCloud::Ptr cloud, PointNormal::Ptr normals, Eigen::Matrix4f tran, RgbdReader* reader);

		virtual void dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader);


		virtual void rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz);

		virtual void freeData();

	private:
		Segmentation seg;
		std::vector<Block> seged_blocks;
		std::vector<SingleTsdfModel> nodes;

		std::vector<std::vector<int>> seged_tags;
		std::vector<std::vector<Eigen::Matrix4f>> extra_Trans;

		//double fractionThre;
		//std::vector<Block> fractions;
	};
}
#endif // MULTITSDF_MODEL_H

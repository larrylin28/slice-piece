#include "RGBDReader.h"
#include "Segmentation.h"
#include "MultiTsdfModel.h"
#include "TsdfModel.h"
#include "MCGenerator.h"

//pcl
#include <pcl/registration/ia_ransac.h>
#include <pcl/visualization/pcl_visualizer.h>


#define GETCOLORR(i) (i%4)*1.0/3
#define GETCOLORG(i) ((i+1)%5)*1.0/4
#define GETCOLORB(i) (((1000-i)+(i/4))%4)*1.0/3

#define ISNAN(x) ((x) != (x)) 

//#define DEBUG_BLOCK_SEGMENT
#define DEBUG_SHOW_BOUND
//#define DEBUG_SHOW_SEG

#define SPECIAL_ID 0

using namespace std;

namespace Rgbd
{	
        MultiTsdfModel::MultiTsdfModel()
		{
			//fractionThre = 0.01;
		}

		MultiTsdfModel::~MultiTsdfModel()
		{
		}

		void MultiTsdfModel::updateTags(RgbdReader* reader)
		{
			for(int i = 0; i < seged_tags.size(); i++)
			{

				for(int j = 0; j < seged_tags[i].size(); j++)
				{
					if(seged_tags[i][j] >= 0 && seged_tags[i][j] < seged_blocks.size())
					{
						seged_tags[i][j] = getbelongID(seged_tags[i][j]);
					}
				}
			}

		
			for(int i = 0; i < seged_tags.size(); i++)
			{
				std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> get = reader->getPointCloud(i);

				pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(get.first);
				pcl::PointCloud<pcl::Normal>::Ptr normals(get.second);
				classifyUntaged(cloud, normals, seged_tags[i]);
			}
		}

		void MultiTsdfModel::classifyUntaged(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags)
		{
			for(int i = 0; i < cloud->size(); i++)
			{
				if(tags[i] < 0 || tags[i] >= seged_blocks.size())
				{
					for(int j = 0;j < seged_blocks.size();j++)
					{
						//belong to seged blocks
						if(getbelongID(j) == j && seged_blocks[j].isInBlock(cloud->at(i)))
						{
							tags[i] = j;
							seged_blocks[j].inliers++;
							break;
						}
					}
				}
			}
		}

		void MultiTsdfModel::showNewRegistration(PointCloud::Ptr cloud, PointNormal::Ptr normals, std::vector<int>& tags, int extraTag, int frameId)
		{
			static boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1;
			if(!vis1)
			{
				static boost::shared_ptr<pcl::visualization::PCLVisualizer> temp(new pcl::visualization::PCLVisualizer ("Reg"));
				vis1 = temp;
				vis1->initCameraParameters ();
				vis1->setBackgroundColor (0.1, 0.1, 0.1);
				vis1->addText("point clouds", 10, 10, "point clouds");
			}

			

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr showCloud(new pcl::PointCloud<pcl::PointXYZRGB>);

			Eigen::Matrix4f tran = extra_Trans[frameId][extraTag];
			cout<<tran;

			for(int i = 0; i < cloud->size(); i++)
			{
				if(tags[i] == extraTag)
				{
					showCloud->push_back(cloud->at(i));
				}
			}

			pcl::PointCloud<pcl::PointXYZRGB>::Ptr transCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
			pcl::transformPointCloud(*showCloud, *transCloud, tran);

			
			pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb_cloud(transCloud);
			stringstream ss;
			ss<<"cloud"<< frameId;
			vis1->addPointCloud<pcl::PointXYZRGB> (transCloud, rgb_cloud, ss.str());
			vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
			vis1->resetStoppedFlag();
			
			while (!vis1->wasStopped())
			{
				//在此处可以添加其他处理
				vis1->spinOnce (100);
				//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
			}
		}

		void MultiTsdfModel::showSegmentation(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, std::vector<Block>& blocks, int tagcount)
		{
			boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1 (new pcl::visualization::PCLVisualizer ("Seg"));

			vis1->initCameraParameters ();
			vis1->setBackgroundColor (0.1, 0.1, 0.1);
			vis1->addText("point clouds", 10, 10, "point clouds");

			vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> showClouds;

			for(int i = 0; i < tagcount; i++)
			{
				showClouds.push_back(pcl::PointCloud<pcl::PointXYZRGB>::Ptr(new pcl::PointCloud<pcl::PointXYZRGB>));
			}

			for(int i = 0; i < cloud->size(); i++)
			{
				if(tags[i] >= 0 && tags[i] < tagcount)
				{
					showClouds[tags[i]]->push_back(cloud->at(i));
				}
			}

			for(int i = 0; i < tagcount; i++)
			{
				//pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> single_color(showClouds[i]);
				double r = GETCOLORR(i);
				double g = GETCOLORG(i);
				double b = GETCOLORB(i);
				pcl::visualization::PointCloudColorHandlerCustom<PointT> single_color(showClouds[i], r*255, g*255, b*255);
				stringstream ss;
				ss<<"cloud"<<i;
				vis1->addPointCloud<pcl::PointXYZRGB> (showClouds[i], single_color, ss.str());
				vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
			}

#ifdef DEBUG_SHOW_BOUND
			static int vertex[8][3] = {{0,1,2}, {0,1,5}, {0,4,2}, {0,4,5},{3,1,2}, {3,1,5}, {3,4,2}, {3,4,5}};
			static int facet[6][4] = {{0,1,3,2}, {0,4,5,1}, {4,6,7,5}, {2,3,7,6}, {1,5,7,3}, {0,2,6,4}};
			for(int i = 0; i < blocks.size(); i++)
			{
				if(blocks[i].belongId >= 0 && blocks[i].belongId != i) continue;
				double r = GETCOLORR(i);
			    double g = GETCOLORG(i);
			    double b = GETCOLORB(i);
				Block& block = blocks[i];
				double vertex_p[8][3];
				Eigen::Matrix4f inv = block.local.inverse();
				for(int j = 0; j < 8 ;j++)
				{
					double coord[3];
					block.getProjectCoord(inv,block.bound.bound[vertex[j][0]], block.bound.bound[vertex[j][1]], block.bound.bound[vertex[j][2]],vertex_p[j]);
				}
				for(int j = 0; j < 6; j++)
				{
					pcl::PointCloud<pcl::PointXYZRGB>::Ptr planeBounds (new pcl::PointCloud<pcl::PointXYZRGB>);
					planeBounds->width    = 4;
					planeBounds->height   = 1;
					planeBounds->is_dense = false;
					planeBounds->points.resize (planeBounds->width * planeBounds->height);
					for(int l = 0;l < 4;l++){
						pcl::PointXYZRGB& pt_s = planeBounds->points[l];
						pt_s.x = vertex_p[facet[j][l]][0];
						pt_s.y = vertex_p[facet[j][l]][1];
						pt_s.z = vertex_p[facet[j][l]][2];
						pt_s.r = 1;
						pt_s.g = 0;
						pt_s.b = 0;
			     
					}
					stringstream ss;
				    ss<<"block"<<i<<","<<j;
					vis1->addPolygon<pcl::PointXYZRGB>(planeBounds,r, g, b,ss.str());
				}
			}
#endif
			
			

			while (!vis1->wasStopped())
			{
				//在此处可以添加其他处理
				vis1->spinOnce (100);
				//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
			}
		}

		int MultiTsdfModel::addNewBlock(PointCloud::Ptr cloud, PointNormal::Ptr normals, Block& b, int tag, std::vector<int>& tags)
		{
			int newtag = -1;
			if(b.coff.seged_id >= 0 && (seged_blocks[b.coff.seged_id].isInBlock(b) || seged_blocks[b.coff.seged_id].canCombiner(b)))
			{
				newtag = b.coff.seged_id;
			}else{
				double minrate = 1.3;
				for(int i = 0;i < seged_blocks.size();i++)
				{
					//belong to seged blocks
					if(seged_blocks[i].isInBlock(b))
					{
						newtag = i;
						break;
					}
					
					double rate = seged_blocks[i].CombinerRate(b);
					if(rate < minrate)
					{
						minrate = rate;
						newtag = i;
						break;
					}
					

				}
			}

			if(newtag >= 0)
			{
				for(int i = 0; i < cloud->size(); i++)
				{
					if(tags[i] == tag)
					{
						seged_blocks[newtag].addInlier(cloud->at(i), normals->at(i));
					}
				}
			}
			else
			{
				newtag = seged_blocks.size();
				seged_blocks.push_back(b);
				
			}
			return newtag;
		}

		int MultiTsdfModel::getbelongID(int blockID)
		{
			int id = blockID;
			while(seged_blocks[id].belongId >= 0 && seged_blocks[id].belongId != id) id = seged_blocks[id].belongId;
			seged_blocks[blockID].belongId = id;
			return id;
		}

		void MultiTsdfModel::blockCluster()
		{
			for(int i = 0;i < seged_blocks.size();i++)
			{
				if(getbelongID(i) == i)
				{
					for(int j = i+1; j < seged_blocks.size(); j++)
					{
						if(getbelongID(j) == j)
						{
							if(seged_blocks[i].isInBlock(seged_blocks[j])) // || seged_blocks[i].canCombiner(seged_blocks[j]))
							{
								seged_blocks[i].combine(seged_blocks[j]);
								seged_blocks[j].belongId = i;
								cout<<"Merged!";
							}
						}
					}
				}
			}

		}

		void MultiTsdfModel::buildBlocks(PointCloud::Ptr cloud, PointNormal::Ptr normals, Eigen::Matrix4f tran, RgbdReader* reader)
		{
			int frameId = seged_tags.size();
			vector<int> tags(cloud->size(), -1);
			vector<int> inittags(cloud->size(), -1);
			if(frameId > 0) reader->initaltags(frameId, seged_tags[frameId -1], inittags, 0.1);

		    vector<Block> blocks;
		    int count = seg.nextSegment(cloud, normals, tags, inittags, blocks, seged_blocks);

			//if(frameId > 0) showSegmentation(cloud, normals, inittags, blocks, seged_blocks.size());

#ifdef DEBUG_SHOW_SEG
			//show segmentation
			showSegmentation(cloud, normals, tags, blocks, count);
#endif

			
			//add blocks
			vector<int> newblocktags(blocks.size(), -1);
			for(int i = 0; i < blocks.size(); i++)
			{
				int newtag = addNewBlock(cloud, normals, blocks[i], i, tags);
				newblocktags[i] = newtag;
				//newblocktags[i] = seged_blocks.size();
				//seged_blocks.push_back(blocks[i]);
			}

			//cluster blocks
			blockCluster();

			for(int i = 0; i < tags.size(); i++)
			{
				if(tags[i] >= 0 && tags[i] < count)
				{
					tags[i] = getbelongID(newblocktags[tags[i]]);
				}
			}
			//classifyUntaged(cloud, normals, tags);

#ifdef DEBUG_BLOCK_SEGMENT
				for(int i = 0; i < seged_blocks.size(); i++)
				{
					cout<<"Block "<<i<<": ";
					seged_blocks[i].output();
				}
				cout<<endl;
#endif

#ifdef DEBUG_SHOW_SEG
			//show segmentation
			showSegmentation(cloud, normals, tags, seged_blocks, seged_blocks.size());
#endif

			seged_tags.push_back(tags);

			
			
			
		}

		void MultiTsdfModel::doNewRegistration(PointCloud::Ptr cloud, PointNormal::Ptr normals, Eigen::Matrix4f tran, RgbdReader* reader)
		{
			int frameId = seged_tags.size() - 1;
			int special_id = SPECIAL_ID;
			std::vector<Eigen::Matrix4f> trs;
			extra_Trans.push_back(trs);
			for(int i = 0; i < seged_blocks.size(); i++)
			{
				if(i == special_id && frameId > 0 && extra_Trans[frameId - 1].size() > special_id)
				{
					Eigen::Matrix4f tr = reader->getExtraSlam(frameId, 0.1, i, seged_tags[frameId - 1], seged_tags[frameId], extra_Trans[frameId-1][i]);
					Eigen::Matrix4f gtr = extra_Trans[frameId - 1][i]*tr;
					extra_Trans[frameId].push_back(gtr);
				}
				else
				{
					extra_Trans[frameId].push_back(Eigen::Matrix4f::Identity());	
				}
			}
			//show new registration
			if(special_id < seged_blocks.size())
			{
				showNewRegistration(cloud, normals, seged_tags[frameId], special_id, frameId);
			}
		}

		void  MultiTsdfModel::MeshGenForBlock(int blockID, RgbdReader* reader, GLMesh* mesh, double devide, double sample_div)
		{
			double r = GETCOLORR(blockID);
			double g = GETCOLORG(blockID);
			double b = GETCOLORB(blockID);
				
			static boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1;
			if(!vis1)
			{
				static boost::shared_ptr<pcl::visualization::PCLVisualizer> temp(new pcl::visualization::PCLVisualizer ("Reg"));
				vis1 = temp;
				vis1->initCameraParameters ();
				vis1->setBackgroundColor (0.1, 0.1, 0.1);
				vis1->addText("point clouds", 10, 10, "point clouds");
			}
			
			if(getbelongID(blockID) != blockID) return;
			

			Block& block = seged_blocks[blockID];
			double extra_bound = devide * 5;
			//TsdfModel tsdf(block.local.inverse(), block.bound.maxX() + extra_bound, block.bound.maxY() + extra_bound, block.bound.maxZ() + extra_bound,
			//	                        block.bound.minX() - extra_bound, block.bound.minY() - extra_bound, block.bound.minZ() - extra_bound,
			//							devide, devide * 10, devide / 10, false);
			//MCGenerator gen(&tsdf);

			int frameCount = reader->getFrameCount();
			Eigen::Matrix4f gtr = Eigen::Matrix4f::Identity();
			std::vector<Eigen::Matrix4f> trs;
			trs.push_back(gtr);
			for(int i = 0; i < frameCount; i++)
			{
				if(i > 0)
				{

					Eigen::Matrix4f tr = reader->getExtraSlam(i, 0.1, blockID, seged_tags[i - 1], seged_tags[i], trs[i-1]);
					Eigen::Matrix4f gtr = trs[i-1]*tr;
					trs.push_back(gtr);
				}

				Eigen::Matrix4f local_tran = trs[i];
				Eigen::Matrix4f global_tran = local_tran * reader->getTrans(i);

				std::pair<pcl::PointCloud<pcl::PointXYZRGB>::Ptr, pcl::PointCloud<pcl::Normal>::Ptr> sub = reader->subPointCloud(local_tran, seged_tags[i], blockID, i);

				pcl::PointCloud<pcl::PointXYZRGB>::Ptr subCloud(sub.first);
				pcl::PointCloud<pcl::Normal>::Ptr subNormals(sub.second);

				pcl::visualization::PointCloudColorHandlerCustom<PointT> single_color(subCloud, r*255, g*255, b*255);
				stringstream ss;
				ss<<"cloud"<<blockID<<","<< i;
				vis1->addPointCloud<pcl::PointXYZRGB> (subCloud, single_color, ss.str());
				vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
				vis1->resetStoppedFlag();

				//tsdf.dataFusion(subCloud, subNormals, global_tran, reader->nHalfXres, reader->nHalfYres, reader->fCoeffX, reader->fCoeffY, reader);  
			}

			
			
			while (!vis1->wasStopped())
			{
				//在此处可以添加其他处理
				vis1->spinOnce (100);
				//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
			}

			/*
			double sizex = block.bound.maxX() - block.bound.minX() + 2 * extra_bound;
			double sizey = block.bound.maxY() - block.bound.minY() + 2 * extra_bound;
			double sizez = block.bound.maxZ() - block.bound.minZ() + 2 * extra_bound;
			int xlen = (sizex / sample_div) + 0.5;
			int ylen = (sizey / sample_div) + 0.5;
			int zlen = (sizez / sample_div) + 0.5;
			gen.GenerateMeshes(block.bound.minX() - extra_bound, block.bound.minY() - extra_bound, block.bound.minZ() - extra_bound,
				               xlen, ylen, zlen, sample_div, mesh);

			tsdf.freeData();
			*/
		}

		void  MultiTsdfModel::showFinalSegment(RgbdReader* reader)
		{

			boost::shared_ptr<pcl::visualization::PCLVisualizer> vis1 (new pcl::visualization::PCLVisualizer ("Reg"));

			vis1->initCameraParameters ();
			vis1->setBackgroundColor (0.1, 0.1, 0.1);
			vis1->addText("point clouds", 10, 10, "point clouds");

			Bound all;
			Block ball;
			int frameCount = reader->getFrameCount();
			for(int i = 0; i < frameCount; i++)
			{
				pcl::PointCloud<pcl::PointXYZRGB>::Ptr showCloud(new pcl::PointCloud<pcl::PointXYZRGB>);
				showCloud->height = 480;
				showCloud->width = 640;
				showCloud->is_dense = false;

				showCloud->points.resize(showCloud->height * showCloud->width);
				pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud = reader->getPointCloud(i).first;
				for(int j = 0; j < cloud->size(); j++)
				{
					showCloud->points[j].x = cloud->points[j].x;
					showCloud->points[j].y = cloud->points[j].y;
					showCloud->points[j].z = cloud->points[j].z;
					if(!ISNAN(cloud->points[j].x))
					{
						all.maxX() = max((double)cloud->points[j].x, all.maxX());
						all.maxY() = max((double)cloud->points[j].y, all.maxY());
						all.maxZ() = max((double)cloud->points[j].z, all.maxZ());
						all.minX() = min((double)cloud->points[j].x, all.minX());
						all.minY() = min((double)cloud->points[j].y, all.minY());
						all.minZ() = min((double)cloud->points[j].z, all.minZ());
					}
					int tag = seged_tags[i][j];
					if(tag >= 0)
					{
						showCloud->points[j].r = GETCOLORR(tag)*255;
						showCloud->points[j].g = GETCOLORG(tag)*255;
						showCloud->points[j].b = GETCOLORB(tag)*255;
					}
				}

				pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb_cloud(showCloud);
				stringstream ss;
				ss<<"cloud"<< i;
				vis1->addPointCloud<pcl::PointXYZRGB> (showCloud, rgb_cloud, ss.str());
				vis1->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, ss.str());
			}

			double total = 0;
			int sum = 0;
			for(int i = 0; i < seged_blocks.size();i++)
			{
				Bound& bound = seged_blocks[i].bound;
				double v = (bound.maxX() - bound.minX()) * (bound.maxY() - bound.minY()) * (bound.maxZ() - bound.minZ());
				total += v;

				double dense = seged_blocks[i].inliers / v;
				sum += seged_blocks[i].inliers;
				cout<<dense<<endl;
			}

			double v = (all.maxX() - all.minX()) * (all.maxY() - all.minY()) * (all.maxZ() - all.minZ());
			cout<<"All: "<<sum / v<<endl;

			cout<<"Total"<<total;
			
			cout<<"All"<<v;

			while (!vis1->wasStopped())
			{
				//在此处可以添加其他处理
				vis1->spinOnce (100);
				//boost::this_thread::sleep (boost::posix_time::microseconds (100000));
			}

		}

		void MultiTsdfModel::dataFusion(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			            Eigen::Matrix4f tran, int nHalfXres, int nHalfYres, double fCoeffX, double fCoeffY, RgbdReader* reader)

		{
			buildBlocks(cloud, normals, tran, reader);

			/*
			for(int i = nodes.size(); i < seged_blocks.size(); i++)
			{
				SingleTsdfModel node(seged_blocks[i]);
				nodes.push_back(node);
			}
			*/
			//doNewRegistration(cloud, normals, tran, reader);

			/*
			int id = SPECIAL_ID;
			int frameId = seged_tags.size() - 1;
			nodes[id].dataFusion(cloud, normals, tran, nHalfXres, nHalfYres, fCoeffX, fCoeffY, reader, id, seged_tags[frameId]);
			*/
		}

		void MultiTsdfModel::rayCast(double ix, double iy, double ax, double ay, double devide,
			                 Eigen::Matrix4f tran, 
			                 float* p_dists, float* p_nx, float* p_ny, float* p_nz)
		{
		}

		void MultiTsdfModel::getSample(float sdevide, float sminX, float sminY, float sminZ, int sXlen, int sYlen, int sZlen, float* absrate)
		{
		}

		void MultiTsdfModel::freeData()
		{
		}
}
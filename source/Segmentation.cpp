#include "Segmentation.h"
#include "Util.h"
#include <map>
#include <queue>

using namespace std;

#define e 2.718281828459
#define DIV 0.003125



namespace Rgbd
{
	
	Block::Block():EL(0),Ea(0),Eb(0),VL(0),Va(0),Vb(0),inliers(0),belongId(-1)
	{
	}

	void Block::clear()
	{
		EL = 0;
		Ea = 0;
		Eb = 0;
		VL = 0;
		Va = 0;
		Vb = 0;
		inliers = 0;
	}

	void Block::getProjectCoord(PointT& p, double coord[3])
	{
		 coord[0] = local(0,0)*p.x+local(0,1)*p.y+local(0,2)*p.z + local(0,3)*1;
	     coord[1] = local(1,0)*p.x+local(1,1)*p.y+local(1,2)*p.z + local(1,3)*1;
	     coord[2] = local(2,0)*p.x+local(2,1)*p.y+local(2,2)*p.z + local(2,3)*1;
	}

	void Block::getProjectCoord(double x, double y, double z, double coord[3])
	{
		 coord[0] = local(0,0)*x+local(0,1)*y+local(0,2)*z + local(0,3)*1;
	     coord[1] = local(1,0)*x+local(1,1)*y+local(1,2)*z + local(1,3)*1;
	     coord[2] = local(2,0)*x+local(2,1)*y+local(2,2)*z + local(2,3)*1;
	}


	void Block::getProjectCoord(Eigen::Matrix4f loc, double x, double y, double z, double coord[3])
	{
		 coord[0] = loc(0,0)*x+loc(0,1)*y+loc(0,2)*z + loc(0,3)*1;
	     coord[1] = loc(1,0)*x+loc(1,1)*y+loc(1,2)*z + loc(1,3)*1;
	     coord[2] = loc(2,0)*x+loc(2,1)*y+loc(2,2)*z + loc(2,3)*1;
	}

	inline double out_range(double l, double r, double x)
	{
		if(x < l) return (l - x);
		else if(x > r) return (x - r);
		else return 0;
	}

	double Block::intersection(Block& b)
	{
		static int vertex[8][3] = {{0,1,2}, {0,1,5}, {0,4,2}, {0,4,5},{3,1,2}, {3,1,5}, {3,4,2}, {3,4,5}};
		Bound newb;
		Eigen::Matrix4f inv = b.local.inverse();
		for(int i = 0; i < 8; i++)
		{
			double coord[3];
			getProjectCoord(inv,b.bound.bound[vertex[i][0]], b.bound.bound[vertex[i][1]], b.bound.bound[vertex[i][2]],coord);
			updateBound(newb, coord[0], coord[1], coord[2]);
		}
		double insX = (min(newb.maxX(), bound.maxX()) - max(newb.minX(), bound.minX()));
		double insY = (min(newb.maxY(), bound.maxY()) - max(newb.minY(), bound.minY()));
		double insZ = (min(newb.maxZ(), bound.maxZ()) - max(newb.minZ(), bound.minZ()));
		double newbVolume = newb.volume();
		if(insX > 0 && insY > 0 && insZ > 0) return insX * insY * insZ / newbVolume;
		else return 0;
	}

	bool Block::canCombiner(Block& b)
	{
		double angthre = 0.85;
		double angle = p_angle(coff, b.coff);
		double rate = CombinerRate(b);
		//if(rate <= 1.0) return true;
		if(angle > angthre && rate < 1.3) return true;
		else return false;
	}

	double Block::CombinerRate(Block& b)
	{
		static int vertex[8][3] = {{0,1,2}, {0,1,5}, {0,4,2}, {0,4,5},{3,1,2}, {3,1,5}, {3,4,2}, {3,4,5}};
		Bound newb;
		Eigen::Matrix4f inv = b.local.inverse();
		for(int i = 0; i < 8; i++)
		{
			double coord[3];
			getProjectCoord(inv,b.bound.bound[vertex[i][0]], b.bound.bound[vertex[i][1]], b.bound.bound[vertex[i][2]],coord);
			updateBound(newb, coord[0], coord[1], coord[2]);
		}
		double combX = (max(newb.maxX(), bound.maxX()) - min(newb.minX(), bound.minX()));
		double combY = (max(newb.maxY(), bound.maxY()) - min(newb.minY(), bound.minY()));
		double combZ = (max(newb.maxZ(), bound.maxZ()) - min(newb.minZ(), bound.minZ()));
		double newbVolume = newb.volume();
		double selfVolume = bound.volume();
		double rate = combX * combY * combZ / (newbVolume + selfVolume);
		return rate;
	}

	double Block::insideRate(PointT& p)
	{
		double inside = isInBlock(p) ? 1.0 : 0.0;
		return 1;
	}

	bool Block::isInBlock(Block& b)
	{
		static int vertex[8][3] = {{0,1,2}, {0,1,5}, {0,4,2}, {0,4,5},{3,1,2}, {3,1,5}, {3,4,2}, {3,4,5}};
		Bound newb;
		Eigen::Matrix4f inv = b.local.inverse();
		for(int i = 0; i < 8; i++)
		{
			double coord[3];
			getProjectCoord(inv,b.bound.bound[vertex[i][0]], b.bound.bound[vertex[i][1]], b.bound.bound[vertex[i][2]],coord);
			updateBound(newb, coord[0], coord[1], coord[2]);
		}
		if(newb.maxX() <= bound.maxX() && newb.minX() >= bound.minX() && newb.maxY() <= bound.maxY() && newb.minY() >= bound.minY() && newb.maxZ() <= bound.maxZ() && newb.minZ() >= bound.minZ())
	       return true;
		else 
		   return false;
	}

	bool Block::isInBlock(PointT& p)
	{
		double coord[3];
		getProjectCoord(p, coord);
		if(coord[0] > bound.minX() && coord[0] < bound.maxX() && coord[1] > bound.minY() && coord[1] < bound.maxY() && coord[2] > bound.minZ() && coord[2] < bound.maxZ()) return true;
		return false;
	}

	double Block::getDistance(PointT& p)
	{
		double coord[3];
		getProjectCoord(p, coord);
		double ox = out_range(bound.minX(), bound.maxX(), coord[0]);
		double oy = out_range(bound.minY(), bound.maxY(), coord[1]);
		double oz = out_range(bound.minZ(), bound.maxZ(), coord[2]);
		
		double v1 = (bound.maxX() - bound.minX()) * (bound.maxY() - bound.minY()) * (bound.maxZ() - bound.minZ());
		double v2 = (bound.maxX() - bound.minX() + ox) * (bound.maxY() - bound.minY() + oy) * (bound.maxZ() - bound.minZ() + oz);
		//return sqrt(v2 - v1);
		return pow(v2-v1,1.0/3.0);
		//return (v2 - v1) / v1;
	}

	double Block::getDistance(Block& b)
	{
		
		return 0;
	}

	void Block::updateBound(Bound& b, double x, double y, double z)
	{
		double coord[3];
		getProjectCoord(x, y, z, coord);
		if(b.minX() <= b.maxX())
		{
			b.minX() = min(b.minX(), coord[0]);
			b.minY() = min(b.minY(), coord[1]);
			b.minZ() = min(b.minZ(), coord[2]);
			b.maxX() = max(b.maxX(), coord[0]);
			b.maxY() = max(b.maxY(), coord[1]);
			b.maxZ() = max(b.maxZ(), coord[2]);
		}
		else
		{
			b.minX() = coord[0];
			b.minY() = coord[1];
			b.minZ() = coord[2];
			b.maxX() = coord[0];
			b.maxY() = coord[1];
			b.maxZ() = coord[2];
		}
	}

	void Block::updateBound(Bound& b, PointT& p, pcl::Normal& n)
	{
		double coord[3];
		getProjectCoord(p, coord);
		if(b.minX() <= b.maxX())
		{
			b.minX() = min(b.minX(), coord[0]);
			b.minY() = min(b.minY(), coord[1]);
			b.minZ() = min(b.minZ(), coord[2]);
			b.maxX() = max(b.maxX(), coord[0]);
			b.maxY() = max(b.maxY(), coord[1]);
			b.maxZ() = max(b.maxZ(), coord[2]);
		}
		else
		{
			b.minX() = coord[0];
			b.minY() = coord[1];
			b.minZ() = coord[2];
			b.maxX() = coord[0];
			b.maxY() = coord[1];
			b.maxZ() = coord[2];
		}
	}
                     
	void Block::addInlier(PointT& p, pcl::Normal& n)
	{	
		//3d space
		//project coord
		updateBound(bound, p, n);

		//color space
		double L, a, b;
		RGB2Lab(p.r, p.g, p.b, L, a, b);
		meanVar(EL, VL, L, inliers);
		meanVar(Ea, Va, a, inliers);
		meanVar(Eb, Vb, b, inliers);

		inliers++;
	}

	void Block::combine(Block& b)
	{
		static int vertex[8][3] = {{0,1,2}, {0,1,5}, {0,4,2}, {0,4,5},{3,1,2}, {3,1,5}, {3,4,2}, {3,4,5}};
		Bound newb;
		Eigen::Matrix4f inv = b.local.inverse();
		for(int i = 0; i < 8; i++)
		{
			double coord[3];
			getProjectCoord(inv,b.bound.bound[vertex[i][0]], b.bound.bound[vertex[i][1]], b.bound.bound[vertex[i][2]],coord);
			updateBound(newb, coord[0], coord[1], coord[2]);
		}
		bound.maxX() = max(newb.maxX(), bound.maxX());
		bound.maxY() = max(newb.maxY(), bound.maxY());
		bound.maxZ() = max(newb.maxZ(), bound.maxZ());
		bound.minX() = min(newb.minX(), bound.minX());
		bound.minY() = min(newb.minY(), bound.minY());
		bound.minZ() = min(newb.minZ(), bound.minZ());
		inliers += b.inliers;
	}
	
	void Block::updateCore()
	{
		ELab[0] = EL;
		ELab[1] = Ea;
		ELab[2] = Eb;
		VLab[0] = VL / inliers;
		VLab[1] = Va / inliers;
		VLab[2] = Vb / inliers;
	}

	void Block::output()
	{
		//cout<<"coff:("<<coff.a<<","<<coff.b<<","<<coff.c<<","<<coff.d<<") ";
		cout<<"bound:(("<<bound.minX()<<","<<bound.maxX()<<"),("<<bound.minY()<<","<<bound.maxY()<<"),("<<bound.minZ()<<","<<bound.maxZ()<<")) ";
		cout<<"inliers:"<<inliers<<endl;
	}

	BlockGraph::BlockGraph(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, std::vector<coffCal>& coffs, std::vector<Block>& v_blocks)
		:blocks(coffs.size()),isKeyEdge(NULL),aver_dists(coffs.size(), 0)
		
	{
		for(int i = 0; i < coffs.size(); i++)
		{
			blocks[i].coff = coffs[i];

			bool axisType = true;
			if(axisType)
			{
				//cal rotate matrix
				Eigen::Matrix3f before;
				Eigen::Matrix3f after;

				after << 1, 0 ,0,
						 0, 1, 0,
						 0, 0, 1;

				for(int col = 0; col < 3; col++)
				{
					for(int row = 0; row < 3; row++)
					{
						before(row, col) = blocks[i].coff.vec[col][row];
					}
				}

				Eigen::Matrix3f rotate = after * before.inverse();
			
				blocks[i].local << 1, 0, 0, 0,
								   0, 1, 0, 0,
								   0, 0, 1, 0,
								   0, 0 ,0, 1;

				for(int ri = 0; ri < 3; ri++)
					for(int rj = 0; rj < 3; rj++)
						blocks[i].local(ri, rj) = rotate(ri, rj);
			}
			else
			{
				//cal rotate matrix
				double before[3];
				double after[3];
				before[0] = coffs[i].getCoff().a;
				before[1] = coffs[i].getCoff().b;
				before[2] = coffs[i].getCoff().c;

				after[0] = 0;
				after[1] = 0;
				after[2] = 1;

				double rotate[4][4];
				Calculation4d(before, after, rotate);

				for(int ri = 0; ri < 4; ri++)
					for(int rj = 0; rj < 4; rj++)
						blocks[i].local(ri, rj) = rotate[ri][rj];
			}
			
		}

		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < blocks.size())
			{
				blocks[tags[i]].addInlier(cloud->at(i), normals->at(i));
			}
		}
		for(int i = 0; i < coffs.size(); i++)
		{
			blocks[i].updateCore();
		}

		//for(int i = 0; i < v_blocks.size(); i++)
		//{
		//	virtual_blocks.push_back(v_blocks[i]);
		//}

		//initial space distance cache
		spaceDistanceCache = new double*[blocks.size()];
		for(int i = 0; i < blocks.size(); i++)
		{
			spaceDistanceCache[i] = new double[blocks.size()];
			for(int j = 0; j < blocks.size(); j++)
			{
				spaceDistanceCache[i][j] = 0;

			}
		}

		
		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < blocks.size())
			{
				for(int j = 0; j < blocks.size(); j++)
				{
					if(tags[i] != j)
					{
						updateDistanceCache(cloud->at(i), normals->at(i), tags[i], j);
					}
				}
			}
		}

		maxDis = 0;
		minDis = 10000000.0;
		for(int i = 0; i < blocks.size(); i++)
		{
			for(int j = 0; j < blocks.size(); j++)
			{
				spaceDistanceCache[i][j] /= blocks[i].inliers;
				aver_dists[i] += spaceDistanceCache[i][j];
				maxDis = max(maxDis, spaceDistanceCache[i][j]);
				minDis = min(minDis, spaceDistanceCache[i][j]);
			}
			aver_dists[i] /= blocks.size();
		}

		
	}

	BlockGraph::~BlockGraph()
	{
		for(int i = 0; i < blocks.size(); i++)
		{
			delete[] spaceDistanceCache[i];
		}
		delete[] spaceDistanceCache;
	}

	void BlockGraph::getBlocks(std::vector<Block>& b)
	{
		for(int i = 0;i < blocks.size(); i++)
		{
			b.push_back(blocks[i]);
		}
	}

	bool BlockGraph::connectedWithKeyEdgeD(CutTreeNode* root)
	{
		assert(isKeyEdge != NULL);
		//bfs search if connected
		queue<int> search;
		vector<bool> reached(root->data->vertex_count, false);
		int reach_count = 0;
		search.push(0);
		reached[0] = true;
		while(true)
		{
			int preNode = -1;		
			while(!search.empty())
			{
				int t = search.front();
				search.pop();
				reach_count++;
				for(int i = 0;i < root->data->vertex_count; i++)
				{
					if(!reached[i] &&  root->data->haveEdge(t,i))
					{
						if(isKeyEdge[root->data->getVertex(t)][root->data->getVertex(i)])
						{
							reached[i] = true;
							search.push(i);
						}
						else if(isKeyEdge[root->data->getVertex(i)][root->data->getVertex(t)])
						{
							preNode = i;
						}
					}
				}
			}
			if(reach_count == root->data->vertex_count) return true;
			if(preNode < 0) return false;
			search.push(preNode);
		    reached[preNode] = true;
		}
		return reach_count == root->data->vertex_count;
	}

	bool BlockGraph::connectedWithKeyEdge(CutTreeNode* root)
	{
		assert(isKeyEdge != NULL);
		//bfs search if connected
		queue<int> search;
		vector<bool> reached(root->data->vertex_count, false);
		int reach_count = 0;
		search.push(0);
		reached[0] = true;
		while(!search.empty())
		{
			int t = search.front();
			search.pop();
			reach_count++;
			for(int i = 0;i < root->data->vertex_count; i++)
			{
				if(!reached[i] && isKeyEdge[root->data->getVertex(t)][root->data->getVertex(i)] && root->data->haveEdge(t,i))
				{
					reached[i] = true;
					search.push(i);
				}
			}
		}
		return reach_count == root->data->vertex_count;
	}

	double BlockGraph::cutThreshold(CutTreeNode* root)
	{
		vector<double> dists(root->data->vertex_count, 0);
		vector<double> weights(root->data->vertex_count, 0);
		for(int i = 0; i < root->data->vertex_count; i++)
		{
			for(int j = 0; j < root->data->vertex_count; j++)
			{
				dists[i] += blockSpaceDistanceOneDirection(root->data->getVertex(i),root->data->getVertex(j));
			}
			dists[i] /= root->data->vertex_count;
			weights[i] = 1 / (dists[i] + 0.0001);
		}

		return 18;
	}

	void BlockGraph::splitNode(CutTreeNode* root, double maxCutValue)
	{
		if(root->data->vertex_count <= 1) return;

		vector<int> setS;
		vector<int> setT;
		double cutvalue = root->data->minCut(setS, setT);
		//double tolerateCut = toleration(root);
		bool isConnected = connectedWithKeyEdge(root);
		double thre = cutThreshold(root);

		//if(cutvalue > maxCutValue)
		//{
		     if(isConnected) return;
		//}

		//if(cutvalue > thre) return;
		//double maxC = max(cutvalue, maxCutValue);

		//ouput cut
		
#ifdef DEBUG_PLANE_SEGMENT
			cout<<"CUT:("<<root->data->getVertex(0);
			for(int i = 1; i < root->data->vertex_count; i++)
			{
				cout<<","<<root->data->getVertex(i);
			}
			cout<<")TO=> ("<<root->data->getVertex(setS[0]);
			for(int i = 1; i < setS.size(); i++)
			{
				cout<<","<<root->data->getVertex(setS[i]);
			}
			cout<<")AND("<<root->data->getVertex(setT[0]);
			for(int i = 1; i < setT.size(); i++)
			{
				cout<<","<<root->data->getVertex(setT[i]);
			}
			cout<<") Cut Value:"<<cutvalue<<endl;
#endif
		

		root->split(setS, setT);
		if(root->left) splitNode(root->left, maxCutValue);
		if(root->right) splitNode(root->right,maxCutValue);
	}

	void BlockGraph::updateDistanceCache(PointT& p, pcl::Normal& n, int fromId, int toId)
	{
		Block& to = blocks[toId]; 
		spaceDistanceCache[fromId][toId] += to.getDistance(p);
	}

	int BlockGraph::graphCut(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, vector<int>& newTags, double disthre)
	{
		Graph* root_g = new Graph(blocks.size());
		
		isKeyEdge = new bool*[root_g->vertex_count];
		for(int i = 0; i < root_g->vertex_count; i++)
		{
			isKeyEdge[i] = new bool[root_g->vertex_count];
			for(int j = 0; j < root_g->vertex_count; j++)
			{
                 isKeyEdge[i][j] = false;
			}
		}
		double m_edge = initialEdge(cloud,normals, tags, *root_g, disthre);

#ifdef DEBUG_PLANE_SEGMENT
			for(int i = 0; i < root_g->vertex_count; i++)
			{
				for(int j = i + 1; j < root_g->vertex_count; j++)
				{
					if(root_g->haveEdge(i,j))
					{
						if(isKeyEdge[i][j])
						{
							cout<<"KeyE("<<i<<","<<j<<"): Space:"<<root_g->getWeight(i,j)<<endl;
						}
						else
						{
							cout<<"Edge("<<i<<","<<j<<"): Space:"<<root_g->getWeight(i,j)<<endl;
						}
					}
				
				}
			}
#endif

		CutTreeNode* root = new CutTreeNode(root_g);
		splitNode(root, 0);
		int tagCount = setNewTags(root, newTags);

		for(int i = 0; i < root_g->vertex_count; i++)
		{
			delete[] isKeyEdge[i];
		}
		delete[] isKeyEdge;
		isKeyEdge = NULL;
		delete root;


		return tagCount;
	}

	int BlockGraph::setNewTags(CutTreeNode* root, std::vector<int>& newTags)
	{
		int tagCount = 0;
		vector<int> tags(newTags.size());
		root->setNewTags(tags, tagCount);
		vector<int> proTags(tagCount, -1);
		int count = 0;
		for(int i = 0; i < newTags.size(); i++)
		{
			if(proTags[tags[i]] < 0)
			{
				proTags[tags[i]] = count++;
			}
			newTags[i] = proTags[tags[i]];
		}
		return count;
	}

	bool BlockGraph::KeyEdgeJudgement(int block_a, int block_b)
	{
		if(blocks[block_a].coff.seged_id >= 0 && blocks[block_b].coff.seged_id >= 0 && blocks[block_a].coff.seged_id != blocks[block_b].coff.seged_id) return false;

		double a2b = spaceDistanceCache[block_a][block_b] * blocks[block_a].inliers;
	    //double b2a = spaceDistanceCache[block_b][block_a] * blocks[block_b].inliers;

		return a2b < 2500;
		//return min(a2b, b2a) < 2500;
	}

	bool Compare(pair<int, double>& a, pair<int, double>& b)
	{
		return a.second < b.second;
	}

	double BlockGraph::initialEdge(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, Graph& graph, double disthre)
	{
		int keep_edge = blocks.size() / 10 + 1;
		//vector<double> m_dist(blocks.size(), 0);
		vector<int> edgeCount(blocks.size(), 0);
		for(int i = 0; i < tags.size(); i++)
		{
			int tagi = tags[i];
			if(tagi >= 0 && tagi < blocks.size())
			{
				
				int px = i % cloud->width;
				int py = i / cloud->width;
				PointT& pp = cloud->at(i);

				//neighbour with blocks
				for(int jx = px-1;jx <= px+1;jx++)
				{
					for(int jy = py-1;jy <= py+1;jy++)
					{
						int nid = jy*cloud->width+jx;
						if(jx >= 0 && jx < cloud->width && jy >= 0 && jy < cloud->height && !ISNAN(normals->at(nid).normal_x))
						{
						   PointT& pn = cloud->at(nid);
						   double disp = NORM2(pp.x-pn.x,pp.y-pn.y,pp.z-pn.z);
						   if(disp < disthre)
						   {
							    int tagj = tags[nid];
								if(tagj >= 0 && tagj < blocks.size() && tagj != tagi)
								{
									if(!graph.haveEdge(tagi, tagj))
									{
										graph.addEdge(tagi, tagj, blockEdgeWeight(tagi, tagj, true));
										//graph.addEdge(tagi, tagj, blockEdgeWeightOneDirection(tagi, tagj, true), blockEdgeWeightOneDirection(tagj, tagi, true));
										edgeCount[tagi]++;
										edgeCount[tagj]++;
									}
								}
						   }
						}
					}
				}
			}
		}

		cout<<"size:"<<blocks.size()<<endl;
		for(int i = 0; i < blocks.size(); i++)
		{
			//if(edgeCount[i] < keep_edge)
			//{
				vector<pair<int,double>> dists;
				for(int j = 0; j < blocks.size(); j++)
				{
					if(j == i) continue;
					//pair<int,double> p(j, blockSpaceDistance(i, j));
					pair<int,double> p(j, blockSpaceDistanceOneDirection(i, j));
					dists.push_back(p);
				}
				sort(dists.begin(), dists.end(), Compare);
				
				for(int k = 0; k < 1; k++)
				{
					int j = dists[k].first;
					if(KeyEdgeJudgement(i, j))
					{
						isKeyEdge[i][j] = true;
						isKeyEdge[j][i] = true;
						
						if(!graph.haveEdge(i, j))
						{
							graph.addEdge(i, j, blockEdgeWeight(i, j, false));
							edgeCount[i]++;
							edgeCount[j]++;

						}
						
					}
				}

				
				//if(edgeCount[i] <= 0)
				//{
					int iter = 0;
					while(edgeCount[i] < keep_edge)
					{
						int j = dists[iter].first;
						if(!graph.haveEdge(i, j))
						{
							graph.addEdge(i, j, blockEdgeWeight(i, j, false));
							//graph.addEdge(i, j, blockEdgeWeightOneDirection(i, j, false), blockEdgeWeightOneDirection(j, i, false));
							edgeCount[i]++;
							edgeCount[j]++;

						}
						iter++;
					}
				//}
			//}
		}

		/*
		for(int i = 0; i < blocks.size(); i++)
		{
			for(int j = 0; j < virtual_blocks.size(); j++)
			{
				double ins = virtual_blocks[j].intersection(blocks[i]);
				if(ins > 0.8)
				{
					int vid = blocks.size() + j;
					graph.addEdge(i, vid, blockEdgeWeight(i, vid, true));
					isKeyEdge[i][vid] = true;
					isKeyEdge[vid][i] = true;
				}
			}
		}
		*/

		return 0;
	}
	

	double BlockGraph::blockSpaceDistance(int a, int b)
	{
		//if(b >= blocks.size()) return spaceDistanceCache[a][b];
		return min(spaceDistanceCache[a][b], spaceDistanceCache[b][a]);
	}

	double BlockGraph::blockSpaceDistanceOneDirection(int a, int b)
	{
		return spaceDistanceCache[a][b];
	}

	double BlockGraph::blockColorDistance(int a, int b)
	{

		double disE = NORM2(blocks[a].ELab[0] - blocks[b].ELab[0], blocks[a].ELab[1] - blocks[b].ELab[1], blocks[a].ELab[2] - blocks[b].ELab[2]);
		double disV = NORM2(blocks[a].VLab[0] - blocks[b].VLab[0], blocks[a].VLab[1] - blocks[b].VLab[1], blocks[a].VLab[2] - blocks[b].VLab[2]);

		return sqrt(disE*disE + disV * disV);
	}

	double BlockGraph::blockEdgeWeight(int a, int b, bool neighbour)
	{
		double space_dis = blockSpaceDistance(a, b);
		double space_weight = 1 / (space_dis + 0.0001);
		//if(!neighbour) space_weight *= 0.3;
		
		return space_weight;  // + color_weight;
	}

	double BlockGraph::blockEdgeWeightOneDirection(int a, int b, bool neighbour)
	{
		double space_dis = blockSpaceDistanceOneDirection(a, b);
		double space_weight = 1 / (space_dis + 0.0001);
		//if(!neighbour) space_weight *= 0.3;
		
		return space_weight;  // + color_weight;
	}

	int Segmentation::nextSegment(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, vector<int>& initTags,
		                          std::vector<Block>& blocks, std::vector<Block>& seged_blocks)
	{
		std::vector<Plane3D> planes;
		double disthre = 0.05;
		double angthre = 0.75;
		dbplane_segment(cloud, normals, 1000, disthre, angthre, disthre, planes, tags);
		int tagCount = planes.size();
		vector<int> newTags_(tagCount);
		vector<int> foundTags(seged_blocks.size(), -1);
		if(seged_blocks.size() > 0)
		{
			tagCount = calSegedTags(cloud, normals, initTags, tags, newTags_, foundTags, planes, seged_blocks);
			
		}else
		{
			for(int i = 0; i < tagCount; i++)
			{
				newTags_[i] = i;
			}
		}

		vector<coffCal> coffs;
		calNewCoffandTag(cloud, normals, tags, tagCount, newTags_, foundTags, coffs);

		//int static loop = 1;
		std::vector<Block> virtual_blocks;
		while(true)
	    //for(int i = 0; i < 0; i++)
		{
			vector<int> newTags(tagCount);
			BlockGraph graph(cloud, normals, tags, coffs, virtual_blocks);
			if(tagCount <= 1)
			{
				graph.getBlocks(blocks);
				break;
			}

			double newTagCount = graph.graphCut(cloud, normals, tags, newTags, disthre);

		    if(newTagCount == tagCount)
			{
				//loop--;
				graph.getBlocks(blocks);
				break;
			}
			tagCount = newTagCount;
			vector<coffCal> newCoffs;
			calNewCoffandTag(cloud, normals, tags, tagCount, newTags, coffs, newCoffs);
			coffs = newCoffs;
		}

		

		
		return tagCount;
		
		//return coffs.size();
	}

	void Segmentation::calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals,
			                  std::vector<int>& tags, int size, std::vector<int>& newTags, 
							  std::vector<int>& foundTags, std::vector<coffCal>& newCoffs)
	{
		for(int i = 0; i < size; i++)
		{
			newCoffs.push_back(coffCal());
		}
		for(int i = 0; i < foundTags.size(); i++)
		{
			if(foundTags[i] >= 0)
			{
				newCoffs[foundTags[i]].seged_id = i;
			}
		}
		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < newTags.size())
			{
				tags[i] = newTags[tags[i]];
				if(tags[i] >= 0)
				{
				   newCoffs[tags[i]].addPoint(cloud->at(i));
				}
			}
		}

		for(int i = 0; i < size; i++)
		{
			Plane3D coff = newCoffs[i].getCoff();
			newCoffs[i].calAxis();
		}
	}

	void Segmentation::calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals,
			                  std::vector<int>& tags, int size, std::vector<int>& newTags,
							  std::vector<coffCal>& oldCoffs, std::vector<coffCal>& newCoffs)
	{
		for(int i = 0; i < size; i++)
		{
			newCoffs.push_back(coffCal());
		}
		for(int i = 0; i < oldCoffs.size(); i++)
		{
			int seged_id = oldCoffs[i].seged_id;
			if(seged_id >= 0)
			{
				newCoffs[newTags[i]].seged_id = seged_id;
			}
		}
		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < newTags.size())
			{
				tags[i] = newTags[tags[i]];
				if(tags[i] >= 0)
				{
				   newCoffs[tags[i]].addPoint(cloud->at(i));
				}
			}
		}

		for(int i = 0; i < size; i++)
		{
			Plane3D coff = newCoffs[i].getCoff();
			newCoffs[i].calAxis();
		}
	}

	int Segmentation::calSegedTags(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& initTags,
			             std::vector<int>& tags, std::vector<int>& newTags, std::vector<int>& foundTags,
			             std::vector<Plane3D>& planes, std::vector<Block>& seged_blocks)
	{
		double inside_threshold = 0.7;
		int** inside = new int*[planes.size()];
		for(int i = 0; i <planes.size(); i++)
		{
			inside[i] = new int[seged_blocks.size()];
			for(int j = 0; j < seged_blocks.size(); j++)
			{
				inside[i][j] = 0;
			}
		}
		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < planes.size())
			{
				if(initTags[i] >= 0) inside[tags[i]][initTags[i]]++;
				else
				{
					for(int j = 0; j< seged_blocks.size(); j++)
					{
						if(seged_blocks[j].isInBlock(cloud->at(i))) inside[tags[i]][j]++;
					}
				}
			}
		}
		
		int foundCount = 0;
		for(int i = 0;i < planes.size(); i++)
		{
			double maxr = 0;
			int maxj = -1;
			for(int j = 0; j < seged_blocks.size(); j++)
			{
				double rate = static_cast<double>(inside[i][j]) / planes[i].inliers;
				if(rate > maxr)
			    {
					maxr = rate;
					maxj = j;
				}
			}

			if(maxr >= inside_threshold)
			{
				if(foundTags[maxj] < 0) foundTags[maxj] = foundCount++;
				newTags[i] = foundTags[maxj];
			}
			else if(planes[i].inliers < 100)
			{
				newTags[i] = -1;
			}
			else
			{
				newTags[i] = foundCount++;
			}
		}
		return foundCount;
	}
}
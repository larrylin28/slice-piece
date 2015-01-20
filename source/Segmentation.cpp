#include "Segmentation.h"
#include "Util.h"
#include <map>
#include <queue>

using namespace std;

#define e 2.718281828459
#define DIV 0.003125



namespace Rgbd
{
	
	Block::Block():EL(0),Ea(0),Eb(0),VL(0),Va(0),Vb(0),inliers(0)
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

	void Block::getProjectCoord(int x, int y, int z, double coord[3])
	{
		 coord[0] = local(0,0)*x+local(0,1)*y+local(0,2)*z + local(0,3)*1;
	     coord[1] = local(1,0)*x+local(1,1)*y+local(1,2)*z + local(1,3)*1;
	     coord[2] = local(2,0)*x+local(2,1)*y+local(2,2)*z + local(2,3)*1;
	}

	inline double out_range(double l, double r, double x)
	{
		if(x < l) return (l - x);
		else if(x > r) return (x - r);
		else return 0;
	}

	double Block::getDistance(PointT& p)
	{
		double coord[3];
		getProjectCoord(p, coord);
		double ox = out_range(bound[0], bound[3], coord[0]);
		double oy = out_range(bound[1], bound[4], coord[1]);
		double oz = out_range(bound[2], bound[5], coord[2]);
		
		double v1 = (bound[3] - bound[0]) * (bound[4] - bound[1]) * (bound[5] - bound[2]);
		double v2 = (bound[3] - bound[0] + ox) * (bound[4] - bound[1] + oy) * (bound[5] - bound[2] + oz);
		return sqrt(v2 - v1);
		//return (v2 - v1) / v1;
	}

	double Block::getDistance(Block& b)
	{
		
		return 0;
	}
                     
	void Block::addInlier(PointT& p, pcl::Normal& n)
	{	
		//3d space
		//project coord
		double coord[3];
		getProjectCoord(p, coord);
		if(inliers > 0)
		{
			bound[0] = min(bound[0], coord[0]);
			bound[1] = min(bound[1], coord[1]);
			bound[2] = min(bound[2], coord[2]);
			bound[3] = max(bound[3], coord[0]);
			bound[4] = max(bound[4], coord[1]);
			bound[5] = max(bound[5], coord[2]);
		}
		else
		{
			bound[0] = coord[0];
			bound[1] = coord[1];
			bound[2] = coord[2];
			bound[3] = coord[0];
			bound[4] = coord[1];
			bound[5] = coord[2];
		}

		//color space
		double L, a, b;
		RGB2Lab(p.r, p.g, p.b, L, a, b);
		meanVar(EL, VL, L, inliers);
		meanVar(Ea, Va, a, inliers);
		meanVar(Eb, Vb, b, inliers);

		inliers++;
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

	BlockGraph::BlockGraph(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, std::vector<Plane3D>& coffs)
		:blocks(coffs.size()),isKeyEdge(NULL)
	{
		for(int i = 0; i < blocks.size(); i++)
		{
			blocks[i].coff = coffs[i];

			//cal rotate matrix
			double before[3];
		    double after[3];
			before[0] = coffs[i].a;
			before[1] = coffs[i].b;
			before[2] = coffs[i].c;

			after[0] = 0;
			after[1] = 0;
			after[2] = 1;

			double rotate[4][4];
			Calculation4d(before, after, rotate);

			for(int ri = 0; ri < 4; ri++)
				for(int rj = 0; rj < 4; rj++)
		            blocks[i].local(ri, rj) = rotate[ri][rj];
		}

		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < blocks.size())
			{
				blocks[tags[i]].addInlier(cloud->at(i), normals->at(i));
			}
		}
		for(int i = 0; i < blocks.size(); i++)
		{
			blocks[i].updateCore();
		}

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
				maxDis = max(maxDis, spaceDistanceCache[i][j]);
				minDis = min(minDis, spaceDistanceCache[i][j]);
			}
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

	double BlockGraph::toleration(CutTreeNode* root)
	{
		double sumVolume = 0;
		double sumInliers = 0;
		for(int i = 0; i < root->data->vertex_count; i++)
		{
			Block &b = blocks[root->data->getVertex(i)];
			sumVolume += (b.bound[3] - b.bound[0]) * (b.bound[4] - b.bound[1]) * (b.bound[5] - b.bound[2]);
			sumInliers  += b.inliers;
		}
		double torerateDis  = sqrt(sumVolume);
		double torerateWeight = 1 / (torerateDis + 0.0001);
		return torerateWeight;
	}

	void BlockGraph::splitNode(CutTreeNode* root, double maxCutValue)
	{
		if(root->data->vertex_count <= 1) return;

		vector<int> setS;
		vector<int> setT;
		double cutvalue = root->data->minCut(setS, setT);
		//double tolerateCut = toleration(root);
		bool isConnected = connectedWithKeyEdge(root);

		//if(cutvalue > maxCutValue)
		//{
			if(isConnected) return;
		//}

		double maxC = max(cutvalue, maxCutValue);

		//ouput cut
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

		root->split(setS, setT);
		if(root->left) splitNode(root->left, maxC);
		if(root->right) splitNode(root->right, maxC);
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

		CutTreeNode* root = new CutTreeNode(root_g);
		splitNode(root, 0);
		int tagCount = 0;
		root->setNewTags(newTags, tagCount);

		for(int i = 0; i < root_g->vertex_count; i++)
		{
			delete[] isKeyEdge[i];
		}
		delete[] isKeyEdge;
		isKeyEdge = NULL;
		delete root;


		return tagCount;
	}

	bool BlockGraph::KeyEdgeJudgement(int block_a, int block_b)
	{
		double a2b = spaceDistanceCache[block_a][block_b] * blocks[block_a].inliers;
		double b2a = spaceDistanceCache[block_b][block_a] * blocks[block_b].inliers;

		return min(a2b, b2a) < 2500;
	}

	bool Compare(pair<int, double>& a, pair<int, double>& b)
	{
		return a.second < b.second;
	}

	double BlockGraph::initialEdge(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags, Graph& graph, double disthre)
	{
		int keep_edge = blocks.size() / 10 + 1;
		double m_edge = 1000000;
		vector<int> edgeCount(blocks.size(), 0);
		for(int i = 0; i < tags.size(); i++)
		{
			int tagi = tags[i];
			if(tagi >= 0 && tagi < blocks.size())
			{
				int px = i % cloud->width;
				int py = i / cloud->width;
				PointT& pp = cloud->at(i);
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
										graph.addEdge(tagi, tagj, blockEdgeWeight(tagi, tagj));
										m_edge = min(m_edge, graph.getWeight(tagi,tagj));
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
					}
				}
				if(edgeCount[i] <= 0)
				{
					int iter = 0;
					while(edgeCount[i] < keep_edge)
					{
						int j = dists[iter].first;
						if(!graph.haveEdge(i, j))
						{
							graph.addEdge(i, j, blockEdgeWeight(i, j) * 0.5);
							m_edge = min(m_edge, graph.getWeight(i,j));
							edgeCount[i]++;
							edgeCount[j]++;
						}
						iter++;
					}
				}
			//}
		}
		return m_edge;
	}
	

	double BlockGraph::blockSpaceDistance(int a, int b)
	{
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

	double BlockGraph::blockEdgeWeight(int a, int b)
	{
		double space_dis = blockSpaceDistance(a, b);

		double space_weight = 1 / (space_dis + 0.0001);

		return space_weight;  // + color_weight;
	}

	int Segmentation::nextSegment(PointCloud::Ptr cloud, PointNormal::Ptr normals, vector<int>& tags)
	{
		std::vector<Plane3D> coffs;
		double disthre = 0.05;
		double angthre = 0.75;
		dbplane_segment(cloud, normals, 1000, disthre, angthre, disthre, coffs, tags);
		int tagCount = coffs.size();

		while(true)
		//for(int i = 0; i < 1; i++)
		{
			if(tagCount <= 1) break;
			vector<int> newTags(tagCount);
			BlockGraph graph(cloud, normals, tags, coffs);
			double newTagCount = graph.graphCut(cloud, normals, tags, newTags, disthre);

		    if(newTagCount == tagCount) break;
			tagCount = newTagCount;
			vector<Plane3D> newCoffs;
			calNewCoffandTag(cloud, normals, tags, tagCount, newTags, newCoffs);
			coffs = newCoffs;
		}
		return tagCount;
		
		//return coffs.size();
	}

	void Segmentation::calNewCoffandTag(PointCloud::Ptr cloud, PointNormal::Ptr normals, 
			                  std::vector<int>& tags, int size, std::vector<int>& newTags, std::vector<Plane3D>& newCoffs)
	{
		vector<coffCal> cal(size, coffCal());
		for(int i = 0; i < tags.size(); i++)
		{
			if(tags[i] >= 0 && tags[i] < newTags.size())
			{
				tags[i] = newTags[tags[i]];
				cal[tags[i]].addPoint(cloud->at(i));
			}
		}
		for(int i = 0; i < size; i++)
		{
			newCoffs.push_back(cal[i].getCoff());
		}
	}

}
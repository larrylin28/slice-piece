#include "MinCut.h"
#include <queue>
#include <iostream>

using namespace std;

namespace Rgbd
{
	double maxWeight = 1000000.0;

	void Graph::addEdge(int v1, int v2, double weight)
	{
		edge[v1][v2] += weight;
		edge[v2][v1] += weight;
	}

	void Graph::swapVertexes(int v1, int v2)
	{
		if(v1 == v2) return;
		swap(vertex[v1], vertex[v2]);
		for(int i = 0; i < vertex_count; i++)
		{
			swap(edge[v1][i], edge[v2][i]);
		}
		for(int i = 0; i < vertex_count; i++)
		{
			swap(edge[i][v1], edge[i][v2]);
		}
	}

	double Graph::minCut(vector<int>& setS, vector<int>& setT)
	{
		double minC = maxWeight * vertex_count;
		bool *visited = new bool[vertex_count];
		int *father = new int[vertex_count];
		double** mp = new double*[vertex_count];
		for(int i = 0; i < vertex_count; i++)
		{
			mp[i] = new double[vertex_count];
		}

		for(int tid = 1; tid < vertex_count; tid++)
		{
			swapVertexes(tid, vertex_count - 1);
			double maxf = maxFlow(visited, father, mp);
			double cutvalue = minCutValue(visited, father);
			if(cutvalue < minC)
			{
				minC = cutvalue;
				setS.clear();
				setT.clear();
				for(int i = 0; i < vertex_count; i++)
				{
					int realid = i;
					if(realid == tid) realid = vertex_count - 1;
					else if(realid ==  vertex_count - 1) realid = tid;
					if(visited[i]) setS.push_back(realid);
					else setT.push_back(realid);
				}
			}
			swapVertexes(tid, vertex_count - 1);
		}

		for(int i = 0; i < vertex_count; i++) delete[] mp[i];
		delete[] mp;
		delete[] visited;
		delete[] father;
		
		return minC;
	}

	double Graph::maxFlow(bool* visited, int* father, double** mp)
	{
		double maxflow = 0;
		for(int i = 0; i < vertex_count; i++)
		{
			for(int j = 0; j < vertex_count; j++)
			{
				mp[i][j] = edge[i][j];
			}
		}
		while(true)
		{
			//һ�δ�ѭ��,�ҵ�һ�����ܵ�����·��
			queue <int> q;
			for(int i = 0; i < vertex_count; i++)
			{
				visited[i] = false;
				father[i] = -1;
			}
			int now;
			visited[0] = true;
			q.push(0);
			while(!q.empty())//�������
			{
				now = q.front();
				q.pop();
				if(now == vertex_count-1) break;
				for(int i = 0; i < vertex_count; i++)
				{
					//ÿ�θ��׽ڵ㶼Ҫ����,Ȩֵ��Ϊ0�ı߾Ͳ�����.
					if(mp[now][i] > 0 && !visited[i])
					{
						father[i] = now;
						visited[i] = true;
						q.push(i);
					}
				}
			}
			//���ܵ�����·��������
			if(!visited[vertex_count-1]) break;
			int u;
			double min = maxWeight;
			for(u = vertex_count-1; u; u = father[u])//�ҳ�Ȩֵ��С�ı�
			{
				if(mp[father[u]][u] < min)
					min = mp[father[u]][u];
			}
			//��ȥ��СȨֵ
			for(u = vertex_count-1; u > 0; u = father[u])
			{
				//ǰ�򻡼�ȥ
				mp[father[u]][u] -= min;
				//���򻡼���
				//����Բ��,��仰�ؼ�
				mp[u][father[u]] += min;
			}
			//��ǰ����·�����ӵ���
			maxflow += min;
		}
		return maxflow;
	}

	double Graph::minCutValue(bool* visited, int* father)
	{
		double cutValue = 0;
		for(int i = 0; i < vertex_count; i++)
		{
			for(int j = i + 1; j < vertex_count; j++)
			{
				if(visited[i] != visited[j])
				{
					cutValue += edge[i][j];
				}
			}
		}
		return cutValue;
	}

	void CutTreeNode::split(vector<int>& setS, vector<int>& setT)
	{
		Graph* left_g = new Graph(setS.size());
		splitGraph(*data, setS, *left_g);
		left = new CutTreeNode(left_g);
		//left->split();

		Graph* right_g = new Graph(setT.size());
		splitGraph(*data, setT, *right_g);
		right = new CutTreeNode(right_g);
		//right->split();
	}

	void CutTreeNode::setNewTags(std::vector<int>& newTags, int& tagCount)
	{
		if(!left && ! right)
		{
			for(int i = 0; i < data->vertex_count; i++)
			{
				newTags[data->getVertex(i)] = tagCount;
			}
			tagCount++;
		}
		
		if(left) left->setNewTags(newTags, tagCount);
		if(right) right->setNewTags(newTags, tagCount);
	}

	void splitGraph(Graph& parent, std::vector<int>& setS, Graph& s)
	{
		for(int i = 0; i < setS.size(); i++)
		{
			s.setVertex(i, parent.getVertex(setS[i]));
		}
		for(int i = 0; i < setS.size(); i++)
		{
			for(int j = i + 1; j < setS.size(); j++)
			{
				s.addEdge(i, j, parent.getWeight(setS[i], setS[j]));
			}
		}
	}

}
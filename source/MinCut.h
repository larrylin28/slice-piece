#ifndef MIN_CUT_H
#define MIN_CUT_H
#include <vector>

namespace Rgbd
{
	
	class Graph
	{
	public:
		int vertex_count;
	public:
		Graph(int c):vertex_count(c)
		{
			vertex = new int[vertex_count];
			edge = new double*[vertex_count];
			for(int i = 0; i < vertex_count; i++)
			{
				vertex[i] = i;
				edge[i] = new double[vertex_count];
				for(int j = 0; j < vertex_count; j++)
				{
					edge[i][j] = 0;
				}
			}
			
		}
		~Graph()
		{
			for(int i = 0; i < vertex_count; i++) delete[] edge[i];
			delete[] edge;
			delete[] vertex;
		}

		bool haveEdge(int v1, int v2)
		{
			return edge[v1][v2] > 0;
		}

		double getWeight(int v1, int v2)
		{
			return edge[v1][v2];
		}

		void setVertex(int v, int id)
		{
			vertex[v] = id;
		}

		int getVertex(int v)
		{
			return vertex[v];
		}

		void addEdge(int v1, int v2, double weight);
		void swapVertexes(int v1, int v2);
		double  minCut(std::vector<int>& setS, std::vector<int>& setT);
	private:
	    double** edge;
		int* vertex;
		double maxFlow(bool* visited, int* father, double** mp);
		double minCutValue(bool* visited, int* father);
	};

	struct CutTreeNode
	{
		Graph* data;
		CutTreeNode* left;
		CutTreeNode* right;

		CutTreeNode(Graph* g):data(g),left(NULL),right(NULL)
		{
		}

		~CutTreeNode()
		{
			delete left;
			delete right;
			delete data;
		}

		void split(std::vector<int>& setS, std::vector<int>& setT);
		void setNewTags(std::vector<int>& newTag, int& tagCount);
	};

	void splitGraph(Graph& parent, std::vector<int>& setS, Graph& s);
}
#endif // MIN_CUT_H

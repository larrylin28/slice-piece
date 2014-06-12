#ifndef GL_WINDOWS_H_
#define GL_WINDOWS_H_

#include <windows.h>		// Header File For Windows
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <vector>
#include <map>
#include <set>
#include <list>
#include <utility>
#include <tuple>
#include <pcl/point_types.h>
// -------------------- OpenMesh   
//#include <OpenMesh/Core/IO/MeshIO.hh>   
//#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh> 

//typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

class GL_Mesh;
struct GLMesh{
	int vertexes[3];
	inline GLMesh(int v1,int v2,int v3)
	{
		vertexes[0] = v1;
		vertexes[1] = v2;
		vertexes[2] = v3;
	}
};

struct GLVertex{
	double v[6];  //x,y,z,nx,ny,nz
	double dist;
	int tran;
	double up;
	double upxyz[3];
	double upcount;
	int arc;
	double degree;
	inline GLVertex():degree(0),arc(0),up(0),upcount(0)
	{
		upxyz[0] = 0;
		upxyz[1] = 0;
		upxyz[2] = 0;
	}

	void DoUpdate();

	void LaplacianUpdate(double v1[6]);

    void addDegree(double v1[6], double v2[6])
	{
		double cosLen1 = ((v1[0]-v[0])*v[3]+(v1[1]-v[1])*v[4]+(v1[2]-v[2])*v[5]);
		double t1[3];
		t1[0] = (v1[0]-v[0]) - cosLen1*v[3];
		t1[1] = (v1[1]-v[1]) - cosLen1*v[4];
		t1[2] = (v1[2]-v[2]) - cosLen1*v[5];

		double cosLen2 = ((v2[0]-v[0])*v[3]+(v2[1]-v[1])*v[4]+(v2[2]-v[2])*v[5]);
		double t2[3];
		t2[0] = (v2[0]-v[0]) - cosLen2*v[3];
		t2[1] = (v2[1]-v[1]) - cosLen2*v[4];
		t2[2] = (v2[2]-v[2]) - cosLen2*v[5];

		double theta = (t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2])/sqrt(t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2])/sqrt(t2[0]*t2[0] + t2[1]*t2[1] + t2[2]*t2[2]);
		double ac = acos(theta);
		degree += ac;
	}
	//std::vector<int> meshes;
};


struct GLMergeVertex{
	int vid;  
	int weight;
};

class GL_Mesh{
public:
  static std::vector<GLVertex>* gl_vertexes;
  static std::map<std::tuple<int,int,int>,GLMergeVertex>* gl_vertex_map;
  static std::map<std::tuple<int,int,int>,GLMesh>* gl_meshes;
  static std::map<std::pair<int,int>,int>* gl_edges;
  static std::vector<Eigen::Matrix4f>* gl_trans;
  //static MyMesh o_meshes;  

  static void updateVertexes()
  {
	  for(std::vector<GLVertex>::iterator it = GL_Mesh::gl_vertexes->begin();it != gl_vertexes->end();it++)
	  {
		  it->DoUpdate();
	  }
  }

  static void addEdge(int a,int b)
  {
	  if(a < b) (*GL_Mesh::gl_edges)[std::pair<int,int>(a,b)]++;
	  else if(a > b) (*GL_Mesh::gl_edges)[std::pair<int,int>(b,a)]++;
  }

  static bool findEdge(int a,int b)
  {
	std::map<std::pair<int,int>,int>::iterator it;
	if(a < b) it = GL_Mesh::gl_edges->find(std::pair<int,int>(a,b));
	else if(a > b) it = GL_Mesh::gl_edges->find(std::pair<int,int>(b,a));
    
	if(it !=  GL_Mesh::gl_edges->end()) return true;
	else return false;
  }

  static bool findMesh(int a,int b,int c)
  {
	 if (a>b) std::swap(a, b);
	if (b>c) 
	{ 
		std::swap(b, c);
		if (a > b) std::swap(a, b);
	}
	std::map<std::tuple<int,int,int>,GLMesh>::iterator it = GL_Mesh::gl_meshes->find(std::tuple<int,int,int>(a,b,c));
    
	if(it !=  GL_Mesh::gl_meshes->end()) return true;
	else return false;
  }

  static int findEdgeDegree(int a,int b)
  {
	std::map<std::pair<int,int>,int>::iterator it;
	if(a < b) it = GL_Mesh::gl_edges->find(std::pair<int,int>(a,b));
	else if(a > b) it = GL_Mesh::gl_edges->find(std::pair<int,int>(b,a));
	else return 0;
    
	if(it !=  GL_Mesh::gl_edges->end()) return it->second;
	else return 0;
  }
};


#endif
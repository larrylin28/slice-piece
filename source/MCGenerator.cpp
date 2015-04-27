#include "MCGenerator.h"
#include "TriangleMesh.h"

namespace Rgbd
{
	MCGenerator::MCGenerator(TsdfModel* model) : sdf(model)
	{
	}

	MCGenerator::~MCGenerator()
	{
	}

	void MCGenerator::GenerateMeshes(float ix, float iy, float iz, int Xlen, int Ylen, int Zlen, float div,	GLMesh* mesh){
		CIsoSurface<float> *gen = new CIsoSurface<float>();
		long size = Xlen*Ylen*Zlen;
		float* absrate = new float[size];
		sdf->getSample(div, ix, iy, iz, Xlen, Ylen, Zlen, absrate);
		gen->GenerateSurface(absrate,50,Xlen-1,Ylen-1,Zlen-1,div,div,div);

		Eigen::Matrix4f& tr = sdf->local_tran;

		for(int i = 0; i < gen->m_nVertices; i++)
		{
			GLVertex v;
			v.v[0] = gen->m_ppt3dVertices[i][0] + ix;
			v.v[1] = gen->m_ppt3dVertices[i][1] + iy;
			v.v[2] = gen->m_ppt3dVertices[i][2] + iz;
			v.v[3] = gen->m_pvec3dNormals[i][0];
			v.v[4] = gen->m_pvec3dNormals[i][1];
			v.v[5] = gen->m_pvec3dNormals[i][2];

			GLVertex tv;
			tv.v[0] = tr(0,0)*v.v[0] + tr(0,1)*v.v[1] + tr(0,2)*v.v[2] + tr(0,3);
			tv.v[1] = tr(1,0)*v.v[0] + tr(1,1)*v.v[1] + tr(1,2)*v.v[2] + tr(1,3);
			tv.v[2] = tr(2,0)*v.v[0] + tr(2,1)*v.v[1] + tr(2,2)*v.v[2] + tr(2,3);
			tv.v[3] = tr(0,0)*v.v[3] + tr(0,1)*v.v[4] + tr(0,2)*v.v[5] + tr(0,3);
			tv.v[4] = tr(1,0)*v.v[3] + tr(1,1)*v.v[4] + tr(1,2)*v.v[5] + tr(1,3);
			tv.v[5] = tr(2,0)*v.v[3] + tr(2,1)*v.v[4] + tr(2,2)*v.v[5] + tr(2,3);
			mesh->gl_vertexes.push_back(tv);
		}

		for(int i = 0; i < gen->m_nTriangles; i++)
		{
			TryaddTriangle(gen->m_piTriangleIndices[i*3],gen->m_piTriangleIndices[i*3 + 1],gen->m_piTriangleIndices[i*3  + 2], mesh);
		}
	}

	void MCGenerator::TryaddTriangle(int v1,int v2,int v3, GLMesh* mesh)
	{
		GLTriangle m1(v1,v2,v3);	   
		int a = v1;
		int b = v2;
		int c = v3;
		if (a>b) std::swap(a, b);
		if (b>c) 
		{ 
			std::swap(b, c);
			if (a > b) std::swap(a, b);
		}
		if(mesh->gl_meshes.find(std::tuple<int,int,int>(a,b,c)) != mesh->gl_meshes.end()) return;
		std::pair<std::tuple<int,int,int>,GLTriangle> ins(std::tuple<int,int,int>(a,b,c),m1);		
		mesh->gl_meshes.insert(ins);
			
		mesh->addEdge(v1,v2);
		mesh->addEdge(v2,v3);
		mesh->addEdge(v3,v1);	
	}
}
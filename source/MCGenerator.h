#ifndef MCGENERATOR_H
#define MCGENERATOR_H
//Marching Cube Generator

#include "TsdfModel.h"
#include "CIsoSurface.h"

namespace Rgbd
{
	class GLMesh;
	class GLTriangle;
	class GLVertex;

	class MCGenerator
	{
	public:
		//Sdfmoel
		TsdfModel* sdf;

	public:

		MCGenerator(TsdfModel* model);
		~MCGenerator();

		void GenerateMeshes(float ix, float iy, float iz, int Xlen, int Ylen, int Zlen, float div, GLMesh* mesh);
		void TryaddTriangle(int v1,int v2,int v3, GLMesh* mesh);

	};
}

#endif // MCGENERATOR_H

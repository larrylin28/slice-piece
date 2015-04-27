
#include "RGBDReader.h"
#include "TsdfModel.h"
#include "HashTsdfModel.h"
#include "MultiTsdfModel.h"
#include "SNPGenerator.h"
#include "MCGenerator.h"
#include "Segmentation.h"
#include "NeHeGL.h"
#include "TriangleMesh.h"
#include <fstream>



using namespace Rgbd;


void outputObj(GLMesh* mesh)
{
	
	std::ofstream output("out.obj");

	for(std::vector<GLVertex>::iterator it = mesh->gl_vertexes.begin();it != mesh->gl_vertexes.end();it++)
	{
		GLVertex v = *it;
		output<<"v "<<v.v[0]<<" "<<v.v[1]<<" "<<v.v[2]<<std::endl;
	}
	for(std::vector<GLVertex>::iterator it = mesh->gl_vertexes.begin();it != mesh->gl_vertexes.end();it++)
	{
		GLVertex v = *it;
		output<<"vn "<<v.v[3]<<" "<<v.v[4]<<" "<<v.v[5]<<std::endl;
	}
	for(std::map<std::tuple<int,int,int>,GLTriangle>::iterator it =mesh->gl_meshes.begin();it != mesh->gl_meshes.end();it++)
	{
		GLTriangle m = it->second;
		output<<"f "<<m.vertexes[0]+1<<"//"<<m.vertexes[0]+1<<" "<<m.vertexes[1]+1<<"//"<<m.vertexes[1]+1<<" "<<m.vertexes[2]+1<<"//"<<m.vertexes[2]+1<<std::endl;
	}
	output.close();
}

void runReconstruction()
{
   RgbdReader reader;
   reader.Read("D:/PCL/record.ONI", 0, 10, 10, 0.1);
   //reader.Read("D:/PCL/record1.ONI",20, 10, 43, 0.1);
   //reader.Read("D:/PCL/DATA/1.ONI",200, 10, 10, 0.1);
   //reader.Read("F:/1.ONI",0, 1, 10, 0.1);
   reader.ShowPointCoud();

   
   MultiTsdfModel tsdf;
   int frameCount = reader.getFrameCount();
   for(int i = 0; i < frameCount; i++)
   {
	   reader.putDataToSdf(&tsdf, i);
   }

   tsdf.updateTags(&reader);
   tsdf.showFinalSegment(&reader);

   /*   
   tsdf.updateTags(&reader);

   GLMesh* mesh = new GLMesh();
   for(int i = 0; i < tsdf.blockCount(); i++)
   {
	   tsdf.MeshGenForBlock(i, &reader, mesh, 0.1, 0.1);
   }
   */


  
   

   /*
   double devide = 0.01;
   Eigen::Matrix4f zero;
   zero << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;
   TsdfModel tsdf(zero, 1.2, 1, 4.0, -1.5, -1.0, 0.0, devide, devide * 10, devide / 10, false);
   
    SNPGenerator gen(&tsdf, devide);
    //MCGenerator gen(&tsdf); 
	GLMesh* mesh = new GLMesh();

   int frameCount = reader.getFrameCount();
   Plane3D current_coff;
   for(int i = 0; i < frameCount; i++)
   {
	   reader.putDataToSdf(&tsdf, i);
	   gen.castNextMeshes(reader.getTrans(i));
   }
   
   //gen.GenerateMeshes(-1.5, -1.0, 0.0, 220, 250, 400, 0.01, mesh);
   //outputObj(mesh);
   GLWindows(gen.mesh);
   
   tsdf.freeData();
   */
}
int
 main (int argc, char** argv)
{
   runReconstruction();
   //_CrtDumpMemoryLeaks();
   return 0;
}
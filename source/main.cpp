
#include "RGBDReader.h"
#include "TsdfModel.h"
#include "HashTsdfModel.h"
#include "MultiTsdfModel.h"
#include "SNPGenerator.h"
#include "Segmentation.h"
#include "NeHeGL.h"



using namespace Rgbd;

void runReconstruction()
{
   RgbdReader reader;
   reader.Read("D:/PCL/record.ONI", 0, 10, 1, 0.1);
   reader.ShowPointCoud();


   MultiTsdfModel tsdf;
   int frameCount = reader.getFrameCount();
   for(int i = 0; i < frameCount; i++)
   {
	   reader.putDataToSdf(&tsdf, i);
   }

   /*
   double devide = 0.01;
   Eigen::Matrix4f zero;
   zero << 1, 0, 0, 0,  
			0, 1, 0, 0,  
			0, 0, 1, 0,  
			0, 0, 0, 1;
   TsdfModel tsdf(zero, 1.2, 1, 1.0, -1.5, -1.0, 0.0, devide, devide * 10, devide / 10);
   
   SNPGenerator gen(&tsdf, devide);
   
   int frameCount = reader.getFrameCount();
   Plane3D current_coff;
   for(int i = 0; i < frameCount; i++)
   {
	   reader.putDataToSdf(&tsdf, i);
	   gen.castNextMeshes(reader.getTrans(i));
   }
   
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
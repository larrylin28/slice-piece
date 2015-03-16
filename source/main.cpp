
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
   reader.Read("D:/PCL/record.ONI", 80, 10, 1, 0.1);
   reader.ShowPointCoud();

   /*
   Segmentation seg;
   
   Block b;
   b.local << 0.92240542, -0.056551334, -0.38206053, 0.00000000,
   -0.056551334, 0.95878512, -0.27844760, 0.00000000,
   0.38206053, 0.27844760, 0.88119048, 0.00000000,
   0.00000000, 0.00000000, 0.00000000, 1.0000000;
   b.coff.a = 0.38206051637059907;
   b.coff.b = 0.27844759728408325;
   b.coff.c = 0.88119050006077126;
   b.coff.d = 2.1890566277146868;
   b.bound.bound[0] = -2.3401529788970947;
   b.bound.bound[1] = -1.5933251380920410;
   b.bound.bound[2] = 1.9525583982467651;
   b.bound.bound[3] = 0.17319589853286743;
   b.bound.bound[4] = -0.35530942678451538;
   b.bound.bound[5] = 2.3175694942474365;
   seg.seged_blocks.push_back(b);
   */

   MultiTsdfModel tsdf;
   int frameCount = reader.getFrameCount();
   for(int i = 0; i < frameCount; i++)
   {
	   reader.putDataToSdf(&tsdf, i);
   }
   //HashTsdfModel hsdf;
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
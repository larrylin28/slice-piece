#include <pcl/features/integral_image_normal.h>
#include <pcl/surface/gp3.h>
#include "kernel.cuh"



#include "sdf.h"
#include "math.h"
#include "Vectors.h"
#include <iostream>
#include <fstream>




#define math_e (2.7182818284590452353602874713526624977572)
#define index(x,y,z) ((z*Ylen)+y)*Xlen+x
#define isnan(x) ((x) != (x)) 


DistandNorm DepthCalculator::getDist(double x,double y,double z)
{
	//逆变换获得原始坐标
	Eigen::Vector4f v;
	v[0] = x;
	v[1] = y;
	v[2] = z;
	v[3] = 1;
	Eigen::Matrix4f mc =  tran.inverse();
	Eigen::Vector4f t = mc * v;
	XnPoint3D real[1];
	XnPoint3D proj[1];
	real[0].X = t[0] / 0.001f;
	real[0].Y = t[1] / 0.001f;
	real[0].Z = t[2] / 0.001f;
	depthGenerator->ConvertRealWorldToProjective(1,real,proj);
	int px = (int)proj[0].X;
	int py = (int)proj[0].Y;
	if(px < 0 || px >= 640 || py < 0 || py >= 480){
		DistandNorm ret;
		ret.dist = 100;
		return ret;
	}
	pcl::PointXYZRGB pixel = cloud->points[py*cloud->width+px];
	pcl::Normal pixelN = normals->points[py*normals->width+px];
	if(isnan(pixelN.normal_x)){
		DistandNorm ret;
		ret.dist = 100;
		return ret;
	}
					
	//求点到点（像素）距离	
	XnPoint3D real2[1];
	XnPoint3D proj2[1];
	Eigen::Vector4f v2;
	v2[0] = pixel.x ;
	v2[1] = pixel.y;
	v2[2] = pixel.z;
	v2[3] = 1;
	Eigen::Vector4f t2 = mc * v2;
	real2[0].X = t2[0] / 0.001f;
	real2[0].Y = t2[1] / 0.001f;
	real2[0].Z = t2[2] / 0.001f;
	depthGenerator->ConvertRealWorldToProjective(1,real2,proj2);
	//double d_p = sqrt((x-pixel.x)*(x-pixel.x)+(y-pixel.y)*(y-pixel.y)+(z-pixel.z)*(z-pixel.z));
	double d_p = (proj[0].Z-proj2[0].Z)* 0.001f;
	//求点到平面（近似）距离
	double a[3];
	double b[3];
	a[0] = pixel.x-x;
	a[1] = pixel.y-y;
	a[2] = pixel.z-z;
	b[0] = pixelN.normal_x;
	b[1] = pixelN.normal_y;
	b[2] = pixelN.normal_z;
	double d_n = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    //深度：
	double depth = d_p;

	DistandNorm ret;
	ret.dist = depth;
	ret.n_x = pixelN.normal_x;
	ret.n_y = pixelN.normal_y;
	ret.n_z = pixelN.normal_z;

	return ret;
}

const double TsdfNode::minDiv = 0.005;

TsdfNode::TsdfNode()
{
}

void TsdfNode::Initial(double ix,double iy,double iz,double ax,double ay,double az,double dev)
{
    minX=ix;
    minY=iy;
    minZ=iz;
	maxX=ax;
    maxY=ay;
    maxZ=az;

    devide = dev;
    thre = 0.1;
    epsi = devide/10;

	Xlen = (ax-ix)/devide +0.5;
    Ylen = (ay-iy)/devide +0.5;
    Zlen = (az-iz)/devide +0.5;

	dist = 100;
	n_x = 0;
	n_y = 0;
	n_z = 0;
	weight = 0;

	child = NULL;
}

double TsdfNode::getWeight(double depth){
	double Sigma = 500;

	if(depth < epsi && depth > -epsi) return 1;
	if(depth > thre || depth < -thre) return 0;

	double diff;
	if(depth > 0) diff = depth - epsi;
	else diff = depth + epsi;
	return pow(math_e,-Sigma*diff*diff);
}

//遍历更新
void TsdfNode::Traversal(DepthCalculator* dg)
{
	if(child != NULL)
	{
	    for(int i = 0;i < Xlen;i++){
			double coreX = minX + ((double)i + 0.5)*devide;
			for(int j = 0;j < Ylen;j++){
				double coreY = minY + ((double)j + 0.5)*devide;
				for(int l = 0;l < Zlen;l++){
					double coreZ = minZ + ((double)l + 0.5)*devide;
					DistandNorm dnn = dg->getDist(coreX,coreY,coreZ);
				    UpdateChild(i,j,l,dnn);
					if(child[i][j][l].child != NULL)
					{
						child[i][j][l].Traversal(dg);
					}else if(abs(child[i][j][l].dist) < (devide*2) && devide > minDiv){
						child[i][j][l].split();
						child[i][j][l].Traversal(dg);
					}
			}
		}
	}
	}
}

void TsdfNode::TraversalDraw()
{
	if(child != NULL)
	{
	    for(int i = 0;i < Xlen;i++){
			double coreX = minX + ((double)i + 0.5)*devide;
			for(int j = 0;j < Ylen;j++){
				double coreY = minY + ((double)j + 0.5)*devide;
				for(int l = 0;l < Zlen;l++){
					double coreZ = minZ + ((double)l + 0.5)*devide;
					if(child[i][j][l].child != NULL)
					{
						child[i][j][l].TraversalDraw();
					}else{
						if(abs(child[i][j][l].dist) < minDiv)
						{
							glNormal3f(child[i][j][l].n_x, child[i][j][l].n_y, child[i][j][l].n_z);
							glVertex3f(coreX,coreY, coreZ);
						}
					}
			}
		}
	}
	}
}


void TsdfNode::split()
{
	child = new TsdfNode**[Xlen];
	for(int i = 0;i < Xlen;i++){
		child[i] = new TsdfNode*[Ylen];
		for(int j = 0;j < Ylen;j++){
			child[i][j] = new TsdfNode[Zlen];
			for(int l = 0;l < Zlen;l++){
				child[i][j][l].Initial(minX+i*devide,minY+j*devide,minZ+l*devide,minX+i*devide+devide,minY+j*devide+devide,minZ+l*devide+devide,devide/10);
			}
		}
	}
}

void TsdfNode::UpdateChild(int x,int y,int z,DistandNorm dnn)
{
	if(dnn.dist >= 100) return;
	double wei = getWeight(dnn.dist);
	if(wei > 0){
		child[x][y][z].dist = (child[x][y][z].dist*child[x][y][z].weight+dnn.dist*wei)/(child[x][y][z].weight+wei);
		double nx = (child[x][y][z].n_x*child[x][y][z].weight+dnn.n_x*wei)/(child[x][y][z].weight+wei);
		double ny = (child[x][y][z].n_y*child[x][y][z].weight+dnn.n_y*wei)/(child[x][y][z].weight+wei);
		double nz = (child[x][y][z].n_z*child[x][y][z].weight+dnn.n_z*wei)/(child[x][y][z].weight+wei);
		double n = sqrt(nx*nx+ny*ny+nz*nz);
		child[x][y][z].n_x = nx/n;
		child[x][y][z].n_y = ny/n;
		child[x][y][z].n_z = nz/n;
		child[x][y][z].weight = child[x][y][z].weight+wei;
	}
}

Sdf::Sdf():oldp(NULL)
{
	local_tran<< 1, 0, 0, 0,  
				0, 1, 0, 0,  
				0, 0, 1, 0,  
				0, 0, 0, 1;

	maxX=1.2;
    maxY=1;
    maxZ=4.0;
    minX=-3.5;
    minY=-2.0;
    minZ=-1.0;

    devide = 0.01;
    thre = 0.1;
    epsi = 0.001;

	mu = 2;
	surfaceAng = (M_PI/4);  
    minAng = (M_PI/18);
    maxAng = (2*M_PI/3);

	Xlen = 470;
    Ylen = 300;
    Zlen = 500;

	int size = Xlen*Ylen*Zlen;
    dist = new float[size];
	n_x = new float[size];
	n_y = new float[size];
	n_z = new float[size];
    weight = new float[size];
	for(int i = 0;i < size;i++){
		dist[i] = 100;
		n_x[i] = 0;
		n_y[i] = 0;
		n_z[i] = 0;
		weight[i] = 0;
	}

	gen = NULL;

	float* tran = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			tran[i*4+j] = local_tran(i,j);
		}
	}
	int ret = initialSdf(minX,minY,minZ,Xlen,Ylen,Zlen,devide,epsi,thre,dist,weight,n_x,n_y,n_z,tran);
	
	delete dist;
	delete weight;
	delete n_x;
	delete n_y;
	delete n_z;
}

Sdf::Sdf(Eigen::Matrix4f ltran, double maxx,double maxy, double maxz,double minx,double miny,double minz,double dev,double th,double ep)
	:local_tran(ltran),maxX(maxx),maxY(maxy),maxZ(maxz),minX(minx),minY(miny),minZ(minz),devide(dev),thre(th),epsi(ep),oldp(NULL)
{
	mu = 2;
	surfaceAng = (M_PI/4);  
    minAng = (M_PI/18);
    maxAng = (2*M_PI/3);

	Xlen = (int)((maxx-minx)/devide) + 1;
	Ylen = (int)((maxy-miny)/devide) + 1;
	Zlen = (int)((maxz-minz)/devide) + 1;

	int size = Xlen*Ylen*Zlen;
    dist = new float[size];
	n_x = new float[size];
	n_y = new float[size];
	n_z = new float[size];
    weight = new float[size];
	for(int i = 0;i < size;i++){
		dist[i] = 100;
		n_x[i] = 0;
		n_y[i] = 0;
		n_z[i] = 0;
		weight[i] = 0;
	}

	gen = NULL;

	float* tran = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			tran[i*4+j] = local_tran(i,j);
		}
	}
	int ret = initialSdf(minX,minY,minZ,Xlen,Ylen,Zlen,devide,epsi,thre,dist,weight,n_x,n_y,n_z,tran);

	delete dist;
	delete weight;
	delete n_x;
	delete n_y;
	delete n_z;
}

void Sdf::changeLocal(Eigen::Matrix4f ltran, double maxx,double maxy, double maxz,double minx,double miny,double minz)
{
	maxX = maxx;
    maxY = maxy;
    maxZ = maxz;
    minX = minx;
    minY = miny;
    minZ = minz;

	Xlen = (int)((maxx-minx)/devide) + 1;
	Ylen = (int)((maxy-miny)/devide) + 1;
	Zlen = (int)((maxz-minz)/devide) + 1;

	float* tran = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			tran[i*4+j] = ltran(i,j);
		}
	}

	float* omc = new float[4*4];
	Eigen::Matrix4f olocalc = local_tran.inverse();
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			omc[i*4+j] = olocalc(i,j);
		}
	}

	int ret = changeSdf(minX, minY, minZ, Xlen, Ylen, Zlen, tran, omc);

	local_tran = ltran;
}

double Sdf::getWeight(double depth){
	double Sigma = 500;

	if(depth < epsi && depth > -epsi) return 1;
	if(depth > thre || depth < -thre) return 0;

	double diff;
	if(depth > 0) diff = depth - epsi;
	else diff = depth + epsi;
	return pow(math_e,-Sigma*diff*diff);
}

void Sdf::Traversal(DepthCalculator* dg){
	for(int i = 0;i < Sdf::sdf->Xlen;i++){
		double coreX = Sdf::sdf->minX + ((double)i + 0.5)*Sdf::sdf->devide;
		for(int j = 0;j < Sdf::sdf->Ylen;j++){
			double coreY = Sdf::sdf->minY + ((double)j + 0.5)*Sdf::sdf->devide;
			for(int l = 0;l < Sdf::sdf->Zlen;l++){
				double coreZ = Sdf::sdf->minZ + ((double)l + 0.5)*Sdf::sdf->devide;
				DistandNorm dnn = dg->getDist(coreX,coreY,coreZ);
				double weight =  Sdf::sdf->getWeight(dnn.dist);
				if(weight > 0){
					Update(i,j,l,dnn,weight);
				}
			}
		}
	}
}

void Sdf::Traversal2(DepthCalculator* dg)
{
	int size = dg->cloud->points.size();
	float* ptsx = new float[size];
	float* ptsy = new float[size];
	float* ptsz = new float[size];
	float* norx = new float[size];
	float* nory = new float[size];
	float* norz = new float[size];
	for(int i = 0;i < size;i++)
	{
		ptsx[i] = dg->cloud->points[i].x;
		ptsy[i] = dg->cloud->points[i].y;
		ptsz[i] = dg->cloud->points[i].z;
		norx[i] = dg->normals->points[i].normal_x;
		nory[i] = dg->normals->points[i].normal_y;
		norz[i] = dg->normals->points[i].normal_z;
	}
	Eigen::Matrix4f mcT = dg->tran.inverse();
	float* mc = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			mc[i*4+j] = mcT(i,j);
		}
	}
	int width = dg->cloud->width;
	XnFieldOfView fov;
	dg->depthGenerator->GetFieldOfView(fov);
	XnMapOutputMode outputMode;
	dg->depthGenerator->GetMapOutputMode(outputMode);

	double fXToZ = tan(fov.fHFOV/2)*2;	
	double fYToZ = tan(fov.fVFOV/2)*2;

	XnUInt32 nHalfXres = outputMode.nXRes / 2;	
	XnUInt32 nHalfYres = outputMode.nYRes / 2;

	XnDouble fCoeffX = outputMode.nXRes / fXToZ;
	XnDouble fCoeffY = outputMode.nYRes / fYToZ;

	int ret = updateSdf(ptsx, ptsy, ptsz, norx, nory, norz,size, mc,width,nHalfXres,nHalfYres,fCoeffX,fCoeffY);
	delete[] ptsx;
	delete[] ptsy;
    delete[] ptsz;
	delete[] norx;
	delete[] nory;
	delete[] norz;
	delete[] mc;
}

void Sdf::getData()
{
	int size = Xlen*Ylen*Zlen;
	int ret = getSdf(dist,weight,n_x,n_y,n_z,size);
}

void Sdf::freeData()
{
	int ret = freeSdf();
}

void Sdf::TraversalDraw(){
	for(int i = 0;i < Sdf::sdf->Xlen;i++){
		double coreX = Sdf::sdf->minX + ((double)i + 0.5)*Sdf::sdf->devide;
		for(int j = 0;j < Sdf::sdf->Ylen;j++){
			double coreY = Sdf::sdf->minY + ((double)j + 0.5)*Sdf::sdf->devide;
			for(int l = 0;l < Sdf::sdf->Zlen;l++){
				double coreZ = Sdf::sdf->minZ + ((double)l + 0.5)*Sdf::sdf->devide;
				int id = ((l*Ylen)+j)*Xlen+i;
				if(abs(dist[id]) < devide)
				{
					glNormal3f(n_x[id], n_y[id], n_z[id]);
					glVertex3f(coreX,coreY, coreZ);
				}
			}
		}
	}
}

void Sdf::Update(int x,int y,int z,DistandNorm dnn,double wei){
	int id = ((z*Ylen)+y)*Xlen+x;
	dist[id] = (dist[id]*weight[id]+dnn.dist*wei)/(weight[id]+wei);
	double nx = (n_x[id]*weight[id]+dnn.n_x*wei)/(weight[id]+wei);
	double ny = (n_y[id]*weight[id]+dnn.n_y*wei)/(weight[id]+wei);
	double nz = (n_z[id]*weight[id]+dnn.n_z*wei)/(weight[id]+wei);
	double n = sqrt(nx*nx+ny*ny+nz*nz);
	n_x[id] = nx/n;
	n_y[id] = ny/n;
	n_z[id] = nz/n;
	weight[id] = weight[id]+wei;
}

void Sdf::GenerateMeshes(){
	
	//MC算法
	gen = new CIsoSurface<double>();
	long size = Xlen*Ylen*Zlen;
	double* absdist = new double[size];
	long index = 0;
	for(int l = 0;l < Zlen;l++){
		for(int j = 0;j < Ylen;j++){
			for(int i = 0;i < Xlen;i++){
				int id = ((l*Ylen)+j)*Xlen+i;
				if(weight[id] > 0 ){
					absdist[index] = 0.02/abs(dist[id]);
				}else{
					absdist[index] = 0.02/100;
				}
				index++;
			}
		}
	}
	gen->GenerateSurface(absdist,1,Xlen-1,Ylen-1,Zlen-1,devide,devide,devide);
}

void Sdf::cloudToMeshes(Eigen::Matrix4f tran)
{
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud = castCloud(tran);
	pcl::PointCloud<pcl::Normal>::Ptr normals (new pcl::PointCloud<pcl::Normal>);
	//computeNormal(cloud,normals);
	pcl::IntegralImageNormalEstimation<pcl::PointXYZRGB, pcl::Normal> ne;
	ne.setNormalEstimationMethod (ne.AVERAGE_3D_GRADIENT);
	ne.setMaxDepthChangeFactor(0.03f);
	ne.setNormalSmoothingSize(20.0f);
	ne.setInputCloud(cloud);
	ne.compute(*normals);

	// Concatenate the XYZ and normal fields*
    pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals (new pcl::PointCloud<pcl::PointNormal>);
	pcl::copyPointCloud(*cloud,*cloud_with_normals);
	pcl::copyPointCloud(*normals,*cloud_with_normals);
    //* cloud_with_normals = cloud + normals

    // Create search tree*
    pcl::search::KdTree<pcl::PointNormal>::Ptr tree2 (new pcl::search::KdTree<pcl::PointNormal>);
    tree2->setInputCloud (cloud_with_normals);

	// Initialize objects
	 pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
	 pcl::PolygonMesh triangles;

	 // Set the maximum distance between connected points (maximum edge length)
	gp3.setSearchRadius (0.025);

	// Set typical values for the parameters
	 gp3.setMu (2.5);
	 gp3.setMaximumNearestNeighbors (100);
	 gp3.setMaximumSurfaceAngle(M_PI/4); // 45 degrees
	 gp3.setMinimumAngle(M_PI/18); // 10 degrees
	 gp3.setMaximumAngle(2*M_PI/3); // 120 degrees
	 gp3.setNormalConsistency(false);

	// Get result
	gp3.setInputCloud (cloud_with_normals);
	gp3.setSearchMethod (tree2);
	gp3.reconstruct (triangles);

	 // Additional vertex information
	std::vector<int> parts = gp3.getPartIDs();
	std::vector<int> states = gp3.getPointStates();
	
}
pcl::PointCloud<pcl::PointXYZRGB>::Ptr Sdf::castCloud(Eigen::Matrix4f tran){
	double ix = 10000,iy = 10000,ax = -10000,ay = -10000;

	Eigen::Vector4f v;
	Eigen::Vector4f t;
	Eigen::Matrix4f mc =  tran.inverse();

	//求x,y值区间
	v[0] = minX;
	v[1] = minY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = minY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = maxY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = maxY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = minY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = minY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = maxY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = maxY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	//对区间内进行采样
	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	pcl::PointCloud<pcl::PointXYZRGB>::Ptr ret(new pcl::PointCloud<pcl::PointXYZRGB>);
    ret->height = yl;
	ret->width = xl;
	ret->is_dense = false;

	ret->points.resize(ret->height * ret->width);

	for(int i = 0;i < xl;i++)
	{
		for(int j = 0;j < yl;j++)
		{
			DistandNorm dnn  = getDepth(ix + i*devide,iy + j*devide,tran);
			double cx = (ix + i*devide);
			double cy = (iy + j*devide);
			double cz = dnn.dist;
			if(dnn.dist <= -10000)
			{
				ret->points[j*xl+i].x = 0;
				ret->points[j*xl+i].y = 0;
				ret->points[j*xl+i].z = -10000;
			}else{
				ret->points[j*xl+i].x = tran(0,0)*cx + tran(0,1)*cy + tran(0,2)*cz + tran(0,3);
				ret->points[j*xl+i].y = tran(1,0)*cx + tran(1,1)*cy + tran(1,2)*cz + tran(1,3);
				ret->points[j*xl+i].z = tran(2,0)*cx + tran(2,1)*cy + tran(2,2)*cz + tran(2,3);
			}
			
		}
	}
	return ret;
}

void Sdf::TryAddRightMesh(int a,int b,int c,Eigen::Matrix4f tran)
{
	//if(GL_Mesh::gl_vertexes->at(a).arc < 6 && GL_Mesh::gl_vertexes->at(b).arc < 6 && GL_Mesh::gl_vertexes->at(c).arc < 6
    if( GL_Mesh::findEdgeDegree(a,b) < 2 && GL_Mesh::findEdgeDegree(b,c) < 2 && GL_Mesh::findEdgeDegree(c,a) < 2)
	{
		TryaddMesh(a,b,c,tran);
	}
}

bool onTriangle(int v[3],double x,double y,double z,double th,Eigen::Matrix4f lc)
{
	Eigen::Vector4f t[3];
	for(int i = 0;i < 3;i++)
	{
		GLVertex vv = GL_Mesh::gl_vertexes->at(v[i]);
		//将上一帧的点投影到本帧投影平面上
		Eigen::Vector4f vec;
		vec[0] = vv.v[0];
		vec[1] = vv.v[1];
		vec[2] = vv.v[2];
		vec[3] = 1;
		t[i] = lc*vec; //投影到当前帧坐标
	}
	double s1 = (t[0][0] - x)*(t[0][1]-t[1][1])-(t[0][0]-t[1][0])*(t[0][1] - y);
	double s2 = (t[1][0] - x)*(t[1][1]-t[2][1])-(t[1][0]-t[2][0])*(t[1][1] - y);
	double s3 = (t[2][0] - x)*(t[2][1]-t[0][1])-(t[2][0]-t[0][0])*(t[2][1] - y);
	if((s1*s2 > 0)&&(s2*s3 > 0))
	{
		double zp1 = t[0][2] + (y-t[0][1])*(t[1][2]-t[0][2])/(t[1][1]-t[0][1]);
		double zp2 = t[0][2] + (y-t[0][1])*(t[2][2]-t[0][2])/(t[2][1]-t[0][1]);
		double xp1 = t[0][0] + (y-t[0][1])*(t[1][0]-t[0][0])/(t[1][1]-t[0][1]);
		double xp2 = t[0][0] + (y-t[0][1])*(t[2][0]-t[0][0])/(t[2][1]-t[0][1]);

		double zp = zp1 + (x-xp1)*(zp2-zp1)/(xp2-xp1);
		if(abs(z - zp) < th) return true;
		else return false;
	}
	else if(s1 == 0 && x >= min(t[0][0],t[1][0]) && x <= max(t[0][0],t[1][0]))
	{
		double zp = t[0][2] + (x-t[0][0])*(t[1][2]-t[0][2])/(t[1][0]-t[0][0]);
		if(abs(z - zp) < th) return true;
		else return false;
	}
	else if(s2 == 0 && x >= min(t[1][0],t[2][0]) && x <= max(t[1][0],t[2][0]))
	{
		double zp = t[1][2] + (x-t[1][0])*(t[2][2]-t[1][2])/(t[2][0]-t[1][0]);
		if(abs(z - zp) < th) return true;
		else return false;
	}
	else if(s3 == 0 && x >= min(t[2][0],t[0][0]) && x <= max(t[2][0],t[0][0]))
	{
		double zp = t[2][2] + (x-t[2][0])*(t[0][2]-t[2][2])/(t[0][0]-t[2][0]);
		if(abs(z - zp) < th) return true;
		else return false;
	}
	return false;
}

int** Sdf::castPieceMeshes(int tranID,int** pre,int lastID)
{
	Eigen::Matrix4f& tran = GL_Mesh::gl_trans->at(tranID);
	Eigen::Matrix4f& last = GL_Mesh::gl_trans->at(lastID);
	Eigen::Matrix4f mc =  tran.inverse();
	Eigen::Matrix4f lc =  last.inverse();
	Eigen::Matrix4f localc = local_tran.inverse();
	Eigen::Matrix4f t_in_local = localc*tran;

	//单帧镜头的参数
	double ix = -1.5-devide*0.5;
	double iy = -1-devide*0.5;
	double ax = 1.2;
	double ay = 1;

	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	int ** planePoints = NULL;
	if(xl > 0 ) planePoints = new int*[xl];
	for(int i = 0;i < xl;i++)
	{
		planePoints[i] = new int[yl];
	}

	int* ixid = new int[yl];
	int* axid = new int[yl];
	for(int i = 0;i < yl;i++)
	{
		ixid[i] = xl;
		axid[i] = -1;
	}
	
	if(GL_Mesh::gl_meshes->size() > 0)
	{
		std::map<std::tuple<int,int,int>,GLMesh>::iterator it = GL_Mesh::gl_meshes->begin();
		for(;it != GL_Mesh::gl_meshes->end();it++)
		{
			int v[3];
			v[0] = it->second.vertexes[0];
			v[1] = it->second.vertexes[1];
			v[2] = it->second.vertexes[2];
			coverPointinTriangle(v,planePoints,xl,yl,mc,ix,iy);
		}
	}
	
	/*
	if(pre != NULL)
	{
		
		for(int i = 0;i < xl;i++)
		{
			for(int j = 0;j < yl;j++)
			{
				if(i > 0 && j > 0 && pre[i][j] >= 0)
				{
					 coverPointinMeshes(pre,planePoints,i,j,xl,yl,mc,ix,iy);
				}
			}
		}
		
		
		for(int i = 0;i < xl;i++)
		{
			for(int j = 0;j < yl;j++)
			{
				if(pre[i][j] >= 0)
				{
					GLVertex vv = GL_Mesh::gl_vertexes->at(pre[i][j]);
					//将上一帧的点投影到本帧投影平面上
					Eigen::Vector4f v;
					Eigen::Vector4f t;
					v[0] = vv.v[0];
					v[1] = vv.v[1];
					v[2] = vv.v[2];
					v[3] = 1;
					t = mc*v; //投影到当前帧坐标
					int xid = static_cast<int>((t[0]-ix)/devide+0.5);
					int yid = static_cast<int>((t[1]-iy)/devide+0.5);
					if(xid >= 0 && xid < xl && yid >=0 && yid <yl)
					{
						planePoints[xid][yid] = pre[i][j];
						if(xid < ixid[yid]) ixid[yid] = xid;
						if(xid > axid[yid]) axid[yid] = xid;
					}
				}
			}
		}
		
		for(int i = 0;i < yl;i++)
		{
			if(ixid[i] < xl)
			{
				int last = planePoints[ixid[i]][i];
				for(int j = ixid[i]+1;j < axid[i];j++)
				{
					//if(planePoints[j][i] < 0 && (planePoints[j+1][i] > 0 || (j+2 <= axid[i] && planePoints[j+2][i] > 0))) planePoints[j][i] = last;
					if(planePoints[j][i] < 0) planePoints[j][i] = last;
					else last = planePoints[j][i];
				}
			}
		}
		
	}
   */
	std::set<int> bounds;
	int forbid = GL_Mesh::gl_vertexes->size();

	float* pmc = new float[4*4];
	for(int i = 0;i < 4;i++)
	{
		for(int j = 0;j < 4;j++)
		{
			pmc[i*4+j] = t_in_local(i,j);
		}
	}
	float* p_dists = new float[xl*yl];
	float* p_nx = new float[xl*yl];
	float* p_ny = new float[xl*yl];
	float* p_nz = new float[xl*yl];
	getAllDepths(ix, iy, devide, xl, yl, pmc, p_dists, p_nx,  p_ny,  p_nz);
	for(int i = 0;i < xl;i++)
	{
		for(int j = 0;j < yl;j++)
		{
			/*
			bool newp = false;
			DistandNorm dnn = getDepth(ix + i*devide,iy + j*devide,tran);
			if(dnn.dist > -10000 && !isnan(dnn.n_x)){
				GLVertex v1 = *(new GLVertex());
				putPoint(v1,ix+i*devide,iy+j*devide,dnn,tran);
				v1.tran = tranID;
				
				if(planePoints[i][j] >= 0)
				{
					GLVertex nv = GL_Mesh::gl_vertexes->at(planePoints[i][j]);
					double dis = sqrt((nv.v[0]-v1.v[0])*(nv.v[0]-v1.v[0]) + (nv.v[1]-v1.v[1])*(nv.v[1]-v1.v[1]) + (nv.v[2]-v1.v[2])*(nv.v[2]-v1.v[2]));
					GL_Mesh::gl_vertexes->push_back(v1);
					planePoints[i][j] =  GL_Mesh::gl_vertexes->size()-1;
					newp = true;
				}else{
					GL_Mesh::gl_vertexes->push_back(v1);
					planePoints[i][j] =  GL_Mesh::gl_vertexes->size()-1;
					newp = true;
				}
					
			}else{
				//nowline[j] = -1;
				planePoints[i][j] = -1;
			}
			*/
			
			if(planePoints[i][j] < 0)
			{

				int pid = (j*xl)+i;
				DistandNorm dnn;
				dnn.dist = p_dists[pid];
				dnn.n_x = p_nx[pid];
				dnn.n_y = p_ny[pid];
				dnn.n_z = p_nz[pid];

				//DistandNorm dnn = getDepth(ix + i*devide,iy + j*devide,t_in_local);
				//DistandNorm dnn = getDepth(ix + i*devide,iy + j*devide,tran);
				if(dnn.dist > -10000 && !isnan(dnn.n_x)){
					GLVertex v1 = *(new GLVertex());
					putPoint(v1,ix+i*devide,iy+j*devide,dnn,tran);
					v1.tran = tranID;
					
						GL_Mesh::gl_vertexes->push_back(v1);
						//GL_Mesh::o_meshes.add_vertex(MyMesh::Point(v1.v[0],v1.v[1],v1.v[2]));
						//nowline[j] = GL_Mesh::gl_vertexes->size()-1;
						planePoints[i][j] =  GL_Mesh::gl_vertexes->size()-1;
				}else{
					//nowline[j] = -1;
					planePoints[i][j] = -1;
				}
			}
			
			if(i > 0 && j > 0 && planePoints[i][j] >= 0)
			{
				bool tic = buildMeshesAtR(planePoints,i,j,tran,forbid);
				
				if(tic)
				{
					//if(planePoints[i][j] >= forbid) bounds.insert(planePoints[i][j]);
					if(planePoints[i][j] >=  forbid) bounds.insert(planePoints[i][j]);
					if(planePoints[i-1][j] >=  forbid) bounds.insert(planePoints[i-1][j]);
					if(planePoints[i][j-1] >=  forbid) bounds.insert(planePoints[i][j-1]);
					if(planePoints[i-1][j-1] >= forbid) bounds.insert(planePoints[i-1][j-1]);
				}
				
				//若存在缝合线
				/*
				if(tic)
				{
					std::set<int> nearest;
					GLVertex R = GL_Mesh::gl_vertexes->at(planePoints[i][j]);
					//本投影平面中的近邻
					for(int jx = i-mu;jx <= i+mu;jx++)
					{
						if(jx < 0) continue;
						if(jx >= xl) break;
						for(int jy = j-mu;jy <= j+mu;jy++)
						{
							if(jy < 0) continue;
							if(jx == i && jy == j) continue;
							if(jy >= yl) break;
							if(planePoints[jx][jy] >= 0)
							{
								GLVertex nv = GL_Mesh::gl_vertexes->at(planePoints[jx][jy]);
								double dis = sqrt((nv.v[0]-R.v[0])*(nv.v[0]-R.v[0]) + (nv.v[1]-R.v[1])*(nv.v[1]-R.v[1]) + (nv.v[2]-R.v[2])*(nv.v[2]-R.v[2]));
								if(dis < 2*mu*devide) nearest.insert(planePoints[jx][jy]);
							}
						}
					}
					if(pre != NULL){
						//上个投影平面中的近邻
						Eigen::Vector4f v;
						Eigen::Vector4f t;
						v[0] = R.v[0];
						v[1] = R.v[1];
						v[2] = R.v[2];
						v[3] = 1;
						t = lc*v; //投影到上一帧坐标
						int xid = static_cast<int>((t[0]-ix)/devide+0.5);
						int yid = static_cast<int>((t[1]-iy)/devide+0.5);
						for(int jx = xid - mu;jx <= xid + mu;jx++)
						{
							if(jx < 0) continue;
							if(jx >= xl) break;
							for(int jy = yid-mu;jy <= yid+mu;jy++)
							{
								if(jy < 0) continue;
								if(jy >= yl) break;
								if(pre[jx][jy] >= 0)
								{
									GLVertex nv = GL_Mesh::gl_vertexes->at(pre[jx][jy]);
									double dis = sqrt((nv.v[0]-R.v[0])*(nv.v[0]-R.v[0]) + (nv.v[1]-R.v[1])*(nv.v[1]-R.v[1]) + (nv.v[2]-R.v[2])*(nv.v[2]-R.v[2]));
									if(dis < 2*mu*devide) nearest.insert(pre[jx][jy]);
								}
							}
						}
					}
					nearest.erase(planePoints[i][j]);
					buildMeshesGreedy(nearest,planePoints[i][j],tran);

				}
				*/
				
			} 
		}
	}
	
	GL_Mesh::updateVertexes();
	
	int m = 1;
	if(pre)
	{
		for(int i = 1;i < xl;i++)
		{
			if(pre[i][0] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i+1 < xl && pre[i+1][0] >= 0) except[0] = pre[i+1][0];
				if(i-1 >= 0 && pre[i-1][0] >= 0) except[1] = pre[i-1][0];
				buildTic(i,0,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = 1;i < yl;i++)
		{
			if(pre[xl-1][i] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i+1 < yl && pre[xl-1][i+1]  >= 0) except[0] = pre[xl-1][i+1];
				if(i-1 >= 0 && pre[xl-1][i-1]  >= 0) except[1] = pre[xl-1][i-1];
				buildTic(xl-1,i,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = xl-2;i >= 0;i--)
		{
			if(pre[i][yl-1] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i-1 >= 0 && pre[i-1][yl-1]  >= 0) except[0] = pre[i-1][yl-1];
				if(i+1 < xl && pre[i+1][yl-1]  >= 0) except[1] = pre[i+1][yl-1];
				buildTic(i,yl-1,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = yl-2;i >= 0;i--)
		{
			if(pre[0][i] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i-1 >= 0 && pre[0][i-1]  >= 0) except[0] = pre[0][i-1];
				if(i+1 < yl && pre[0][i+1]  >= 0) except[1] = pre[0][i+1];
				buildTic(0,i,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
	}
	m = 2;
	if(pre)
	{
		for(int i = 1;i < xl;i++)
		{
			if(pre[i][0] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i+1 < xl && pre[i+1][0] >= 0) except[0] = pre[i+1][0];
				if(i-1 >= 0 && pre[i-1][0] >= 0) except[1] = pre[i-1][0];
				buildTic(i,0,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = 1;i < yl;i++)
		{
			if(pre[xl-1][i] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i+1 < yl && pre[xl-1][i+1]  >= 0) except[0] = pre[xl-1][i+1];
				if(i-1 >= 0 && pre[xl-1][i-1]  >= 0) except[1] = pre[xl-1][i-1];
				buildTic(xl-1,i,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = xl-2;i >= 0;i--)
		{
			if(pre[i][yl-1] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i-1 >= 0 && pre[i-1][yl-1]  >= 0) except[0] = pre[i-1][yl-1];
				if(i+1 < xl && pre[i+1][yl-1]  >= 0) except[1] = pre[i+1][yl-1];
				buildTic(i,yl-1,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
		for(int i = yl-2;i >= 0;i--)
		{
			if(pre[0][i] >= 0)
			{
				int except[2];
				except[0]= -1;
				except[1]= -1;
				if(i-1 >= 0 && pre[0][i-1]  >= 0) except[0] = pre[0][i-1];
				if(i+1 < yl && pre[0][i+1]  >= 0) except[1] = pre[0][i+1];
				buildTic(0,i,ix,iy,xl,yl,pre,planePoints,mc,tran, forbid,except,m);
			}
		}
	}
	
	
	for(std::set<int>::iterator it = bounds.begin();it != bounds.end();it++)
	{
		int b = *it;
		std::set<int> nearest;
	    GLVertex& R = GL_Mesh::gl_vertexes->at(b);
		Eigen::Vector4f v;
		v[0] = R.v[0];
		v[1] = R.v[1];
		v[2] = R.v[2];
		v[3] = 1;
		Eigen::Vector4f t1 = mc*v; 
	    int xid1 = static_cast<int>((t1[0]-ix)/devide+0.5);
	    int yid1 = static_cast<int>((t1[1]-iy)/devide+0.5);
	    findNearest(planePoints,xid1,yid1,xl,yl,b,nearest,2,forbid);
		Eigen::Vector4f t2 = lc*v; 
	    int xid2 = static_cast<int>((t2[0]-ix)/devide+0.5);
	    int yid2 = static_cast<int>((t2[1]-iy)/devide+0.5);
	    findNearest(pre,xid2,yid2,xl,yl,b,nearest,2,0);
		nearest.erase(b);
		int except[2];
		except[0]= -1;
		except[1]= -1;
		buildMeshesGreedy(nearest,b,tran,0,except);
	}
	
	return planePoints;
}

void Sdf::buildTic(int px,int py,double ix,double iy,double xl,double yl,int** pre,int** planePoints,Eigen::Matrix4f mc, Eigen::Matrix4f tran,int forbid,int* except,int m)
{
	std::set<int> nearest;
	GLVertex& R = GL_Mesh::gl_vertexes->at(pre[px][py]);
	Eigen::Vector4f v;
	Eigen::Vector4f t;
	v[0] = R.v[0];
	v[1] = R.v[1];
	v[2] = R.v[2];
	v[3] = 1;
	t = mc*v; //投影到当前帧坐标
	int xid = static_cast<int>((t[0]-ix)/devide+0.5);
	int yid = static_cast<int>((t[1]-iy)/devide+0.5);
	findNearest(planePoints,xid,yid,xl,yl,pre[px][py],nearest,m,forbid);
	if(nearest.size() <= 0) return;
	findNearest(pre,px,py,xl,yl,pre[px][py],nearest,m,0);
	nearest.erase(pre[px][py]);
	buildMeshesGreedy(nearest,pre[px][py],tran,forbid,except);
}

void Sdf::findNearestAd(int** planePoints,int xid,int yid,int xl,int yl,int Rid,std::set<int>& nearest,int m,int forbid)
{
}

void Sdf::findNearest(int** planePoints,int xid,int yid,int xl,int yl,int Rid,std::set<int>& nearest,int m,int forbid)
{
	GLVertex R = GL_Mesh::gl_vertexes->at(Rid);
	for(int jx = xid-m;jx <= xid+m;jx++)
	{
		if(jx < 0) continue;
		if(jx >= xl) continue;
		for(int jy = yid-m;jy <= yid+m;jy++)
		{
			if(jy < 0) continue;
			if(jy >= yl) continue;
			if(planePoints[jx][jy] >= forbid)
			{
				GLVertex nv = GL_Mesh::gl_vertexes->at(planePoints[jx][jy]);
				double dis = sqrt((nv.v[0]-R.v[0])*(nv.v[0]-R.v[0]) + (nv.v[1]-R.v[1])*(nv.v[1]-R.v[1]) + (nv.v[2]-R.v[2])*(nv.v[2]-R.v[2]));
				if(dis < 6*m*devide) nearest.insert(planePoints[jx][jy]);
			}
		}
	}
}

void Sdf::TryaddMesh(int v1,int v2,int v3,Eigen::Matrix4f tran)
{
	GLMesh m1 = *(new GLMesh(v1,v2,v3));	   
	if(validMesh(m1,tran)){
		int a = v1;
		int b = v2;
		int c = v3;
		if (a>b) std::swap(a, b);
		if (b>c) 
		{ 
			std::swap(b, c);
			if (a > b) std::swap(a, b);
		}
		if(GL_Mesh::gl_meshes->find(std::tuple<int,int,int>(a,b,c)) != GL_Mesh::gl_meshes->end()) return;
		std::pair<std::tuple<int,int,int>,GLMesh> ins(std::tuple<int,int,int>(a,b,c),m1);		
		GL_Mesh::gl_meshes->insert(ins);
		GL_Mesh::addEdge(v1,v2);
		GL_Mesh::addEdge(v2,v3);
		GL_Mesh::addEdge(v3,v1);
		GLVertex& vv1 = GL_Mesh::gl_vertexes->at(v1);
		GLVertex& vv2 = GL_Mesh::gl_vertexes->at(v2);
		GLVertex& vv3 = GL_Mesh::gl_vertexes->at(v3);
	    vv1.addDegree(vv2.v,vv3.v);
		vv2.addDegree(vv1.v,vv3.v);
		vv3.addDegree(vv1.v,vv2.v);
		GL_Mesh::gl_vertexes->at(v1).arc++;
		GL_Mesh::gl_vertexes->at(v2).arc++;
		GL_Mesh::gl_vertexes->at(v3).arc++;

		//作为一阶邻域，各个顶点互相进行Laplacian滤波
		/*
		vv1.LaplacianUpdate(vv2.v);
		vv1.LaplacianUpdate(vv3.v);

		vv2.LaplacianUpdate(vv1.v);
		vv2.LaplacianUpdate(vv3.v);

		vv3.LaplacianUpdate(vv1.v);
		vv3.LaplacianUpdate(vv2.v);
		*/
		//GL_Mesh::o_meshes.add_face(GL_Mesh::o_meshes.vertex_handle(v1),GL_Mesh::o_meshes.vertex_handle(v2),GL_Mesh::o_meshes.vertex_handle(v3));
	}
			
}

void Sdf::coverPointinMeshes(int** pre,int** pts,int i, int j, int xl,int yl,Eigen::Matrix4f mc,double minx,double miny)
{
	
	if(pre[i-1][j-1] >= 0 && pre[i-1][j] >= 0)
	{
		if(GL_Mesh::findMesh(pre[i][j],pre[i-1][j-1],pre[i-1][j]))
		{
			int v[3];
			v[0] = pre[i][j];
			v[1] = pre[i-1][j-1];
			v[2] = pre[i-1][j];
			coverPointinTriangle(v,pts,xl,yl,mc,minx,miny);
		}
	}
	
	if(pre[i-1][j-1] >= 0 && pre[i][j-1] >= 0)
	{
		if(GL_Mesh::findMesh(pre[i][j],pre[i-1][j-1],pre[i][j-1]))
		{
			int v[3];
			v[0] = pre[i][j];
			v[1] = pre[i-1][j-1];
			v[2] = pre[i][j-1];
			coverPointinTriangle(v,pts,xl,yl,mc,minx,miny);
		}
	}
	
}

bool inTriangle(Eigen::Vector4f t[3],double x,double y)
{
	double s1 = (t[0][0] - x)*(t[0][1]-t[1][1])-(t[0][0]-t[1][0])*(t[0][1] - y);
	double s2 = (t[1][0] - x)*(t[1][1]-t[2][1])-(t[1][0]-t[2][0])*(t[1][1] - y);
	double s3 = (t[2][0] - x)*(t[2][1]-t[0][1])-(t[2][0]-t[0][0])*(t[2][1] - y);
	if((s1*s2 > 0)&&(s2*s3 > 0)) return true;
	if(s1 == 0 && x >= min(t[0][0],t[1][0]) && x <= max(t[0][0],t[1][0])) return true;
	if(s2 == 0 && x >= min(t[1][0],t[2][0]) && x <= max(t[1][0],t[2][0])) return true;
	if(s3 == 0 && x >= min(t[2][0],t[0][0]) && x <= max(t[2][0],t[0][0])) return true;
	return false;
}

bool LessV(Eigen::Vector4f& a, Eigen::Vector4f& b)
{
	return (a[0] == b[0])? a[1] < b[1] : a[0] < b[0];
}

void Sdf::coverPointinTriangle(int v[3],int** pts, int xl,int yl,Eigen::Matrix4f mc,double minx,double miny)
{
	int ix,iy,ax,ay;
	Eigen::Vector4f t[3];
	for(int i = 0;i < 3;i++)
	{
		GLVertex vv = GL_Mesh::gl_vertexes->at(v[i]);
		//将上一帧的点投影到本帧投影平面上
		Eigen::Vector4f vec;
		vec[0] = vv.v[0];
		vec[1] = vv.v[1];
		vec[2] = vv.v[2];
		vec[3] = 1;
		t[i] = mc*vec; //投影到当前帧坐标
	}
	int x0 = static_cast<int>((t[0][0]-minx)/devide);
	int y0 = static_cast<int>((t[0][1]-miny)/devide);
	ix = x0;
	iy = y0;
	ax = x0+1;
	ay = y0+1;
	for(int i = 1;i < 3;i++)
	{
		int xi = static_cast<int>((t[i][0]-minx)/devide);
	    int yi = static_cast<int>((t[i][1]-miny)/devide);
		if(xi < ix) ix = xi;
		if(yi < iy) iy = yi;
		if(xi+1 > ax) ax = xi+1;
		if(yi+1 > ay) ay = yi+1;
	}
	for(int i = ix;i <= ax;i++)
	{
		if(i < 0) continue;
		if(i >= xl) break;
		for(int j = iy; j<= ay;j++)
		{
			if(j < 0) continue;
		    if(j >= yl) break;
			double x = minx + i*devide;
			double y = miny + j*devide;
			if(inTriangle(t,x,y))
			{
				//if(!LessV(t[0],t[1]) && !!LessV(t[0],t[2]))  pts[i][j] = v[0];
				//else if(!LessV(t[1],t[0]) && !!LessV(t[1],t[2])) pts[i][j] = v[1];
				//else pts[i][j] = v[2];
				double dis1 = sqrt((t[0][0]-x)*(t[0][0]-x)+(t[0][1]-y)*(t[0][1]-y));
				double dis2 = sqrt((t[1][0]-x)*(t[1][0]-x)+(t[1][1]-y)*(t[1][1]-y));
				double dis3 = sqrt((t[2][0]-x)*(t[2][0]-x)+(t[2][1]-y)*(t[2][1]-y));
				if(dis1 <= dis2 && dis1 <= dis3) pts[i][j] = v[0];
				else if(dis2 <= dis1 && dis2 <= dis3) pts[i][j] = v[1];
				else pts[i][j] = v[2];
				//pts[i][j] = v[0];
			}
		}
	}
}

//返回值 是否存在缝合线
inline bool Sdf::buildMeshesAtR(int** points,int i,int j,Eigen::Matrix4f tran,int forbid)
{
	bool canbuild = false;
	//存储下三角面片
	if(points[i-1][j-1] >= 0 && points[i-1][j] >= 0)
	{
		if(points[i-1][j-1] >= forbid || points[i-1][j] >= forbid || points[i][j] >= forbid){
		
			//作为一阶邻域，各个顶点互相进行Laplacian滤波

			GLVertex& vv1 = GL_Mesh::gl_vertexes->at(points[i-1][j-1]);
			GLVertex& vv2 = GL_Mesh::gl_vertexes->at(points[i-1][j]);
			GLVertex& vv3 = GL_Mesh::gl_vertexes->at(points[i][j]);

			vv1.LaplacianUpdate(vv2.v);
			vv1.LaplacianUpdate(vv3.v);

			vv2.LaplacianUpdate(vv1.v);
			vv2.LaplacianUpdate(vv3.v);

			vv3.LaplacianUpdate(vv1.v);
			vv3.LaplacianUpdate(vv2.v);
			if(points[i-1][j-1] >= forbid && points[i-1][j] >= forbid && points[i][j] >= forbid)
			{
				TryaddMesh(points[i][j],points[i-1][j],points[i-1][j-1],tran);
			}
			else
			{
				//TryaddMesh(points[i][j],points[i-1][j],points[i-1][j-1],tran);
				//TryAddRightMesh(points[i][j],points[i-1][j],points[i-1][j-1],tran);
				canbuild = true;
			}
		}
	}

	//存储上三角面片
	if(points[i-1][j-1] >= 0 && points[i][j-1] >= 0)
	{
		if(points[i-1][j-1] >= forbid || points[i][j-1] >= forbid || points[i][j] >= forbid){

			//作为一阶邻域，各个顶点互相进行Laplacian滤波
			GLVertex& vv1 = GL_Mesh::gl_vertexes->at(points[i-1][j-1]);
			GLVertex& vv2 = GL_Mesh::gl_vertexes->at(points[i][j-1]);
			GLVertex& vv3 = GL_Mesh::gl_vertexes->at(points[i][j]);

			vv1.LaplacianUpdate(vv2.v);
			vv1.LaplacianUpdate(vv3.v);

			vv2.LaplacianUpdate(vv1.v);
			vv2.LaplacianUpdate(vv3.v);

			vv3.LaplacianUpdate(vv1.v);
			vv3.LaplacianUpdate(vv2.v);

			if(points[i-1][j-1] >= forbid && points[i][j-1] >= forbid && points[i][j] >= forbid)
			{
				TryaddMesh(points[i][j],points[i-1][j-1], points[i][j-1],tran);
			}
			else
			{
				//TryaddMesh(points[i][j],points[i-1][j-1], points[i][j-1],tran);
				//TryAddRightMesh(points[i][j],points[i-1][j-1], points[i][j-1],tran);
				canbuild = true;
			}
		}
	}

	//返回值 是否存在缝合线
	return canbuild;
}

//判断线段端点在直线的哪一侧.
double direction(ProjectPoint p1,ProjectPoint p2,ProjectPoint p3)
{
	return (p2.x-p1.x)*(p3.y-p1.y)-(p3.x-p1.x)*(p2.y-p1.y);
}

//判断与直线共线的点是否在线段上
bool onSegment(ProjectPoint p1,ProjectPoint p2,ProjectPoint p3)
{
	double left,right;

	if(p1.x<=p2.x)
	{
		left=p1.x;
		right=p2.x;
	}
	else
	{
		left=p2.x;
		right=p1.x;
	}

	if(left<=p3.x&&p3.x<=right)
	{
		return true;
	}
	else
	{
		return false;
	}
}




bool segmentsIntersect(ProjectPoint p1,ProjectPoint p2,ProjectPoint p3,ProjectPoint p4)
{
	double d1,d2,d3,d4; 
	/*
		以下2句作用：以p3p4所在直线为标准,判断线段p1p2是否跨越p3p4所在直线
		即判断p1在p3p4的一边，p2在p3p4另一边
	*/
	d1=direction(p3,p4,p1);
	d2=direction(p3,p4,p2);
	/*
		以下2句作用：同上，不过这次判断的是p3p4是否跨越了p1p2所在直线
	*/
	d3=direction(p1,p2,p3);
	d4=direction(p1,p2,p4);

	//互相跨越,必相交
	if(d1*d2<0&&d3*d4<0)
	{
		return true;
	}
	/*
		否则,看是否端点在另一个线段上
	*/
	else if(0==d1&&onSegment(p3,p4,p1)==true)
	{
		return true;
	}
	else if(0==d2&&onSegment(p3,p4,p2)==true)
	{
		return true;
	}
	else if(0==d3&&onSegment(p1,p2,p3)==true)
	{
		return true;
	}
	else if(0==d4&&onSegment(p1,p2,p4)==true)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool comparePP(ProjectPoint& a, ProjectPoint& b)
{
	return (abs(a.theta - b.theta) < 0.00000001)?a.radius < b.radius:a.theta<b.theta;
}

void Sdf::buildMeshesGreedy(std::set<int> nearest,int Rid,Eigen::Matrix4f tran,int forbid, int* except)
{
	if(nearest.size() <=0 ) return;
	//投影平面y轴
	GLVertex& R = GL_Mesh::gl_vertexes->at(Rid);
	ProjectPoint tempR;
	tempR.x = 0;
	tempR.y = 0;
	double mcos = (-tran(0,1)*R.v[3]-tran(1,1)*R.v[4]-tran(2,1)*R.v[5]);
	double my[3];
	my[0] = - tran(0,1) - mcos*R.v[3];
	my[1] = - tran(1,1) - mcos*R.v[4];
	my[2] = - tran(2,1) - mcos*R.v[5];
	double mylen = sqrt(my[0]*my[0]+my[1]*my[1]+my[2]*my[2]);
	my[0] /= mylen;
	my[1] /= mylen;
	my[2] /= mylen;
	//投影平面x轴
	double mx[3];
	mx[0] = my[1]*R.v[5]-my[2]*R.v[4];
	mx[1] = my[2]*R.v[3]-my[0]*R.v[5];
	mx[2] = my[0]*R.v[4]-my[1]*R.v[3];
	double mxlen = sqrt(mx[0]*mx[0]+mx[1]*mx[1]+mx[2]*mx[2]);
	mx[0] /= mxlen;
	mx[1] /= mxlen;
	mx[2] /= mxlen;

	//计算投影坐标并按极坐标排序
	std::vector<ProjectPoint> nearPro;
	for(std::set<int>::iterator it = nearest.begin();it != nearest.end();it++)
	{
		ProjectPoint pp;
		pp.vertexID = *it;
		pp.v = GL_Mesh::gl_vertexes->at(*it);
		//投影到近平面
		double distance = sqrt((pp.v.v[0]-R.v[0])*(pp.v.v[0]-R.v[0])+(pp.v.v[1]-R.v[1])*(pp.v.v[1]-R.v[1])+(pp.v.v[2]-R.v[2])*(pp.v.v[2]-R.v[2]));
		double cosAng = ((pp.v.v[0]-R.v[0])*R.v[3]+(pp.v.v[1]-R.v[1])*R.v[4]+(pp.v.v[2]-R.v[2])*R.v[5])/distance;
		//待完善：角度阈值
		//if(cosAng > M_PI/4) continue;
		//投影后的向量
		double t[3];
		t[0] = (pp.v.v[0]-R.v[0]) - cosAng*distance*R.v[3];
		t[1] = (pp.v.v[1]-R.v[1]) - cosAng*distance*R.v[4];
		t[2] = (pp.v.v[2]-R.v[2]) - cosAng*distance*R.v[5];
		pp.x = t[0]*mx[0]+t[1]*mx[1]+t[2]*mx[2];
		pp.y = t[0]*my[0]+t[1]*my[1]+t[2]*my[2];
		double tlen = sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
		pp.radius = tlen;
	    double the = acos(pp.x/tlen);
		if(pp.y >= 0) pp.theta = the;
		else pp.theta = 2*M_PI - the;

		
		//if(pp.v.degree >= 2*M_PI) pp.deleted = true;
	    if(pp.vertexID < forbid) pp.deleted = true;
		else pp.deleted = false;
		//pp.radius = sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
		//pp.costheta = (t[0]*mx[0]+t[1]*mx[1]+t[2]*mx[2])/pp.radius;
		nearPro.push_back(pp);
	}
	
	std::sort(nearPro.begin(),nearPro.end(),comparePP);
	std::set<std::pair<int,int>> crossEdges;
	std::vector<bool> block1;
	for(int i = 0;i < nearPro.size();i++)
	{
		block1.push_back(false);
		if(i > 0 && (abs(nearPro[i].theta - nearPro[i-1].theta) < 0.00000001)) nearPro[i].deleted = true;
	}
	//可见性检测
	
	for(int i = 0;i < nearPro.size();i++)
	{
		for(int j = i+1;j < nearPro.size();j++)
		{
			if(GL_Mesh::findEdge(nearPro[i].vertexID,nearPro[j].vertexID))
			{
				//判断夹角方向
				double theta = nearPro[j].theta-nearPro[i].theta;
				if(theta < 0) theta += 2*M_PI;
				if(GL_Mesh::findMesh(Rid,nearPro[i].vertexID,nearPro[j].vertexID))//GL_Mesh::findEdge(nearPro[i].vertexID,Rid) && GL_Mesh::findEdge(nearPro[j].vertexID,Rid))
				{
					//if(theta < 0) theta += 2*M_PI;
			        R.degree += theta;
					//阻断i、j间的所有点
					if(theta < M_PI)
					{
						block1[i] = true;
						for(int k = i+1;k <j;k++)
						{
							nearPro[k].deleted = true;
							block1[k] = true;
						}
					}else{
						block1[j] = true;
						for(int k = 0;k <i;k++)
						{
							nearPro[k].deleted = true;
							block1[k] = true;
						}
						for(int k = j+1;k <nearPro.size();k++)
						{
							nearPro[k].deleted = true;
							block1[k] = true;
						}
					}
				}else{
				    //阻断i、j间所有在比边ij远的点
					if(theta < M_PI)
					{
						for(int k = i+1;k <j;k++)
						{
							if(segmentsIntersect(nearPro[i],nearPro[j],tempR,nearPro[k])){
								nearPro[k].deleted = true;
							}
						}
					}else{
						for(int k = 0;k <i;k++)
						{
							if(segmentsIntersect(nearPro[i],nearPro[j],tempR,nearPro[k])){
								nearPro[k].deleted = true;
							}
						}
						for(int k = j+1;k <nearPro.size();k++)
						{
							if(segmentsIntersect(nearPro[i],nearPro[j],tempR,nearPro[k])){
								nearPro[k].deleted = true;
							}
						}
					}
				}
			}
		}
	}
	
	//获得所有可见点
	std::vector<ProjectPoint> canSee;
	std::vector<bool> block2;
	bool block = false;
	for(int i = 0;i < nearPro.size();i++)
	{
		if(!nearPro[i].deleted || nearPro[i].vertexID == except[0] || nearPro[i].vertexID == except[1]){
			canSee.push_back(nearPro[i]);
			block2.push_back(block);
			block = false;
		}
		if(block1[i]) block = true;
	}
	if(block && block2.size() > 0) block2[0] = true;
	//搭建三角面片
	for(int i = 0;i < canSee.size();i++)
	{
		int a = i;
		int b = i+1;
		if(i == canSee.size()-1) b = 0;
		//if(GL_Mesh::findEdge(canSee[a],Rid) && GL_Mesh::findEdge(canSee[b],Rid)){
		if(block2[b]){
			//三角形已存在
			//double theta = canSee[b].theta-canSee[a].theta;
			//if(theta < 0) theta += 2*M_PI;
			//R.degree += theta;
		}else{
			//搭建三角形
			double theta = canSee[b].theta-canSee[a].theta;
			if(theta < 0) theta += 2*M_PI;
			if(theta < M_PI){
				TryAddRightMesh(Rid,canSee[a].vertexID,canSee[b].vertexID,tran);
				R.degree += theta;
			}
		}
	}
}
void Sdf::castIsoMeshes2(std::vector<Eigen::Matrix4f> trans)
{
	//单帧镜头的参数
	double ix = -1.5;
	double iy = -1;
	double ax = 1.2;
	double ay = 1;
	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	double searchx = 0.0;
	double searchy = 0.0;


}
void Sdf::castIsoMeshes()
{
	//单帧镜头的参数
	double ix = -1.5;
	double iy = -1;
	double ax = 1.2;
	double ay = 1;

	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	int** points = NULL;
	for(int i = 0;i < GL_Mesh::gl_trans->size();i++)
	{
		int** newp = NULL;
		if(i > 0)
		  newp = castPieceMeshes(i,points,i-1);
		else
		  newp = castPieceMeshes(i,points,i);
		if(points != NULL)
		{
			for(int i = 0;i <xl;i++) delete[] points[i];
			delete[] points;
		}
		points = newp;
	}
}

void Sdf::castNextMeshes(int i)
{
	//单帧镜头的参数
	double ix = -1.5;
	double iy = -1;
	double ax = 1.2;
	double ay = 1;

	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	int** newp = NULL;
	if(i > 0)
		newp = castPieceMeshes(i,oldp,i-1);
	else
		newp = castPieceMeshes(i,oldp,i);
	if(oldp != NULL)
	{
		for(int i = 0;i <xl;i++) delete[] oldp[i];
		delete[] oldp;
	}
	oldp = newp;
}
void Sdf::castMeshes(Eigen::Matrix4f tran){
	double ix = 10000,iy = 10000,ax = -10000,ay = -10000;

	Eigen::Vector4f v;
	Eigen::Vector4f t;
	Eigen::Matrix4f mc =  tran.inverse();

	//求x,y值区间
	v[0] = minX;
	v[1] = minY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = minY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = maxY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = minX;
	v[1] = maxY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = minY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = minY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = maxY;
	v[2] = minZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	v[0] = maxX;
	v[1] = maxY;
	v[2] = maxZ;
	v[3] = 1;
	t = mc * v;
	if(t[0] < ix) ix = t[0];
	if(t[0] > ax) ax = t[0];
	if(t[1] < iy) iy = t[1];
	if(t[1] > ay) ay = t[1];

	//cout<<ix<<","<<iy<<endl;
	//对区间内进行采样
	int xl = (ax-ix) / devide;
	int yl = (ay-iy) / devide;

	//ix = -1.5;
	//iy = -1;
	//xl = 270;
	//yl = 200;
	//缓存结构，存储本行与上一行的表面点
	//int* lastline = NULL;
	//int* nowline = NULL;
	int ** planePoints = NULL;
	if(xl > 0 ) planePoints = new int*[xl];
	for(int i = 0;i < xl;i++)
	{
		planePoints[i] = new int[yl];
		//nowline = new int[yl];
		for(int j = 0;j < yl;j++)
		{
			DistandNorm dnn = getDepth(ix + i*devide,iy + j*devide,tran);
			if(dnn.dist > -10000 && !isnan(dnn.n_x)){
				GLVertex v1 = *(new GLVertex());
				putPoint(v1,ix+i*devide,iy+j*devide,dnn,tran);
				int idx = (v1.v[0] - minX)/devide;
				int idy = (v1.v[1] - minY)/devide;
				int idz = (v1.v[2] - minZ)/devide;
				std::map<std::tuple<int,int,int>,GLMergeVertex>::iterator it = GL_Mesh::gl_vertex_map->find(std::tuple<int,int,int>(idx,idy,idz));
				if(it!= GL_Mesh::gl_vertex_map->end())
				{
					GLVertex& vr = GL_Mesh::gl_vertexes->at(it->second.vid);
					for(int lp = 0;lp < 6;lp++)
					{
						vr.v[lp] = (vr.v[lp]*it->second.weight + v1.v[lp])/(it->second.weight+1);
					}
					it->second.weight++;
					planePoints[i][j] =  it->second.vid;
				}else{
					GL_Mesh::gl_vertexes->push_back(v1);
					//nowline[j] = GL_Mesh::gl_vertexes->size()-1;
					planePoints[i][j] =  GL_Mesh::gl_vertexes->size()-1;
					GLMergeVertex mv;
					mv.vid = planePoints[i][j];
					mv.weight = 1;
					std::pair<std::tuple<int,int,int>,GLMergeVertex> ins(std::tuple<int,int,int>(idx,idy,idz),mv);
					GL_Mesh::gl_vertex_map->insert(ins);
				}
			}else{
				//nowline[j] = -1;
				planePoints[i][j] = -1;
			}
			//nowline[j] = getDepth(ix + i*devide,iy + j*devide,tran);
			if(i > 0 && j > 0 && planePoints[i][j] >= 0)
			{
				buildMeshesAtR(planePoints,i,j,tran);
				//buildMeshesGreedy(planePoints,i,j,xl,yl,tran);
			} 
		}
	}

	GL_Mesh::updateVertexes();
}

inline void Sdf::buildMeshesAtR(int** points,int i,int j,Eigen::Matrix4f tran)
{
	//存储下三角面片
	if(points[i-1][j-1] >= 0 && points[i-1][j] >= 0)
	{
		TryaddMesh(points[i][j],points[i-1][j],points[i-1][j-1],tran);
	}

	//存储上三角面片
	if(points[i-1][j-1] >= 0 && points[i][j-1] >= 0)
	{
		TryaddMesh(points[i][j],points[i-1][j-1], points[i][j-1],tran);
	}
}

void Sdf::putPoint(GLVertex& m,double x,double y,DistandNorm dnn,Eigen::Matrix4f tran)
{
	Eigen::Matrix4f mc =  tran;
	m.v[0] = mc(0,0)*x + mc(0,1)*y + mc(0,2)*dnn.dist + mc(0,3);
	m.v[1] = mc(1,0)*x + mc(1,1)*y + mc(1,2)*dnn.dist + mc(1,3);
	m.v[2] = mc(2,0)*x + mc(2,1)*y + mc(2,2)*dnn.dist + mc(2,3);
	m.v[3] = dnn.n_x;
	m.v[4] = dnn.n_y;
	m.v[5] = dnn.n_z;
	m.dist = dnn.dist;
}

bool Sdf::validMesh(GLMesh& m,Eigen::Matrix4f tran)
{
	GLVertex v[3];
	v[0] = GL_Mesh::gl_vertexes->at(m.vertexes[0]);
	v[1] = GL_Mesh::gl_vertexes->at(m.vertexes[1]);
	v[2] = GL_Mesh::gl_vertexes->at(m.vertexes[2]);
    double th = devide*10;
	double d1 = sqrt((v[1].v[0] - v[0].v[0])*(v[1].v[0] - v[0].v[0])+(v[1].v[1] - v[0].v[1])*(v[1].v[1] - v[0].v[1])+(v[1].v[2] - v[0].v[2])*(v[1].v[2] - v[0].v[2]));
	double d2 = sqrt((v[2].v[0] - v[0].v[0])*(v[2].v[0] - v[0].v[0])+(v[2].v[1] - v[0].v[1])*(v[2].v[1] - v[0].v[1])+(v[2].v[2] - v[0].v[2])*(v[2].v[2] - v[0].v[2]));
	double d3 = sqrt((v[1].v[0] - v[2].v[0])*(v[1].v[0] - v[2].v[0])+(v[1].v[1] - v[2].v[1])*(v[1].v[1] - v[2].v[1])+(v[1].v[2] - v[2].v[2])*(v[1].v[2] - v[2].v[2]));
	
	if(d1 > th) return false;
	if(d2 > th) return false;
	if(d3 > th) return false;

	//double angle1 =  ((v[1].v[0] - v[0].v[0])*(v[2].v[0] - v[0].v[0])+(v[1].v[1] - v[0].v[1])*(v[2].v[1] - v[0].v[1])+(v[1].v[2] - v[0].v[2])*(v[2].v[2] - v[0].v[2]))/d1/d2;
	//double angle2 =  ((v[0].v[0] - v[1].v[0])*(v[2].v[0] - v[1].v[0])+(v[0].v[1] - v[1].v[1])*(v[2].v[1] - v[1].v[1])+(v[0].v[2] - v[1].v[2])*(v[2].v[2] - v[1].v[2]))/d1/d3;
	//double angle3 =  ((v[1].v[0] - v[2].v[0])*(v[0].v[0] - v[2].v[0])+(v[1].v[1] - v[2].v[1])*(v[0].v[1] - v[2].v[1])+(v[1].v[2] - v[2].v[2])*(v[0].v[2] - v[2].v[2]))/d2/d3;

	//if(angle1 > 0.995) return false;
	//if(angle2 > 0.995) return false;
	//if(angle3 > 0.995) return false;
	
	double a[3];
	double b[3];
	double c[3];
	double d[3];
	a[0] = v[1].v[0] - v[0].v[0];
	a[1] = v[1].v[1] - v[0].v[1];
	a[2] = v[1].v[2] - v[0].v[2];
	b[0] = v[2].v[0] - v[0].v[0];
	b[1] = v[2].v[1] - v[0].v[1];
	b[2] = v[2].v[2] - v[0].v[2];

	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
	double clen = sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);

	d[0] = tran(0,2);
	d[1] = tran(1,2);
	d[2] = tran(2,2);
	double dlen = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

	double angle = (c[0]*d[0]+c[1]*d[1]+c[2]*d[2])/clen/dlen;

 	if(abs(angle) < 0.05) return false;

	//d[0] = v[0].v[3];
	//d[1] = v[0].v[4];
	//d[2] = v[0].v[5];
	//dlen = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);

    //angle = (c[0]*d[0]+c[1]*d[1]+c[2]*d[2])/clen/dlen;

 	//if(abs(angle) < 0.05) return false;
	
	return true;
}

DistandNorm Sdf::getDepth(double x,double y,Eigen::Matrix4f tran)
{
	DistandNorm ret;
	ret.dist = -10000;
	//求z值区间
	Eigen::Matrix4f mc =  tran;
	double min = -10000,max = 10000;
	if(mc(0,2) > 0)
	{
		double left = (minX - mc(0,0)*x - mc(0,1)*y - mc(0,3))/mc(0,2);
		if(left > min) min = left;
		double right = (maxX - mc(0,0)*x - mc(0,1)*y - mc(0,3))/mc(0,2);
		if(right < max) max = right;
	}else if(mc(0,2) < 0)
	{
		double left = (minX - mc(0,0)*x - mc(0,1)*y - mc(0,3))/mc(0,2);
		if(left < max) max = left;
		double right = (maxX - mc(0,0)*x - mc(0,1)*y - mc(0,3))/mc(0,2);
		if(right > min) min = right;
	}
	if(mc(1,2) > 0)
	{
		double left = (minY - mc(1,0)*x - mc(1,1)*y - mc(1,3))/mc(1,2);
		if(left > min) min = left;
		double right = (maxY - mc(1,0)*x - mc(1,1)*y - mc(1,3))/mc(1,2);
		if(right < max) max = right;
	}else if(mc(1,2) < 0)
	{
		double left = (minY - mc(1,0)*x - mc(1,1)*y - mc(1,3))/mc(1,2);
		if(left < max) max = left;
		double right = (maxY - mc(1,0)*x - mc(1,1)*y - mc(1,3))/mc(1,2);
		if(right > min) min = right;
	}
	if(mc(2,2) > 0)
	{
		double left = (minZ - mc(2,0)*x - mc(2,1)*y - mc(2,3))/mc(2,2);
		if(left > min) min = left;
		double right = (maxZ - mc(2,0)*x - mc(2,1)*y - mc(2,3))/mc(2,2);
		if(right < max) max = right;
	}else if(mc(2,2) < 0)
	{
		double left = (minZ - mc(2,0)*x - mc(2,1)*y - mc(2,3))/mc(2,2);
		if(left < max) max = left;
		double right = (maxZ - mc(2,0)*x - mc(2,1)*y - mc(2,3))/mc(2,2);
		if(right > min) min = right;
	}

	//z存在取值
	if(min <= max)
	{
		//获得第一个格子
		double fx = mc(0,0)*x + mc(0,1)*y + mc(0,2)*min + mc(0,3);
		double fy = mc(1,0)*x + mc(1,1)*y + mc(1,2)*min + mc(1,3);
		double fz = mc(2,0)*x + mc(2,1)*y + mc(2,2)*min + mc(2,3);

		double originx = mc(0,0)*x + mc(0,1)*y + mc(0,3);
		double originy = mc(1,0)*x + mc(1,1)*y + mc(1,3);
		double originz = mc(2,0)*x + mc(2,1)*y + mc(2,3);
		
		double d = 0;
		int ox = -1;
		int oy = -1;
		int oz = -1;
		double oldd = 0;

		double step = devide; //步长

		while(d <= max-min)
		{
			//扫描
			int xid = (fx-minX)/devide;
			int yid = (fy-minY)/devide;
			int zid = (fz-minZ)/devide;

			int id = ((zid*Ylen)+yid)*Xlen+xid;
			if(xid < 0 || xid >= Xlen || yid < 0 || yid >= Ylen || zid < 0 || zid >= Zlen)
			{
				ox = -1;
				oldd = 0;
			}
			else if(weight[id] <= 0)
			{
				ox = -1;
				oldd = 0;
			}
			else if(dist[id] <= 0)
			{
				ox = xid;
				oy = yid;
				oz = zid;
				oldd = d;
			}else if(dist[id] > 0 && ox >= 0){
				

				int oid = ((oz*Ylen)+oy)*Xlen+ox;
				//将两个格子的中心分别投影到扫描线上
				double d1 = getProjectd(minX+ox*devide,minY+oy*devide,minZ+oz*devide,originx,originy,originz,n_x[oid],n_y[oid],n_z[oid],mc(0,2),mc(1,2),mc(2,2));
                if(abs(oldd-d1) > devide) d1 = oldd; 
				double d2 = getProjectd(minX+xid*devide,minY+yid*devide,minZ+zid*devide,originx,originy,originz,n_x[id],n_y[id],n_z[id],mc(0,2),mc(1,2),mc(2,2));
				if(abs(d-d2) > devide) d2 = d;
				//double d1 = oldd;
				//double d2 = d;

				double rate = (dist[id]-dist[oid])*(-dist[oid]);
				double zerocross = d1+ rate*(d2-d1);
				

				/*
				double d1 = oldd;
				double d2 = d;
				double f1 = getTriInter(originx+oldd*mc(0,2),originy+oldd*mc(1,2),originz+oldd*mc(2,2));
				double f2 = getTriInter(fx,fy,fz);

				double rate = f1/(f2-f1);
				double zerocross = d1- rate*(d2-d1);
				*/

				ret.dist = min+zerocross;
				ret.n_x = n_x[id];
				ret.n_y = n_y[id];
				ret.n_z = n_z[id];
				break;
			}else{
				ox = -1;
				oldd = 0;
			}

			d += step;
			fx+=mc(0,2)*step;
			fy+=mc(1,2)*step;
			fz+=mc(2,2)*step;
		}	
	}
	return ret;
}

double Sdf::getProjectd(double x,double y,double z,double ox,double oy,double oz,double nx,double ny,double nz,double mx,double my,double mz)
{
	return ((x-ox)*nx + (y-oy)*ny + (z-oz)*nz)/(nx*mx+ny*my+nz*mz);
}


void Sdf::outputObj()
{
	
	std::ofstream output("out.obj");

	for(std::vector<GLVertex>::iterator it = GL_Mesh::gl_vertexes->begin();it != GL_Mesh::gl_vertexes->end();it++)
	{
		GLVertex v = *it;
		output<<"v "<<v.v[0]<<" "<<v.v[1]<<" "<<v.v[2]<<std::endl;
	}
	for(std::vector<GLVertex>::iterator it = GL_Mesh::gl_vertexes->begin();it != GL_Mesh::gl_vertexes->end();it++)
	{
		GLVertex v = *it;
		output<<"vn "<<v.v[3]<<" "<<v.v[4]<<" "<<v.v[5]<<std::endl;
	}
	for(std::map<std::tuple<int,int,int>,GLMesh>::iterator it = GL_Mesh::gl_meshes->begin();it != GL_Mesh::gl_meshes->end();it++)
	{
		GLMesh m = it->second;
		output<<"f "<<m.vertexes[0]+1<<"//"<<m.vertexes[0]+1<<" "<<m.vertexes[1]+1<<"//"<<m.vertexes[1]+1<<" "<<m.vertexes[2]+1<<"//"<<m.vertexes[2]+1<<std::endl;
	}
	output.close();

	/*
	try
	{
		if ( !OpenMesh::IO::write_mesh(GL_Mesh::o_meshes, "output.obj") )
		{
			std::cerr << "Cannot write mesh to file 'output.off'" << std::endl;
		}
	}
	catch( std::exception& x )
	{
		std::cerr << x.what() << std::endl;
	}*/
}

void Sdf::computeNormal(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pts,pcl::PointCloud<pcl::Normal>::Ptr normals)
{
   pcl::PointCloud<pcl::PointXYZRGB>::Ptr tangent_x(new pcl::PointCloud<pcl::PointXYZRGB>);
   pcl::PointCloud<pcl::PointXYZRGB>::Ptr tangent_y(new pcl::PointCloud<pcl::PointXYZRGB>);

   pcl::PointCloud<pcl::PointXYZRGB>::Ptr integral_x(new pcl::PointCloud<pcl::PointXYZRGB>);
   pcl::PointCloud<pcl::PointXYZRGB>::Ptr integral_y(new pcl::PointCloud<pcl::PointXYZRGB>);

   //x,y方向切向量图
   tangent_x->height = pts->height;
   tangent_x->width = pts->width;
   tangent_x->is_dense = false;
   tangent_x->points.resize(tangent_x->height * tangent_x->width);
   pcl::PointXYZRGB lastV,nowV;
   for(int i = 0;i < pts->height;i++){
	   for(int j = 0;j < pts->width;j++){
		   lastV.x = 0;
		   lastV.y = 0;
		   lastV.z = 0;
		   nowV.x = 0;
		   nowV.y = 0;
		   nowV.z = 0;
		   if(j > 0){
			   lastV.x = pts->at(j,i).x - pts->at(j-1,i).x;
		       lastV.y = pts->at(j,i).y - pts->at(j-1,i).y;
		       lastV.z = pts->at(j,i).z - pts->at(j-1,i).z;
		   }
		   if(j < pts->width-1){
			   nowV.x = pts->at(j+1,i).x - pts->at(j,i).x;
		       nowV.y = pts->at(j+1,i).y - pts->at(j,i).y;
		       nowV.z = pts->at(j+1,i).z - pts->at(j,i).z;
		   }
		   tangent_x->at(j,i).x = lastV.x + nowV.x;
		   tangent_x->at(j,i).y = lastV.y + nowV.y;
		   tangent_x->at(j,i).z = lastV.z + nowV.z;
	   }
   }

   tangent_y->height = pts->height;
   tangent_y->width = pts->width;
   tangent_y->is_dense = false;
   tangent_y->points.resize(tangent_y->height * tangent_y->width);
    for(int i = 0;i < pts->width;i++){
	   for(int j = 0;j < pts->height;j++){
		   lastV.x = 0;
		   lastV.y = 0;
		   lastV.z = 0;
		   nowV.x = 0;
		   nowV.y = 0;
		   nowV.z = 0;
		   if(j > 0){
			   lastV.x = pts->at(i,j).x - pts->at(i,j-1).x;
		       lastV.y = pts->at(i,j).y - pts->at(i,j-1).y;
		       lastV.z = pts->at(i,j).z - pts->at(i,j-1).z;
		   }
		   if(j < pts->height-1){
			   nowV.x = pts->at(i,j+1).x - pts->at(i,j).x;
		       nowV.y = pts->at(i,j+1).y - pts->at(i,j).y;
		       nowV.z = pts->at(i,j+1).z - pts->at(i,j).z;
		   }
		   tangent_y->at(i,j).x = lastV.x + nowV.x;
		   tangent_y->at(i,j).y = lastV.y + nowV.y;
		   tangent_y->at(i,j).z = lastV.z + nowV.z;
	   }
   }

   
   integral_x->height = pts->height;
   integral_x->width = pts->width;
   integral_x->is_dense = false;
   integral_x->points.resize(integral_x->height * integral_x->width);
   integral_y->height = pts->height;
   integral_y->width = pts->width;
   integral_y->is_dense = false;
   integral_y->points.resize(integral_y->height * integral_y->width);
   for(int i = 0;i < pts->width;i++){
	   for(int j = 0;j < pts->height;j++){
		   if(i == 0 || j == 0){
			   integral_x->at(i,j).x = tangent_x->at(i,j).x;
			   integral_x->at(i,j).y = tangent_x->at(i,j).y;
			   integral_x->at(i,j).z = tangent_x->at(i,j).z;
			   integral_y->at(i,j).x = tangent_y->at(i,j).x;
			   integral_y->at(i,j).y = tangent_y->at(i,j).y;
			   integral_y->at(i,j).z = tangent_y->at(i,j).z;
		   }else{
			   integral_x->at(i,j).x = integral_x->at(i-1,j).x + integral_x->at(i,j-1).x - integral_x->at(i-1,j-1).x + tangent_x->at(i,j).x;
			   integral_x->at(i,j).y = integral_x->at(i-1,j).y + integral_x->at(i,j-1).y - integral_x->at(i-1,j-1).y + tangent_x->at(i,j).y;
			   integral_x->at(i,j).z = integral_x->at(i-1,j).z + integral_x->at(i,j-1).z - integral_x->at(i-1,j-1).z + tangent_x->at(i,j).z;
			   integral_y->at(i,j).x = integral_y->at(i-1,j).x + integral_y->at(i,j-1).x - integral_y->at(i-1,j-1).x + tangent_y->at(i,j).x;
			   integral_y->at(i,j).y = integral_y->at(i-1,j).y + integral_y->at(i,j-1).y - integral_y->at(i-1,j-1).y + tangent_y->at(i,j).y;
			   integral_y->at(i,j).z = integral_y->at(i-1,j).z + integral_y->at(i,j-1).z - integral_y->at(i-1,j-1).z + tangent_y->at(i,j).z;
		   }	   
	   }
   }

   normals->height = pts->height;
   normals->width = pts->width;
   normals->is_dense = false;
   normals->points.resize(normals->height * normals->width);
   for(int i = 0;i < pts->width;i++){
	   for(int j = 0;j < pts->height;j++){
		   int left = i - 5;
		   if(left < 0) left = 0;
		   int up = j - 5;
		   if(up < 0) up = 0;
		   int right = i + 5;
		   if(right > pts->width-1) right = pts->width-1;
		   int down = j + 5;
		   if(down > pts->height-1) down = pts->height-1;

		   double x1 =  integral_x->at(right,down).x - integral_x->at(left,down).x - integral_x->at(right,up).x + integral_x->at(left,up).x;
		   double y1 =  integral_x->at(right,down).y - integral_x->at(left,down).y - integral_x->at(right,up).y + integral_x->at(left,up).y;
		   double z1 =  integral_x->at(right,down).z - integral_x->at(left,down).z - integral_x->at(right,up).z + integral_x->at(left,up).z;
		   double x2 =  integral_y->at(right,down).x - integral_y->at(left,down).x - integral_y->at(right,up).x + integral_y->at(left,up).x;
		   double y2 =  integral_y->at(right,down).y - integral_y->at(left,down).y - integral_y->at(right,up).y + integral_y->at(left,up).y;
		   double z2 =  integral_y->at(right,down).z - integral_y->at(left,down).z - integral_y->at(right,up).z + integral_y->at(left,up).z;

		   double nx = y1*z2 - y2*z1;
		   double ny = z1*x2 - z2*x1;
		   double nz = x1*y2 - x2*y1;

		   double n = sqrt(nx*nx+ny*ny+nz*nz);

		   normals->at(i,j).normal_x = nx/n;
		   normals->at(i,j).normal_y = ny/n;
		   normals->at(i,j).normal_z = nz/n;
	   }
   }
}
//pcl::PointXYZ Sdf::getZeroCross(double a,double b,double c)
//{
	//return NULL;
//}

Sdf* Sdf::sdf = NULL;
TsdfNode* Sdf::tsdf = NULL;
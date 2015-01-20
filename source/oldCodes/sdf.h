#ifndef SIGNED_DISTANCE_FUNCTION_H
#define SIGNED_DISTANCE_FUNCTION_H

#include "CIsoSurface.h"
#include <XnCppWrapper.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include "glwindows.h"


struct DistandNorm
{
	double dist;
	double n_x;
	double n_y;
	double n_z;
};

class DepthCalculator
{
  public:
	  pcl::PointCloud<pcl::Normal>::Ptr normals;
	  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud;
	  xn::DepthGenerator* depthGenerator;
	  Eigen::Matrix4f tran;

      DistandNorm getDist(double x,double y,double z);
};



class TsdfNode{
	
public:
    TsdfNode();

	//最小分裂阈值
	const static double minDiv;

    //坐标位置
	double minX;
	double minY;
	double minZ;

	double maxX;
	double maxY;
	double maxZ;
	 
	//空间分割
	int Xlen;
	int Ylen;
	int Zlen;
	double devide;

	//最小距离阈值
	double epsi;
	//最大距离阈值
	double thre;

	//本身的距离值
	double dist;
	//本身的法向
	double n_x;
	double n_y;
	double n_z;
	//本身的权重值
	double weight;

	//体素存储结构
	TsdfNode*** child;

	//初始化
	void Initial(double ix,double iy,double iz,double ax,double ay,double az,double dev);
	//分裂，当距离小于分裂阈值时则分裂
	void split();
	//遍历更新
	void Traversal(DepthCalculator* dg);
	void TraversalDraw();

	double getWeight(double depth);
    void UpdateChild(int x,int y,int z,DistandNorm dnn);
};

class Sdf{
public:
  Sdf();
  Sdf(Eigen::Matrix4f ltran, double maxx,double maxy, double maxz,double minx,double miny,double minz,double dev,double th,double ep);
  void changeLocal(Eigen::Matrix4f ltran, double maxx,double maxy, double maxz,double minx,double miny,double minz);

  //SDF的旋转矩阵
  Eigen::Matrix4f local_tran;

  //生成的网格
  typedef boost::shared_ptr<GL_Mesh> MeshPtr;
  MeshPtr mesh;

  //空间范围
  double maxX;
  double maxY;
  double maxZ;
  double minX;
  double minY;
  double minZ;
  int Xlen;
  int Ylen;
  int Zlen;

  //分割参数
  double devide;
  double thre;
  double epsi;

  //网格构建参数
  int mu;
  double surfaceAng;  
  double minAng;
  double maxAng;


  //体素存储结构
  float* dist;
  float* n_x;
  float* n_y;
  float* n_z;
  float* weight;

  //临时存储提取点云
  int** oldp;

  double getWeight(double depth);
  void Update(int x,int y,int z,DistandNorm dnn,double weight);
  void Traversal(DepthCalculator* dg);
  void Traversal2(DepthCalculator* dg);
  void getData();
  void freeData();
  void TraversalDraw();
  void GenerateMeshes();

  pcl::PointCloud<pcl::PointXYZRGB>::Ptr Sdf::castCloud(Eigen::Matrix4f tran);
  void cloudToMeshes(Eigen::Matrix4f tran);
  void computeNormal(pcl::PointCloud<pcl::PointXYZRGB>::Ptr pts,pcl::PointCloud<pcl::Normal>::Ptr normals);

  int** castPieceMeshes(int tranID,int** pre,int lastID);
  void castNextMeshes(int i);
  void castIsoMeshes();
  void castIsoMeshes2(std::vector<Eigen::Matrix4f> trans);
  void castMeshes(Eigen::Matrix4f tran);
  void findNearest(int** planePoints,int xid,int yid,int xl,int yl,int Rid,std::set<int>& nearest,int m,int forbid);
  void findNearestAd(int** planePoints,int xid,int yid,int xl,int yl,int Rid,std::set<int>& nearest,int m,int forbid);
  void buildTic(int px,int py,double ix,double iy,double xl,double yl,int** pre,int** planePoints,Eigen::Matrix4f mc, Eigen::Matrix4f tran,int forbid,int* except,int m);
  inline void buildMeshesAtR(int** points,int i,int j,Eigen::Matrix4f tran);
  inline bool buildMeshesAtR(int** points,int i,int j,Eigen::Matrix4f tran,int forbid);
  void coverPointinMeshes(int** pre,int** pts,int i, int j, int xl,int yl,Eigen::Matrix4f mc,double minx,double miny);
  void coverPointinTriangle(int v[3],int** pts, int xl,int yl,Eigen::Matrix4f mc,double minx,double miny);
  void buildMeshesGreedy(std::set<int> nearest,int Rid,Eigen::Matrix4f tran,int forbid, int* except);
  DistandNorm getDepth(double x,double y,Eigen::Matrix4f tran);
  void putPoint(GLVertex& m,double x,double y,DistandNorm dnn,Eigen::Matrix4f tran);
  void TryaddMesh(int v1,int v2,int v3,Eigen::Matrix4f tran);
  void TryAddRightMesh(int a,int b,int c,Eigen::Matrix4f tran);
  bool validMesh(GLMesh& m,Eigen::Matrix4f tran);
  double getProjectd(double x,double y,double z,double ox,double oy,double oz,double nx,double ny,double nz,double mx,double my,double mz);
  void outputObj();

  CIsoSurface<double>* gen;

  static Sdf* sdf;
  static TsdfNode* tsdf;
};

struct ProjectPoint{
	int vertexID;
	GLVertex v;
	double x; //投影平面x坐标
	double y; //投影平面y坐标
	double theta; //投影平面极角
	double radius;
	//double costheta;
	bool deleted;
};

#endif
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

	//��С������ֵ
	const static double minDiv;

    //����λ��
	double minX;
	double minY;
	double minZ;

	double maxX;
	double maxY;
	double maxZ;
	 
	//�ռ�ָ�
	int Xlen;
	int Ylen;
	int Zlen;
	double devide;

	//��С������ֵ
	double epsi;
	//��������ֵ
	double thre;

	//����ľ���ֵ
	double dist;
	//����ķ���
	double n_x;
	double n_y;
	double n_z;
	//�����Ȩ��ֵ
	double weight;

	//���ش洢�ṹ
	TsdfNode*** child;

	//��ʼ��
	void Initial(double ix,double iy,double iz,double ax,double ay,double az,double dev);
	//���ѣ�������С�ڷ�����ֵʱ�����
	void split();
	//��������
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

  //SDF����ת����
  Eigen::Matrix4f local_tran;

  //���ɵ�����
  typedef boost::shared_ptr<GL_Mesh> MeshPtr;
  MeshPtr mesh;

  //�ռ䷶Χ
  double maxX;
  double maxY;
  double maxZ;
  double minX;
  double minY;
  double minZ;
  int Xlen;
  int Ylen;
  int Zlen;

  //�ָ����
  double devide;
  double thre;
  double epsi;

  //���񹹽�����
  int mu;
  double surfaceAng;  
  double minAng;
  double maxAng;


  //���ش洢�ṹ
  float* dist;
  float* n_x;
  float* n_y;
  float* n_z;
  float* weight;

  //��ʱ�洢��ȡ����
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
	double x; //ͶӰƽ��x����
	double y; //ͶӰƽ��y����
	double theta; //ͶӰƽ�漫��
	double radius;
	//double costheta;
	bool deleted;
};

#endif
#include "g2o/core/jacobian_workspace.h"
#include "g2o/stuff/macros.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/types/slam3d/edge_se3_prior.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include "g2o/types/slam3d/parameter_se3_offset.h"
#include "g2o/types/slam3d/vertex_se3.h"
//#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/pcg/linear_solver_pcg.h"
#include "g2o/solvers/dense/linear_solver_dense.h"

#include "slam3D.h"

using namespace std;
using namespace g2o;

typedef g2o::BlockSolver< g2o::BlockSolverTraits<6, 3> >  SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearCSparseSolver;
typedef g2o::LinearSolverDense<SlamBlockSolver::PoseMatrixType> SlamLinearDenseSolver;
//typedef g2o::LinearSolverCholmod<SlamBlockSolver::PoseMatrixType> SlamLinearCholmodSolver;
typedef g2o::LinearSolverPCG<SlamBlockSolver::PoseMatrixType> SlamLinearPCGSolver;
double zeroest[7] = {0.0,0.0,0.0,0.0,0.0,0.0,1.0};

GraphManager::GraphManager(){
	optimizer = new g2o::SparseOptimizer();
	
	SlamBlockSolver* solver = NULL;

	SlamLinearCSparseSolver* linearSolver = new SlamLinearCSparseSolver();
    linearSolver->setBlockOrdering(false);
    solver = new SlamBlockSolver(linearSolver);
	OptimizationAlgorithmDogleg * algo = new g2o::OptimizationAlgorithmDogleg(solver);
    optimizer->setAlgorithm(algo);

	optimizer->setVerbose(true);
	/*
    ParameterSE3Offset* sensorOffset = new ParameterSE3Offset;
	sensorOffset->setOffset(...);
	sensorOffset->setId(0);
	optimizer->addParameter(sensorOffset);
	*/
}
GraphManager::~GraphManager(){
	optimizer->clear();
	delete optimizer;
}

void GraphManager::addVertex(int id){
	
    VertexSE3* robot =  new VertexSE3;
    robot->setId(id);
	robot->setEstimateDataImpl(zeroest);
    optimizer->addVertex(robot);
}

void GraphManager::addEdge(int from,int to,Eigen::Matrix4f* trans,double variance,int inliers){
	EdgeSE3* odometry = new EdgeSE3;
    odometry->vertices()[0] = optimizer->vertex(from);
    odometry->vertices()[1] = optimizer->vertex(to);
	g2o::SE3Quat meanest = getEstimate(trans->cast<double>());
    odometry->setMeasurement(meanest);
    odometry->setInformation(Eigen::Matrix<double,6,6>::Identity()*inliers/variance); //
    optimizer->addEdge(odometry);
}

Eigen::Matrix4f* GraphManager::getEuleMatrix(double* est){

	float x = est[3];
	float y = est[4];
	float z = est[5];
	float w = est[6];

	float x2 = x * x;
	float y2 = y * y;
	float z2 = z * z;
	float xy = x * y;
	float xz = x * z;
	float yz = y * z;
	float wx = w * x;
	float wy = w * y;
	float wz = w * z;

	Eigen::Matrix4f* ret = new Eigen::Matrix4f;
	ret->data()[0] = 1.0f - 2.0f * (y2 + z2);
	ret->data()[1]  = 2.0f * (xy - wz);
	ret->data()[2]  = 2.0f * (xz + wy);
	ret->data()[3]  = 0.0f;

	ret->data()[4]  = 2.0f * (xy + wz);
	ret->data()[5] = 1.0f - 2.0f * (x2 + z2);
	ret->data()[6] = 2.0f * (yz - wx);
	ret->data()[7] = 0.0f;

	ret->data()[8] = 2.0f * (xz - wy);
	ret->data()[9] = 2.0f * (yz + wx);
	ret->data()[10] = 1.0f - 2.0f * (x2 + y2);
	ret->data()[11] = 0.0f;

	ret->data()[12] = est[0];
	ret->data()[13] = est[1];
	ret->data()[14] = est[2];
	ret->data()[15] = 1.0f;

	return ret;
}

g2o::SE3Quat GraphManager::getEstimate(const Eigen::Matrix4d& m){
  Eigen::Affine3d eigen_transform(m);
  Eigen::Quaterniond eigen_quat(eigen_transform.rotation());
  Eigen::Vector3d translation(m(0, 3), m(1, 3), m(2, 3));
  g2o::SE3Quat result(eigen_quat, translation);

  return result;
/*
   double w,x,y,z;
   float tr, s, q[4];
   int i, j, k;
   
   int nxt[3] = {1, 2, 0 };
   // 计算矩阵轨迹
   tr = (*m)(0,0) + (*m)(1,1) + (*m)(2,2);
   
   // 检查矩阵轨迹是正还是负
   if(tr>0.0f)
   {
    s = sqrt(tr + 1.0f);
    w = s / 2.0f;
    s = 0.5f / s;
    x = ((*m)(2,1) - (*m)(1,2)) * s;
    y = ((*m)(0,2) - (*m)(2,0)) * s;
    z = ((*m)(1,0) - (*m)(0,1)) * s;
   }
   else
   {
    // 轨迹是负
    // 寻找m11 m22 m33中的最大分量
    i = 0;
    if((*m)(1,1)>(*m)(0,0)) i = 1;
    if((*m)(2,2)>(*m)(i,i)) i = 2;
    j = nxt[i];
    k = nxt[j];
    
    s = sqrt(((*m)(i,i) - ((*m)(j,j) + (*m)(k,k))) + 1.0f);
    q[i] = s * 0.5f;
    if( s!= 0.0f) s = 0.5f / s;
    q[3] = ((*m)(k,j) - (*m)(j,k)) * s;
    q[j] = ((*m)(i,j) + (*m)(j,i)) * s;
    q[k] = ((*m)(i,k) + (*m)(k,i)) * s;
    x = q[0];
    y = q[1];
    z = q[2];
    w = q[3];	
   }
    double* ret = new double[7];
	ret[0] = m->data()[12];
	ret[1] = m->data()[13];
	ret[2] = m->data()[14];
	ret[3] = x;
	ret[4] = y;
	ret[5] = z;
	ret[6] = w;
	return ret;
	*/
}

void GraphManager::initializeOptimization(){
	//optimizer->save("tutorial_before.g2o");

	// prepare and run the optimization
	// fix the first robot pose to account for gauge freedom
	VertexSE3* firstRobotPose = dynamic_cast<VertexSE3*>(optimizer->vertex(0));
	firstRobotPose->setFixed(true);
	optimizer->setVerbose(true);

	cerr << "Optimizing" << endl;
	optimizer->initializeOptimization();
	optimizer->optimize(100);
	cerr << "done." << endl;

	//optimizer->save("tutorial_after.g2o");
}

void GraphManager::getJacobian(){
  VertexSE3 v1;
  v1.setId(0); 

  VertexSE3 v2;
  v2.setId(1); 

  EdgeSE3 e;
  e.setVertex(0, &v1);
  e.setVertex(1, &v2);
  e.setInformation(Eigen::Matrix<double,6,6>::Identity());

  JacobianWorkspace jacobianWorkspace;
  JacobianWorkspace numericJacobianWorkspace;
  numericJacobianWorkspace.updateSize(&e);
  numericJacobianWorkspace.allocate();
}



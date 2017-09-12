/*
 *  Copyright 2013 CNRS
 *
 *  Nirmal Giftsun
 */



#include "sot_collision.hh"
#include <sot/core/debug.hh>

namespace ml = maal::boost;

using namespace std;
using namespace boost::assign;
using namespace boost::tuples;
using namespace fcl;
using namespace ml;
using namespace dynamicgraph::sotcollision;
using namespace dynamicgraph::sot;
using namespace ::dynamicgraph::command;
using namespace dynamicgraph;


DYNAMICGRAPH_FACTORY_ENTITY_PLUGIN(SotCollision, "SotCollision");

SotCollision::SotCollision(const std::string& inName) :
  Entity(inName),
proximitySensorDistanceSIN(NULL,"SotCollision("+name+")::input(vector)::proximitySensorDistance"),
proximitySensorPoseSIN(NULL,"SotCollision("+name+")::input(vector)::proximitySensorPose"),
jointJacobianSIN(NULL,"SotCollision("+name+")::input(matrix)::jointJacobian"),
jointTransformationSIN(NULL,"SotCollision("+name+")::input(matrixHomo)::jointTransformation"),
collisionDistance( boost::bind(&SotCollision::computeimdVector,this,_1,_2),sotNOSIGNAL,"sotCollision("+name+")::output(vector)::"+std::string("collisionDistance")),
collisionJacobian(boost::bind(&SotCollision::computeimdJacobian,this,_1,_2),sotNOSIGNAL,"sotCollision("+name+")::output(matrix)::"+std::string("collisionJacobian")),
closestPoints( boost::bind(&SotCollision::computeClosestPoints,this,_1,_2),collisionJacobian,"sotCollision("+name+")::output(vector)::"+std::string("closestPoints")),
celltr( boost::bind(&SotCollision::computecelltr,this,_1,_2),collisionJacobian,"sotCollision("+name+")::output(vector)::"+std::string("celltr"))
{
			signalRegistration(proximitySensorDistanceSIN);
      signalRegistration(proximitySensorPoseSIN);
			signalRegistration(jointTransformationSIN);
			signalRegistration(jointJacobianSIN);            
			signalRegistration(collisionJacobian);
			signalRegistration(collisionDistance);
      signalRegistration(closestPoints);    
     signalRegistration(celltr);  

      collisionDistance.addDependency(proximitySensorDistanceSIN);
      collisionJacobian.addDependency(proximitySensorDistanceSIN);
      collisionJacobian.addDependency(jointTransformationSIN);
      collisionJacobian.addDependency(jointJacobianSIN);
      closestPoints.addDependency(collisionJacobian);
      celltr.addDependency(collisionJacobian);

      num_skinsensors =0;
      ln.resize(1, 3);

      UNIT_ROTATION(0,0) = 1; UNIT_ROTATION(0,1) = 0; UNIT_ROTATION(0,2)=0;
      UNIT_ROTATION(1,0) = 0; UNIT_ROTATION(1,1) = 1; UNIT_ROTATION(1,2)=0;
      UNIT_ROTATION(2,0) = 0; UNIT_ROTATION(2,1) = 0; UNIT_ROTATION(2,2)=1;
      UNIT_ROTATION(0,3) = 0;  UNIT_ROTATION(1,3) =0 ; UNIT_ROTATION(2,3)=0;

      using namespace ::dynamicgraph::command;
      std::string docstring;
      docstring =
          "\n"
          "    Set dimension\n"
          "      Input: float\n";
      addCommand("setNumSkinSensors",
                 makeDirectSetter(*this, &num_skinsensors, docstring));
}
 

SotCollision::~SotCollision()
{
}


dg::Matrix& SotCollision::computeClosestPoints(dg::Matrix& res, int time)
{
   res  = closestPointsAcc;
   return res;
}

MatrixHomogeneous& SotCollision::computecelltr(MatrixHomogeneous& res, int time)
{
   res  = ctr;
   return res;
}

dg::Vector& SotCollision::computeimdVector(dg::Vector & res, int time)
{
  dg::Vector distance = proximitySensorDistanceSIN(time);
  res.resize( num_skinsensors);
  ln.resize(num_skinsensors,3);
  for (int p = 0; p < distance.size() ; p++)
  {
    res(p) = distance(p);
    ln(p, 0) = p;ln(p, 1) = 0;
    ln(p, 2) = distance(p);
  }
  return res;
}
dg::Matrix &SotCollision::computeimdJacobian(Matrix & res, int time)
{
  closestPointsAcc.resize(6,num_skinsensors); 
  dg::Matrix cellPose = proximitySensorPoseSIN(time);
  //std::cout <<"cp" <<std::endl<<cellPose<<std::endl; 
  dg::Matrix jointJacobian = jointJacobianSIN(time);
  MatrixHomogeneous jointTransformation = jointTransformationSIN(time);
  const int cols = 12;//jointJacobian.cols();
  res.resize(num_skinsensors, cols);  
  dg::Matrix l;l.resize(1,3);
  //Vector position, quat;

  // Eigen::Vector3d position(cellPose(0,0),cellPose(1,0),cellPose(2,0));
  // MatrixHomogeneous cPtrGlobal; 
  // cPtrGlobal.setIdentity();
  // cPtrGlobal.translate(position);
  // Eigen::Quaterniond q(cellPose(6,0), cellPose(3,0), cellPose(4,0), cellPose(5,0)); 
  // cPtrGlobal.rotate(q);
  // //cPtrGlobal.rotate(Eigen::AngleAxisd(cellPose(3,0), VectorRotation::UnitX()) * Eigen::AngleAxisd(cellPose(4,0), VectorRotation::UnitY()) * Eigen::AngleAxisd(cellPose(5,0), VectorRotation::UnitZ()));
 //  ctr = cPtrGlobal;
    // linear velocity and angular velocity
  dg::Matrix Jiv(3, cols), Jiw(3, cols);
  Jiv = jointJacobian.block(0, 0, 3, cols);
  Jiw = jointJacobian.block(3, 0, 3, cols);
  //std::cout <<"jj" <<std::endl<<jointJacobian<<std::endl<<"jiv" <<std::endl<<Jiv<<std::endl<<"jiw" <<std::endl<<Jiw<<std::endl; 
  for (int r = 0; r < num_skinsensors; r++)
  {
    l= ln.block(r,0,1,3);
    //std::cout <<"normal = " <<l<<std::endl;
    Eigen::Vector3d position(cellPose(0,r),cellPose(1,r),cellPose(2,r));
    MatrixHomogeneous cPtrGlobal, cPtrLocal,cRange,cPtrLocalj; 
    cPtrGlobal.setIdentity();
    cPtrGlobal.translate(position);	
      Eigen::Quaterniond q(cellPose(6,r), cellPose(3,r), cellPose(4,r), cellPose(5,r)); 
      cPtrGlobal.rotate(q);
   //std::cout <<"pos" <<std::endl<<cPtrGlobal<<std::endl; 

    //cPtrLocal = jointTransformation.inverse() * cPtrGlobal;   
    cPtrLocal = cPtrGlobal;  
    ctr = cPtrGlobal;
    
    Eigen::Vector3d trz(0,0,l(2));
    cRange.setIdentity();
    cRange.translate(trz);
  
    cPtrLocalj = cPtrLocal * cRange;
    //std::cout <<"crl" <<std::endl<<cPtrLocalj ;
    Eigen::Vector3d n;
    n = cPtrLocalj.translation()-cPtrLocal.translation();
    //std::cout <<"sensor" <<std::endl<<cPtrLocalj.translation()<<std::endl<<"cell" <<cPtrLocal.translation()<<std::endl; 
    //std::cout <<"normal" <<std::endl<<(cPtrLocalj.translation()-cPtrLocal.translation())<<std::endl; 
    //std::cout <<"time= " <<time<<std::endl;
    closestPointsAcc(0,r) = cPtrLocal.translation()(0);closestPointsAcc(1,r) = cPtrLocal.translation()(1);closestPointsAcc(2,r) = cPtrLocal.translation()(2);
    closestPointsAcc(3,r) = cPtrLocalj.translation()(0);closestPointsAcc(4,r) = cPtrLocalj.translation()(1);closestPointsAcc(5,r) = cPtrLocalj.translation()(2);
    //closestPointsAcc(0,r) = 0;closestPointsAcc(1,r) = 1;closestPointsAcc(2,r) = 2;
    //closestPointsAcc(3,r) = 3;closestPointsAcc(4,r) = 4;closestPointsAcc(5,r) = 5;

    dg::Matrix cellPointSkewMatrix;
    cellPointSkewMatrix.resize(3, 3);
    cellPointSkewMatrix(0, 0) = 0;
    cellPointSkewMatrix(0, 1) = -cPtrLocal(2, 3);
    cellPointSkewMatrix(0, 2) = cPtrLocal(1, 3);
    cellPointSkewMatrix(1, 0) = cPtrLocal(2, 3);
    cellPointSkewMatrix(1, 1) = 0;
    cellPointSkewMatrix(1, 2) = -cPtrLocal(0, 3);
    cellPointSkewMatrix(2, 0) = -cPtrLocal(1, 3);
    cellPointSkewMatrix(2, 1) = cPtrLocal(0, 3);
    cellPointSkewMatrix(2, 2) = 0;
    
    dg::Matrix cellPointGradient;
    cellPointGradient = Jiv - (cellPointSkewMatrix * Jiw); 
    dg::Matrix cj;  
    n.normalize(); 
    cj = -n.transpose() * cellPointGradient;  
    for (int d = 0; d < cols; d++)
    {
      res(r, d) = cj(0, d);
    }
  }

  return res;
}


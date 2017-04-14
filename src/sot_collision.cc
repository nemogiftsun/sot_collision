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
proximitySensorSIN(NULL,"SotCollision("+name+")::input(vector)::proximitySensor"),
fclmodelupdateSINTERN(boost::bind(&SotCollision::updatefclmodels,this,_1,_2),sotNOSIGNAL,"sotCollision("+name+")::intern(dummy)::fclmodelupdateSINTERN" ),
collisionDistance( boost::bind(&SotCollision::computeimdVector,this,_1,_2),
fclmodelupdateSINTERN,
"sotCollision("+name+")::output(vector)::"+std::string("collisionDistance")),
collisionJacobian(boost::bind(&SotCollision::computeimdJacobian,this,_1,_2),
fclmodelupdateSINTERN,
"sotCollision("+name+")::output(matrix)::"+std::string("collisionJacobian")),
collisionModelState( boost::bind(&SotCollision::computeCollisionModelState,this,_1,_2),
fclmodelupdateSINTERN,
"sotCollision("+name+")::output(vector)::"+std::string("collisionModelState")),
closestPointi( boost::bind(&SotCollision::computeClosestPointi,this,_1,_2),
fclmodelupdateSINTERN,
"sotCollision("+name+")::output(vector)::"+std::string("closestPointi")),
closestPointj( boost::bind(&SotCollision::computeClosestPointj,this,_1,_2),
fclmodelupdateSINTERN,
"sotCollision("+name+")::output(vector)::"+std::string("closestPointj"))
{
			signalRegistration(proximitySensorSIN);
			signalRegistration(collisionJacobian);
			signalRegistration(collisionDistance);
      signalRegistration(collisionModelState);
      signalRegistration(closestPointj);
      signalRegistration(closestPointi);	

      collisionDistance.addDependency(proximitySensorSIN);
      collisionDistance.addDependency(fclmodelupdateSINTERN);
      collisionJacobian.addDependency(fclmodelupdateSINTERN);
      collisionModelState.addDependency(fclmodelupdateSINTERN);
      closestPointi.addDependency(fclmodelupdateSINTERN);
      closestPointj.addDependency(fclmodelupdateSINTERN);

      dimension = 0;
      num_collisionpairs = 0;
 
      UNIT_ROTATION(0,0) = 1; UNIT_ROTATION(0,1) = 0; UNIT_ROTATION(0,2)=1;
      UNIT_ROTATION(1,0) = 0; UNIT_ROTATION(1,1) = 1; UNIT_ROTATION(1,2)=0;
      UNIT_ROTATION(2,0) = 0; UNIT_ROTATION(2,1) = 0; UNIT_ROTATION(2,2)=1;
      UNIT_ROTATION(0,3) = 0;  UNIT_ROTATION(1,3) =0 ; UNIT_ROTATION(2,3)=0;
 
      using namespace ::dynamicgraph::command;
      std::string docstring;
      docstring =
        "\n"
        "    Get Collision Checks\n"
        "\n";
      addCommand(std::string("getcheck"),
	         new ::dynamicgraph::command::Getter<SotCollision, bool>
	         (*this, &SotCollision::getcheck, docstring));
      docstring =
        "\n"
        "    Get dimension of the collision model\n"
        "\n";
      addCommand(std::string("getdimension"),
	         new ::dynamicgraph::command::Getter<SotCollision, double>
	         (*this, &SotCollision::getdimension, docstring));	 
      

      addCommand(std::string("getfclstate"),
	         new ::dynamicgraph::command::Getter<SotCollision, dg::Matrix>
	         (*this, &SotCollision::getfclstate, docstring));	 	
  
 
////// new important and useful commands

      addCommand(std::string("createcollisionlink"),
	         makeCommandVoid4(*this,&SotCollision::createcollisionlink,docstring));

      addCommand(std::string("createcollisionpair"),
	         makeCommandVoid2(*this,&SotCollision::createcollisionpair,docstring));

}
 

SotCollision::~SotCollision()
{
}

int& SotCollision::updatefclmodels(int& dummy,int time)
{



// update the fcl models from current transformation matrix
    for(int i=0;i<dimension;++i)
    {       
       fclstate.resize(dimension,10);

       MatrixHomogeneous T = body_transformation_input[i]->access(time);

       
       Matrix3f rot(T(0,0),T(0,1),T(0,2),T(1,0),T(1,1),T(1,2),T(2,0),T(2,1),T(2,2));
       Vec3f trans(T(0,3),T(1,3),T(2,3));
       
       transform_links[i]= Transform3f(rot,trans);
       transform_links[i] = transform_links[i]* transform_joint_links[i];
       
       Vec3f tpart = transform_links[i].getTranslation();
       Quaternion3f quat = transform_links[i].getQuatRotation();

       fclstate(i,0) = tpart[0];
       fclstate(i,1) = tpart[1]; 
       fclstate(i,2) = tpart[2];
      
       fclstate(i,3) = quat.getX();  
       fclstate(i,4) = quat.getY();
       fclstate(i,5) = quat.getZ();
       fclstate(i,6) = quat.getW(); 

       if(link_info_map[i].link_type == "capsule")
       {

               fclstate(i,7) = capsule_links[link_info_map[i].index].radius; 
               fclstate(i,8) = capsule_links[link_info_map[i].index].lz;  
               fclstate(i,9) = 0;  
       }
       else if(link_info_map[i].link_type == "box")
       {
               fclstate(i,7) = box_links[link_info_map[i].index].side[0]; 
               fclstate(i,8) = box_links[link_info_map[i].index].side[1];   
               fclstate(i,9) = box_links[link_info_map[i].index].side[2];     
       }
   
    }

// calculate inter model distance 
// ONLY INTER CAPSULE DISTANCE COMPUTED. A GENERALIZATION WILL BE INTRODUCED LATER 

    // define solver
    GJKSolver_libccd solver;
    GJKSolver_indep solver1;

    
    // variables for fcl models
    FCL_REAL dist;
    Vec3f l1 , l2;

    // compute the imdm
    //  This matrix is a vector of dimension * dimension 
    // (distance,closestpointinthecorrespondingrowlink)

    for(int i=0;i<dimension;++i)
    {

				//Capsule linka (capsule_links[i].radius, capsule_links[i].lz);

				Transform3f linka_transform  = transform_links[i];

				imdm.at(i).at(i) = make_tuple(0,(0,0,0));
       
				for(int j=i+1;j<dimension;++j)
				{                

				    if(link_info_map[i].link_type == "capsule")
				    {

								if(link_info_map[j].link_type == "capsule")
								{
                	  check = solver1.shapeDistance(capsule_links[link_info_map[i].index],transform_links[i], capsule_links[link_info_map[j].index],transform_links[j], &dist, &l1, &l2);

								}
                else if(link_info_map[j].link_type == "box")
								{
                	  check = solver1.shapeDistance(capsule_links[link_info_map[i].index],transform_links[i], box_links[link_info_map[j].index],transform_links[j], &dist, &l1, &l2);
								}

				    }
				    else if(link_info_map[i].link_type == "box")
				    {

								if(link_info_map[j].link_type == "capsule")
								{
                	  check = solver1.shapeDistance(box_links[link_info_map[i].index],transform_links[i], capsule_links[link_info_map[j].index],transform_links[j], &dist, &l1, &l2);

								}
                else if(link_info_map[j].link_type == "box")
								{
                         
                         //normalized distance vector
                         //normalized distance vector
                     if (link_info_map[j].link_location == "external"){
                            dg::Vector distance = proximitySensorSIN(time);
                            dist = (distance(i));
                            //l1 = (0,0,0.0);
                            //l2 = (0,0,((dist/20))); 
			    Matrix3f rotation;
			    //l1 = l1 + transform_links[i].transform(l1);
			    //l2 = l2 + transform_links[i].transform(l2);
                            Vec3f lchange(0,0,0);
			    rotation.setEulerZYX(0,0,0);
                            Transform3f tli=transform_links[i]; Transform3f tlj=transform_links[i];
                            Transform3f temp(rotation,lchange);
                            l1 = (tli*temp).getTranslation();
			    //Vec3f linit(0,0,0);
                            Vec3f lchange1(0,0,((dist)));
                            //lchange = (0,0,-0.1);
                            //temp.setTranslation(lchange);   
                            Transform3f temp1(rotation,lchange1); 
                            //Vec3f lchange = (0,0,-dist);
                            //transform_links[i].setTranslation(lchange);                  
                            //l2 = (tlj*temp1).getTranslation(); 
                            l2 = lchange1;
                            //l1 = linit;
                            
                            //std::cout << "transform x" << transform_links[i].getTranslation()[0] <<" y "<<transform_links[i].getTranslation()[1]<<" z "<<transform_links[i].getTranslation()[2]<<std::endl;    
			    //std::cout << "transform rotation x " << transform_links[i].getQuatRotation().getX() << " rotation y " << transform_links[i].getQuatRotation().getY() << "transform rotation z " << transform_links[i].getQuatRotation().getZ()<< "transform rotation w " << transform_links[i].getQuatRotation().getW() << std::endl;    
   			    //std::cout << "l1 = "  <<l1 <<std::endl;
     			    //std::cout << "l2 = "  <<l2 <<std::endl;


                     }
                     else{
                	  check = solver1.shapeDistance(box_links[link_info_map[i].index],transform_links[i], box_links[link_info_map[j].index],transform_links[j], &dist, &l1, &l2);

                           l1 = transform_links[i].transform(l1);
                           
                           l2 = transform_links[j].transform(l2); }                                        
								}

				    }

                        //std::cout << dist<<std::endl;
		        imdm.at(i).at(j) = make_tuple(dist,make_tuple(l1[0],l1[1],l1[2]));
		        imdm.at(j).at(i) = make_tuple(dist,make_tuple(l2[0],l2[1],l2[2]));     
		    }      
     }
   
//compute inter model jacobian


    return dummy;
}


dg::Matrix& SotCollision::computeimdJacobian(Matrix& res, int time)
{

    closestPointI.resize(3,num_collisionpairs);
    closestPointJ.resize(3,num_collisionpairs);


for (int p=0; p < num_collisionpairs;p++)
{

  std::string body0 = collision_pairs[p][0];
  //std::cout << body0 <<std::endl;
  std::string body1 = collision_pairs[p][1];
  //std::cout << body1 <<std::endl;
  // update fcl models
  fclmodelupdateSINTERN(time);
  
  // init the right value from map
  int i = fcl_body_map[body0];
  int j = fcl_body_map[body1];

  //normalized line segment
  dg::Matrix ln,la;
  ln.resize(1,3); 
  la.resize(1,3);


  // position input
  MatrixHomogeneous Mi = body_transformation_input[i]->access(time);
  MatrixHomogeneous Mj = body_transformation_input[j]->access(time);
  dg::Matrix Ji = body_jacobian_input[i]->access(time);
  dg::Matrix Jj = body_jacobian_input[j]->access(time);

  const int cols = Ji.cols();
  // resizing result
  res.resize(num_collisionpairs,cols);
  //std::cout << res <<std::endl;

  // compute collision jacobian distance
     // l1 and l2 vectors
  Vec3f l1(imdm[i][j].get<1>().get<0>(),imdm[i][j].get<1>().get<1>(),imdm[i][j].get<1>().get<2>());
  Vec3f l2(imdm[j][i].get<1>().get<0>(),imdm[j][i].get<1>().get<1>(),imdm[j][i].get<1>().get<2>());
  Vec3f l = l1-l2;
  l.normalize();

     closestPointI(0,p) = l1[0];closestPointI(1,p) = l1[1];closestPointI(2,p) = l1[2];
     closestPointJ(0,p) = l2[0];closestPointJ(1,p) = l2[1];closestPointJ(2,p) = l2[2];
     //std::cout << "cjl1 = "  <<l1 <<std::endl;
     //std::cout << "cjl2 = "  <<l2 <<std::endl;
     ln(0,0)= l2[0]; ln(0,1)=l2[1] ; ln(0,2)=l2[2];
     la(0,0)= l1[0]; la(0,1)=l1[1] ; la(0,2)=l1[2];
     // closest point computation in object a and b

     MatrixHomogeneous aclosestpointtransform, bclosestpointtransform = UNIT_ROTATION;
      //l1(0,0)= 0; ln(0,1)=0.045 ; ln(0,2)=0;
     aclosestpointtransform(0,3) = l1[0];  aclosestpointtransform(1,3) =l1[1] ; aclosestpointtransform(2,3)=l1[2];
     //aclosestpointtransform(0,3) = 0;  aclosestpointtransform(1,3) =0.045 ; aclosestpointtransform(2,3)=0;
     bclosestpointtransform(0,3) = l2[0];  bclosestpointtransform(1,3) =l2[1] ; bclosestpointtransform(2,3)=l2[2];

     MatrixHomogeneous aclosestpointlocalframe, bclosestpointlocalframe;

     // find closest point transform in object a with respect the local joint frame
     aclosestpointlocalframe = Mi.inverse()*aclosestpointtransform;
     //std::cout << "localframea = "  <<aclosestpointlocalframe(0,3)<<","<<aclosestpointlocalframe(1,3)<<","<<aclosestpointlocalframe(2,3) <<std::endl;
     bclosestpointlocalframe = Mj.inverse()*bclosestpointtransform;

     // construct antisymmetric matrices of the closest point for both a and b objects
     dg::Matrix ssmaclosestpoint, ssmbclosestpoint;
     ssmaclosestpoint.resize(3,3);
     ssmbclosestpoint.resize(3,3);

     ssmaclosestpoint(0,0) = 0 ; ssmaclosestpoint(0,1) = -aclosestpointlocalframe(2,3); ssmaclosestpoint(0,2) = aclosestpointlocalframe(1,3) ;
     ssmaclosestpoint(1,0) = aclosestpointlocalframe(2,3) ; ssmaclosestpoint(1,1) = 0; ssmaclosestpoint(1,2) = -aclosestpointlocalframe(0,3) ;
     ssmaclosestpoint(2,0) = -aclosestpointlocalframe(1,3); ssmaclosestpoint(2,1) = aclosestpointlocalframe(0,3); ssmaclosestpoint(2,2) = 0;
    /*
    ssmaclosestpoint(0,0) = aclosestpointlocalframe(0,3) ; ssmaclosestpoint(0,1) = 0; ssmaclosestpoint(0,2) = 0 ;
     ssmaclosestpoint(1,0) = 0 ; ssmaclosestpoint(1,1) = aclosestpointlocalframe(1,3); ssmaclosestpoint(1,2) = 0 ;
     ssmaclosestpoint(2,0) = 0; ssmaclosestpoint(2,1) = 0; ssmaclosestpoint(2,2) = aclosestpointlocalframe(2,3);*/


     ssmbclosestpoint(0,0) = 0 ; ssmbclosestpoint(0,1) = -bclosestpointlocalframe(2,3); ssmbclosestpoint(0,2) = bclosestpointlocalframe(1,3) ;
     ssmbclosestpoint(1,0) = bclosestpointlocalframe(2,3) ; ssmbclosestpoint(1,1) = 0; ssmbclosestpoint(1,2) = -bclosestpointlocalframe(0,3) ;
     ssmbclosestpoint(2,0) = -bclosestpointlocalframe(1,3); ssmbclosestpoint(2,1) = bclosestpointlocalframe(0,3); ssmbclosestpoint(2,2) = 0;
     //std::cout << "ssma = "  <<ssmaclosestpoint;

     // linear velocity and angular velocity for both joints 
     dg::Matrix Jiv(3,cols),Jiw(3,cols),Jjv(3,cols),Jjw(3,cols);


     //Ji.extract(0,0,3,cols,Jiv);
     //Ji.extract(3,0,3,cols,Jiw);
     Jiv = Ji.block(0,0,3,cols);
     Jiw = Ji.block(2,0,3,cols);

     //std::cout << "Ji= "  <<Ji<<std::endl;
     //std::cout << "jiv = "  <<Jiv<<std::endl;
     //std::cout << "jiw = "  <<Jiw<<std::endl;

     if(link_info_map[j].link_location == "internal")
     {
	Jjv = Jj.block(0,0,3,cols);
        Jjw = Jj.block(2,0,3,cols);
        //Jj.extract(0,0,3,cols,Jjv);
        //Jj.extract(2,0,3,cols,Jjw);
     }

     // point gradient computation
     dg::Matrix abpointgradient;

     if(link_info_map[j].link_location == "internal")
     {
     abpointgradient  = (Jiv-(ssmaclosestpoint*Jiw)) - (Jjv-(ssmbclosestpoint*Jjw));
     }
     else
     {
      abpointgradient  = Jiv - (ssmaclosestpoint*Jiw);    
      //abpointgradient  = (Jiv+(la*Jiw)); 
     }
     //compute collision jacobian
     dg::Matrix cj;
     //ln 3*1
     //abpoint 3*3
     sotDEBUG(1)<<"normal = " <<ln<<endl;
     //std::cout <<"normal = " <<ln<<std::endl;
     
     cj = -ln*abpointgradient;
     for (int d=0;d<cols;d++)
     {
      res(p,d) = cj(0,d); 
     }
     
   } 
   return res;
}

dg::Matrix& SotCollision::computeCollisionModelState(dg::Matrix& res, int time)
{
   fclmodelupdateSINTERN(time);
   res  = fclstate;
   //std::cout << "fclstate x" << fclstate(0,0) <<" y "<<fclstate(0,1)<<" z"<<fclstate(0,2)<<std::endl;
   //std::cout << "fclstate rot x" << fclstate(0,3) <<" rot y "<<fclstate(0,4)<<" rot z"<<fclstate(0,5)<<" rot w"<<fclstate(0,6)<<std::endl;
   return res;
}

dg::Matrix& SotCollision::computeClosestPointi(dg::Matrix& res, int time)
{
   fclmodelupdateSINTERN(time);
   res  = closestPointI;
   return res;
}

dg::Matrix& SotCollision::computeClosestPointj(dg::Matrix& res, int time)
{
   fclmodelupdateSINTERN(time);
   res  = closestPointJ;
   return res;
}


dg::Vector& SotCollision::computeimdVector(dg::Vector& res, int time)
{

  fclmodelupdateSINTERN(time);
  res.resize(num_collisionpairs);
for (int p=0; p < num_collisionpairs;p++)
{
  std::string body0 = collision_pairs[p][0];
  std::string body1 = collision_pairs[p][1];
  int i = fcl_body_map[body0];
  int j = fcl_body_map[body1];

	res(p)= imdm[i][j].get<0>();
}
  //std::cout << "The imdvector "<<i<<"is "<<res;

  return res;
}


dynamicgraph::SignalPtr< MatrixHomogeneous,int >& SotCollision::createPositionSignalIN( const std::string& signame)
{


  dynamicgraph::SignalPtr< MatrixHomogeneous,int> * sig
    = new dynamicgraph::SignalPtr< MatrixHomogeneous,int>
    (NULL,"sotCollision("+name+")::input(matrix)::"+signame);

  genericSignalRefs.push_back( sig );
  signalRegistration( *sig );

 

  return *sig;  
}

dynamicgraph::SignalPtr< dg::Matrix,int >& SotCollision::createJacobiansignalIN( const std::string& signame)
{


  dynamicgraph::SignalPtr< dg::Matrix,int> * sig
    = new dynamicgraph::SignalPtr< dg::Matrix,int>
    (NULL,"sotCollision("+name+")::input(matrix)::"+signame);

  genericSignalRefs.push_back( sig );
  signalRegistration( *sig );


  return *sig;  
}



void SotCollision::createcollisionpair(const std::string& body0,const std::string& body1)
{

    num_collisionpairs = num_collisionpairs+1;
    collision_pairs.resize(num_collisionpairs);
    for(int i = 0; i<num_collisionpairs;i++)
    {  
     collision_pairs[i].resize(2); 
    }
    collision_pairs[num_collisionpairs-1][0] = body0;
    collision_pairs[num_collisionpairs-1][1] = body1;

    // init the right value from map
    int i = fcl_body_map[body0];
    int j = fcl_body_map[body1];

    collisionDistance.addDependency( *body_transformation_input[i] );
    collisionDistance.addDependency( *body_transformation_input[j] );

    collisionJacobian.addDependency( *body_jacobian_input[i] );
    collisionJacobian.addDependency( *body_jacobian_input[j] );

}


void SotCollision::createcollisionlink(const std::string& bodyname, const std::string& bodytype,const std::string& bodynaturelocation,const dg::Vector& bodydescription)
{

// set the dimension for fcl checks
    dimension  = dimension+1;
    // resize the model parameters and the transformations 
    capsule_links.reserve(dimension);
    transform_links.resize(dimension);
    // joint-link transformation matrix to be in sync with rviz as orientation representation of ros is different with sot-dynamics
    transform_joint_links.resize(dimension);
    // resize intermodel distance matrix
		imdm.resize(dimension);
		for(int i=0;i<dimension;++i)
		{
		    imdm[i].resize(dimension);
		}  

    // create input signals
    //    - pose signal from sot dynamic
    // create output signals    
    //    - inter fcl model distance signal for each body 
    //    - jacobian signal for each body

    Matrix3f rotation;

      // create input position/jacbian signals which is supposedly to be plugged to sot dynamic position/jacobian signals
        int i = dimension - 1;
        std::string signame_in = bodyname;
        fcl_body_map.insert(std::pair<std::string,int>(signame_in,i) ); 
        body_transformation_input[i] = &createPositionSignalIN(signame_in); 
        body_jacobian_input[i] = &createJacobiansignalIN(std::string("J")+signame_in); 

        //ShapeBase collision_link;
        if(bodytype == "capsule")
        {
		      Capsule capsule (bodydescription(0),bodydescription(1));
		      capsule_links.push_back(capsule);
          link_info link;
          link.link_type = "capsule";
          link.index = capsule_links.size()-1;
          link.link_location = bodynaturelocation;
          link_info_map.insert(std::pair<int,link_info>(i,link)); 
          
        }
        else if(bodytype == "box")
        {
          Box box(bodydescription(0),bodydescription(1),bodydescription(2));
          box_links.push_back(box);
          link_info link;
          link.link_type = "box";
          link.index = box_links.size()-1;
          link.link_location = bodynaturelocation;
          link_info_map.insert(std::pair<int,link_info>(i,link)); 
        }
        
        //create transforms
        Vec3f translation(bodydescription(3),bodydescription(4),bodydescription(5));
        rotation.setEulerZYX(bodydescription(8),bodydescription(7),bodydescription(6));
        Transform3f joint_link_transform(rotation,translation);

        // push joint link transforms extracted from body description
        transform_joint_links[i] = joint_link_transform;             
        // initialize the fcl body transforms
        //transform_links.push_back(joint_link_transform);

      
      fclmodelupdateSINTERN.addDependency( *body_transformation_input[i] );

      fclmodelupdateSINTERN.addDependency(*body_jacobian_input[i]);

}


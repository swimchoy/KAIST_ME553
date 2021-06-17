//
// Created by suyoung on 21. 6. 17..
//


// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include "raisim/RaisimServer.hpp"
#include "exercise_9_20204577.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.setGravity({0,0,-9.81}); //// use this gravity in your ABA!!!!

  // kinova
  auto gantry = world.addArticulatedSystem("/home/suyoung/Workspace/Lectures/ME553/rsc/twolink/urdf/twolink.urdf");
  gantry->setName("gantry");
  server.focusOn(gantry);

  // kinova configuration
  Eigen::VectorXd gc(gantry->getGeneralizedCoordinateDim()), gv(gantry->getDOF()), gf(gantry->getDOF());
  gc << 1.0, M_PI_2+0.5+0.5;
  gv << 0.1, 0.2;
  gf << 0.1, 0.2;

//  gc << 0.0, 0.1, 0.2, 0.0, 1.0, 0.0, 0.0, 0.0, 0.2;
//  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.5, 0.4;
//  gf << 0.1, 0.2, 0.3, 0.4, 0.4, 0.3, 0.2, 0.3;

  gantry->setState(gc, gv);
  gantry->setGeneralizedForce(gf);
  world.integrate1();

  std::cout<<gantry->getMassMatrix().e()<<std::endl;
  std::cout<<gantry->getNonlinearities().e()<<std::endl;

//  Eigen::MatrixXd J;
//  J.setZero(3, 3);
//  gantry->getDenseFrameJacobian("ball_fixed", J);
//  std::cout<<J<<std::endl;
//
//  std::cout<<(J * gantry->getMassMatrix().e().inverse() * J.transpose()).inverse() << std::endl;

//  if((getMassMatrixUsingCRBA(gc, kinova) - kinova->getMassMatrix().e()).norm() < 1e-8)
//    std::cout<<"CRBA passed"<<std::endl;
//  else
//    std::cout<<"CRBA failed"<<std::endl;
//
//  if((getNonlinearitiesUsingRNE(gc, gv, kinova) - kinova->getNonlinearities().e()).norm() < 1e-8)
//    std::cout<<"RNE passed "<<std::endl;
//  else
//    std::cout<<"RNE failed"<<std::endl;
//
//  Eigen::MatrixXd Minv = kinova->getInverseMassMatrix().e();
//  Eigen::VectorXd b = kinova->getNonlinearities().e();
//  Eigen::VectorXd acc_raisim = Minv * (gf-b);
//
//  if ((getGaUsingABA(gc, gv, gf) - acc_raisim).norm() < 1e-8)
//    std::cout<<"ABA passed"<<std::endl;
//  else
//    std::cout<<"ABA failed"<<std::endl;
}

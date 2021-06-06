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
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/kinova/urdf/kinova_no_fingers.urdf");
  kinova->setName("kinova");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd gc(kinova->getGeneralizedCoordinateDim()), gv(kinova->getDOF()), gf(kinova->getDOF());
  gc << 0.0, 0.1, 0.2, 0.0, 2.0, 0.2;
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  gf << 0.1, 0.2, 0.3, 0.4, 0.4, 0.3;

  kinova->setState(gc, gv);
  kinova->setGeneralizedForce(gf);
  world.integrate1();

  Eigen::MatrixXd Minv = kinova->getInverseMassMatrix().e();
  Eigen::VectorXd b = kinova->getNonlinearities().e();
  Eigen::VectorXd acc_raisim = Minv * (gf-b);

  if ((getGaUsingABA(gc, gv, gf) - acc_raisim).norm() < 1e-8)
    std::cout<<"ABA passed"<<std::endl;
  else
    std::cout<<"ABA failed"<<std::endl;
}

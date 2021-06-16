#include "raisim/RaisimServer.hpp"
#include "exercise_8_20204577.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();

  // anymal
  auto anymal = world.addArticulatedSystem(binaryPath.getDirectory() + "../rsc/anymal/urdf/anymal_oneleg.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);

  // anymal configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(anymal->getDOF());
  Eigen::VectorXd gf(anymal->getDOF());

//  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8;
//  gv << 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1;
//  gf << 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1;

  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, 1.0, 0., 0., 0.;
  gv << 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1, 0.2, 0.3;
  gf << 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.1, 0.1, 0.2, 0.3;

  anymal->setState(gc, gv);
  world.integrate1();

  if((getMassMatrixUsingCRBA(gc) - anymal->getMassMatrix().e()).norm() < 1e-8)
    std::cout<<"CRBA passed"<<std::endl;
  else
    std::cout<<"CRBA failed"<<std::endl;

  if((getNonlinearitiesUsingRNE(gc, gv) - anymal->getNonlinearities().e()).norm() < 1e-8)
    std::cout<<"RNE passed "<<std::endl;
  else
    std::cout<<"RNE failed"<<std::endl;

  Eigen::MatrixXd Minv = anymal->getInverseMassMatrix().e();
  Eigen::VectorXd b = anymal->getNonlinearities().e();
  Eigen::VectorXd acc_raisim = Minv * (gf-b);

  if ((getGaUsingABA(gc, gv, gf) - acc_raisim).norm() < 1e-8)
    std::cout<<"ABA passed"<<std::endl;
  else
    std::cout<<"ABA failed"<<std::endl;
  

  return 0;
}

#include "raisim/RaisimServer.hpp"
#include "exercise_6_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.0001);

  // kinova
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/kinova/urdf/kinova_no_fingers_single_collision.urdf");
  kinova->setName("kinova");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd gc(kinova->getGeneralizedCoordinateDim()), gv(kinova->getDOF());
  gc << 0.0, 2.56, -1.5, 0.0, 2.0, 0.0;
  gv.setZero();
  kinova->setState(gc, gv);

  /// this function updates internal variables in raisim (such as the mass matrix and nonlinearities)
  world.integrate1();

  // visualization
  server.launchServer();

  /// visualize to debug!
  for (int i=0; i<2000000; i++) {
    /// your code here
    kinova->setGeneralizedForce(getGeneralizedForce(kinova));
    world.integrate2();
    if(!world.getContactProblem()->empty()) {
      if(std::pow(world.getContactProblem()->at(0).imp_i.norm()/world.getTimeStep() - 10,2) <1e-5)
        std::cout<<"passed"<<std::endl;
      else
        std::cout<<"failed. pushing at "<<world.getContactProblem()->at(0).imp_i.norm()/world.getTimeStep()<<std::endl;
    }
    /// you get grade if you pass for 10 consecutive runs

    world.integrate1();
  }

  server.killServer();
}

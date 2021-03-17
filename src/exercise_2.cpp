// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include "raisim/RaisimServer.hpp"
#include "exercise_2_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.0001);

  // anymal
  auto anymal = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/anymal/urdf/anymal.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);

  // anymal configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(anymal->getDOF());

  gc << 0, 0, 10.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  gv << 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4;
  anymal->setState(gc, gv);

  // visualization
  server.launchServer();
  raisim::Vec<3> footVel, footAngVel;
  for (int i=0; i<2000000; i++) {
    std::this_thread::sleep_for(std::chrono::microseconds(1000));
    world.integrate1();

    anymal->getFrameVelocity("LF_ADAPTER_TO_FOOT", footVel);
    anymal->getFrameAngularVelocity("LF_ADAPTER_TO_FOOT", footAngVel);
//    std::cout<<"foot velocity "<<footVel<<std::endl;

    if((footVel.e() - getFootLinearVelocity(gc, gv)).norm() < 1e-10) {
      std::cout<<"the linear velocity is correct "<<std::endl;
    } else {
      std::cout<<"the linear velocity is not correct "<<std::endl;
    }

    if((footAngVel.e() - getFootAngularVelocity(gc, gv)).norm() < 1e-10) {
      std::cout<<"the angular velocity is correct "<<std::endl;
    } else {
      std::cout<<"the angular velocity is not correct "<<std::endl;
    }

    world.integrate2();
  }

  server.killServer();
}

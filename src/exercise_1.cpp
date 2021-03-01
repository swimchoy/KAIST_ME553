// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include "raisim/RaisimServer.hpp"

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world;
  raisim::RaisimServer server(&world);
  world.addGround();

  // kinova
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "\\rsc\\kinova\\urdf\\kinova.urdf");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd jointNominalConfig(kinova->getGeneralizedCoordinateDim());
  jointNominalConfig << 0.0, 2.76, -1.57, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0;
  kinova->setGeneralizedCoordinate(jointNominalConfig);
  kinova->setName("kinova");

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.15);
  debugSphere->setColor(1,0,0,1);
  debugSphere->setPosition(1,1,0);

  // compute the position of the end-effector
  raisim::rpyToRotMat_extrinsic()

  // visualization
  server.launchServer();
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}

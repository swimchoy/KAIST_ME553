// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include "raisim/World.hpp"
#include "raisim/RaisimServer.hpp"

int main(int argc, char *argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  /// create raisim world
  raisim::World world;
  raisim::RaisimServer server(&world);
  world.addGround(-2);

  // monkey
  std::string monkeyFile = binaryPath.getDirectory() + "\\rsc\\monkey.obj";
  auto monkey = world.addMesh(monkeyFile, 1.0, raisim::Mat<3, 3>::getIdentity(), {0, 0, 0});
  monkey->setName("monkey");

  // debug sphere
  auto debugSphere = server.addVisualSphere("debug_sphere", 0.15);
  debugSphere->setColor(0, 1, 0, 1);
  Eigen::Vector3d spherePosition{0.0, -0.2, 0.8};
  debugSphere->setPosition(spherePosition);

  // rotation
  raisim::Mat<3, 3> rot_ws;
  raisim::angleAxisToRotMat({1, 0, 0}, M_PI_2, rot_ws);
  Eigen::Matrix3d rot_ws_e = rot_ws.e();

  // rotate the debug sphere
  Eigen::Vector3d newSpherePosition = rot_ws_e * spherePosition;
  debugSphere->setPosition(newSpherePosition);

  // rotate the monkey
  monkey->setOrientation(rot_ws_e);

  // launch raisim server
  server.launchServer();

  for (int i = 0; i < 100000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}

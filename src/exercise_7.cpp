#include "raisim/RaisimServer.hpp"
#include "exercise_7_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  // graphics objects
  auto sphere1 = server.addVisualSphere("sphere", 0.5, 0, 1, 0);
  sphere1->setPosition(0,0,0.5);
  auto sphere2 = server.addVisualSphere("sphere2", 0.7, 0, 1, 1);
  sphere2->setPosition(0.1, 0.1, 3);

  // visualization
  server.launchServer();
  SimulationClass simulation_class;

  /// visualize to debug!
  for (int i=0; i<2000000; i++) {
    simulation_class.integrate();
    simulation_class.setPosition(sphere1, sphere2);
    std::this_thread::sleep_for(std::chrono::microseconds(1000));
  }

  server.killServer();
}

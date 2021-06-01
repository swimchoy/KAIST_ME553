#include "raisim/RaisimServer.hpp"
#include "exercise_7_20204577.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround(0., "ground");
  world.setTimeStep(0.001);

  // graphics objects
  auto sphere1 = server.addVisualSphere("sphere", 0.5, 0, 1, 0);
  sphere1->setPosition(0,0,0.5);
  auto sphere2 = server.addVisualSphere("sphere2", 0.7, 0, 1, 1);
  sphere2->setPosition(0.1, 0.1, 3);

  // visualization
  server.launchServer();
  SimulationClass simulation_class;

  /// to debug
  raisim::Vec<3> g {0, 0, -9.81};
  world.setGravity(g);

  auto debugsphere1 = world.addSphere(0.5, 0.1 * 4/3 * M_PI * std::pow(0.5,3), "ball1");
  auto debugsphere2 = world.addSphere(0.7, 0.1 * 4/3 * M_PI * std::pow(0.7,3), "ball2");

  debugsphere1->setPosition(0, 0, 0.5);
  debugsphere2->setPosition(0.1, 0.1, 3);

  world.setMaterialPairProp("ground", "ball1", 0.8, 0.0, 0.01);
  world.setMaterialPairProp("ground", "ball2", 0.8, 0.0, 0.01);
  world.setMaterialPairProp("ball1", "ball2", 0.8, 0.0, 0.01);


  /// visualize to debug!
  for (int i=0; i<2000000; i++) {

    simulation_class.integrate();
    simulation_class.setPosition(sphere1, sphere2);


    world.integrate();
    std::this_thread::sleep_for(std::chrono::microseconds(1000));
  }

  server.killServer();
}

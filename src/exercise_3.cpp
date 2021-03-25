#include "raisim/RaisimServer.hpp"
#include "exercise_3_STUDENTID.hpp"


int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setGravity({0,0,0});

  // kinova
  auto kinova = world.addArticulatedSystem(binaryPath.getDirectory() + "/rsc/kinova/urdf/kinova_no_fingers.urdf");
  kinova->setName("kinova");
  server.focusOn(kinova);

  // kinova configuration
  Eigen::VectorXd gc(kinova->getGeneralizedCoordinateDim()), gv(kinova->getDOF());
  gc << 0.0, 2.56, -1.5, 0.0, 2.0, 0.0;
  kinova->setGeneralizedCoordinate(gc);
  kinova->setIntegrationScheme(raisim::ArticulatedSystem::IntegrationScheme::EULER);
  raisim::Vec<3> posDes{0.5, -0.1, 0.4}, posCur, posDiff;
  raisim::Vec<4> quatDes{0, 1, 0, 0}, quatCur, quatDiff;
  raisim::Mat<3,3> rotCur;

  // visualization
  server.launchServer();

  for (int i=0; i<2000000; i++) {
    kinova->getState(gc, gv);
    world.integrate1();

    /// DO NOT USE OR MODIFY THIS CODE. /////////////////////////////////////////////
    kinova->getFramePosition("kinova_joint_end_effector",posCur);
    kinova->getFrameOrientation("kinova_joint_end_effector",rotCur);
    raisim::rotMatToQuat(rotCur, quatCur);
    posDiff = posDes - posCur;
    quatDiff = quatDes - quatCur;
    if(posDiff.e().norm() < 1e-5 && quatDiff.e().norm() < 1e-4) {
      std::cout<<"passed "<<std::endl;
    }
    ////////////////////////////////////////////////////////////////////////////////

    kinova->setGeneralizedVelocity(getVelocityCommand(gc, posDes.e(), quatDes.e()));
    world.integrate2();
    std::this_thread::sleep_for(std::chrono::microseconds(1000));
  }

  server.killServer();
}

#pragma once

struct robotDef {
  robotDef() {
    W_r_W2COM.setZero(7);
    W_J_COM.setZero(7);
    W_Jp_COM.setZero(7);
    W_Ja_COM.setZero(7);
    W_R.setZero(7);
    mass.setZero(7);
    inertia.setZero(7);
    for(auto& iner: inertia)
      iner.setZero();
  }

  // constants
  std::vector<mass> mass; // mass of each body
  std::vector<Matrix3d> B_I_COM; // inertia of each body
  std::vector<Vector3d> B_r_COM; // position of the com from the body frame
  std::vector<Vector3d> B_p; // joint axes in the body frame
  std::vector<Vector3d> B_r_j;
  std::vector<Matrix3d> B_R_j;

  // variables
  std::vector<Vector3d> W_p; // joint axes in the body frame
  std::vector<Vector3d> W_r_W2COM; // world to com position
  std::vector<Matrix<double,6,6>> W_J_COM; // world Jacobian of the COM
  std::vector<Matrix<double,3,6>> W_Jp_COM; // world position Jacobian of the COM
  std::vector<Matrix<double,3,6>> W_Ja_COM; // world angular Jacobian of the COM
  std::vector<Matrix3d> W_R; // orientation of each body
} robotDef;

void constructRobotDef(robotDef& robot) {
  robot.mass = {0.46784, 0.7477, 0.99, 0.6763, 0.463, 1.327, 0.01}

  robot.inertia[0].diagonal << 0.000951270861568, 0.000951270861568, 0.000374272;
  robot.inertia[1].diagonal << 0.00152031725204, 0.00152031725204, 0.00059816;
  robot.inertia[2].diagonal << 0.010502207991, 0.000792, 0.010502207991;
  robot.inertia[3].diagonal << 0.00142022431908, 0.000304335, 0.00142022431908;
  robot.inertia[4].diagonal << 0.0004321316048, 0.0004321316048, 9.26e-05;
  robot.inertia[5].diagonal << 0.0004321316048, 9.26e-05, 0.0004321316048;
  robot.inertia[6].diagonal << 0.0004403232387, 0.0004403232387, 0.0007416;

  for()

  joint_base.rpy = {0, 0, 0};
  joint_base.xyz = {0, 0, 0};
  joint_base.axis = {0, 0, 0};

  link_base.mass = 0.46784;
  link_base.rpy = {0, 0, 0};
  link_base.xyz = {0, 0, 0.1255};
  link_base.inertia = {0.000951270861568, 0.000951270861568, 0.000374272};

  joint_1.rpy = {0, 3.14159265359, 0};
  joint_1.xyz = {0, 0, 0.15675};
  joint_1.axis = {0, 0, 1};

  link_1.mass = 0.7477;
  link_1.rpy = {0, 0, 0};
  link_1.xyz = {0, -0.002, -0.0605};
  link_1.inertia = {0.00152031725204, 0.00152031725204, 0.00059816};

  joint_2.rpy = {-1.57079632679, 0, 3.14159265359};
  joint_2.xyz = {0, 0.0016, -0.11875};
  joint_2.axis = {0, 0, 1};

  link_2.mass = 0.99;
  link_2.rpy = {0, 0, 0};
  link_2.xyz = {0, -0.2065, -0.01};
  link_2.inertia = {0.010502207991, 0.000792, 0.010502207991};

  joint_3.rpy = {0, 3.14159265359, 0};
  joint_3.xyz = {0, -0.410, 0};
  joint_3.axis = {0, 0, 1};

  link_3.mass = 0.6763;
  link_3.rpy = {0, 0, 0};
  link_3.xyz = {0, 0.081, -0.0086};
  link_3.inertia = {0.00142022431908, 0.000304335, 0.00142022431908};

  joint_4.rpy = {-1.57079632679, 0, 3.14159265359};
  joint_4.xyz = {0, 0.2073, -0.0114};
  joint_4.axis = {0, 0, 1};

  link_4.mass = 0.463;
  link_4.rpy = {0, 0, 0};
  link_4.xyz = {0, 0.0028848942, -0.0541932613};
  link_4.inertia = {0.0004321316048, 0.0004321316048, 9.26e-05};

  joint_5.rpy = {1.57079632679, 0, 3.14159265359};
  joint_5.xyz = {0, 0, -0.10375};
  joint_5.axis = {0, 0, 1};

  link_5.mass = 0.463;
  link_5.rpy = {0, 0, 0};
  link_5.xyz = {0, 0.0497208855, -0.0028562765};
  link_5.inertia = {0.0004321316048, 9.26e-05, 0.0004321316048};

  joint_6.rpy = {-1.57079632679, 0, 3.14159265359};
  joint_6.xyz = {0, 0.10375, 0};
  joint_6.axis = {0, 0, 1};

  link_6.mass = 1.327;
  link_6.rpy = {0, 0, 0};
  link_6.xyz = {0, 0, -0.06};
  link_6.inertia = {0.0004403232387, 0.0004403232387, 0.0007416};

  joint_end.rpy = {3.14159265359, 0, 0};
  joint_end.xyz = {0, 0, -0.1600};
  joint_end.axis = {0, 0, 0};

  link_end.mass = 0.01;
  link_end.rpy = {0, 0, 0};
  link_end.xyz = {0, 0, 0};
  link_end.inertia = {0.01, 0.01, 0.01};
}

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {


    return Eigen::MatrixXd::Ones(6,6);
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return Eigen::VectorXd::Ones(6);
}


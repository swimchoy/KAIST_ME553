#pragma once

Eigen::Matrix3d FixedFrameRPY(const Eigen::Vector3d &rpy);
Eigen::Matrix3d rotation_Z(const float angle);

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition(const Eigen::VectorXd &gc) {
  /// given
  Eigen::MatrixXd xyz_set, rpy_set;
  xyz_set.setZero(7, 3), rpy_set.setZero(7, 3);
  xyz_set << 0.0, 0.0, 0.15675,
      0.0, 0.0016, -0.11875,
      0.0, -0.410, 0.0,
      0.0, 0.2073, -0.0114,
      0.0, 0.0, -0.10375,
      0.0, 0.10375, 0.0,
      0.0, 0.0, -0.1600;
  rpy_set << 0.0, 3.14159265359, 0.0,
      -1.57079632679, 0.0, 3.14159265359,
      0.0, 3.14159265359, 0.0,
      -1.57079632679, 0.0, 3.14159265359,
      1.57079632679, 0.0, 3.14159265359,
      -1.57079632679, 0.0, 3.14159265359,
      3.14159265359, 0.0, 0.0;

  Eigen::Vector3d positionVector;
  Eigen::Matrix3d rotationMat;

  positionVector.setZero();
  rotationMat.setIdentity();
  positionVector = positionVector + xyz_set.row(0).transpose();

  for (int idx = 0; idx < 6; ++idx) {
    rotationMat = rotationMat * FixedFrameRPY(rpy_set.row(idx)) * rotation_Z(gc(idx));
    positionVector = positionVector + rotationMat * xyz_set.row(idx + 1).transpose();
  };

  return positionVector; /// replace this
}

Eigen::Matrix3d rotation_X(const float angle) {
  Eigen::Matrix3d R;
  R << 1, 0, 0, 0, cos(angle), -sin(angle), 0, sin(angle), cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Y(const float angle) {
  Eigen::Matrix3d R;
  R << cos(angle), 0, sin(angle), 0, 1, 0, -sin(angle), 0, cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Z(const float angle) {
  Eigen::Matrix3d R;
  R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
  return R;
}

Eigen::Matrix3d FixedFrameRPY(const Eigen::Vector3d &rpy) {
  return rotation_Z(rpy(2)) * rotation_Y(rpy(1)) * rotation_X(rpy(0));
}
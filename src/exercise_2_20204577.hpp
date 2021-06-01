#pragma once
#include <vector>

Eigen::Matrix3d rotation_X(const double angle);
Eigen::Matrix3d rotation_Y(const double angle);
Eigen::Matrix3d rotation_Z(const double angle);
Eigen::Matrix3d FixedFrameRPY(const Eigen::Vector3d &rpy);
Eigen::Matrix3d QuaternionToRotMat(const Eigen::Vector4d &q);

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity(const Eigen::VectorXd &gc, const Eigen::VectorXd &gv) {
  Eigen::Vector3d footVel; footVel.setZero();
  Eigen::MatrixXd xyz_set; xyz_set.setZero(5, 3);
  xyz_set << 0.277, 0.116, 0.0,
      0.0635, 0.041, 0.0,
      0.0, 0.109, -0.25,
      0.1, -0.02, 0.0,
      0.0, 0.0, -0.32125;

  Eigen::Vector3d axisX, axisY, axisZ;
  axisX.setZero(); axisY.setZero(); axisZ.setZero();
  axisX(0) = 1; axisY(1) = 1; axisZ(2) = 1;

  std::vector<Eigen::Matrix3d> rotationMat;
  rotationMat.push_back(QuaternionToRotMat(gc.segment(3, 4)));
  rotationMat.push_back(rotation_X(gc(7)));
  rotationMat.push_back(rotation_Y(gc(8)));
  rotationMat.push_back(rotation_Y(gc(9)));

  Eigen::MatrixXd relativePositionVectors;  relativePositionVectors.setZero(6, 3);
  relativePositionVectors.row(0) = gc.head(3).transpose();
  relativePositionVectors.row(1) = (rotationMat[0] * xyz_set.row(0).transpose()).transpose();
  relativePositionVectors.row(2) = (rotationMat[0] * rotationMat[1] * xyz_set.row(1).transpose()).transpose();
  relativePositionVectors.row(3) = (rotationMat[0] * rotationMat[1] * rotationMat[2] * xyz_set.row(2).transpose()).transpose();
  relativePositionVectors.row(4) = (rotationMat[0] * rotationMat[1] * rotationMat[2] * rotationMat[3] * xyz_set.row(3).transpose()).transpose();
  relativePositionVectors.row(5) = (rotationMat[0] * rotationMat[1] * rotationMat[2] * rotationMat[3] * xyz_set.row(4).transpose()).transpose();

  std::vector<Eigen::Vector3d> endEffectorPosVectors; endEffectorPosVectors.clear();
  Eigen::Vector3d Zero3d; Zero3d.setZero();
  for (int idx = 0; idx < 6; ++idx) {
    endEffectorPosVectors.push_back(relativePositionVectors.colwise().sum());
    relativePositionVectors.row(idx) = Zero3d;
  }

  Eigen::Vector3d baseAngularVel=gv.segment(3,3);
  footVel += gv.head(3);
  footVel += baseAngularVel.cross(endEffectorPosVectors[1]);
  footVel += (gv(6) * rotationMat[0] * rotationMat[1] * axisX).cross(endEffectorPosVectors[2]);
  footVel += (gv(7) * rotationMat[0] * rotationMat[1] * rotationMat[2] * axisY).cross(endEffectorPosVectors[3]);
  footVel += (gv(8) * rotationMat[0] * rotationMat[1] * rotationMat[2] * rotationMat[3] * axisY).cross(endEffectorPosVectors[4]);

  return footVel; /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity(const Eigen::VectorXd &gc, const Eigen::VectorXd &gv) {
  Eigen::Vector3d footAngularVel; footAngularVel.setZero();

  Eigen::Vector3d axisX, axisY, axisZ;
  axisX.setZero(); axisY.setZero(); axisZ.setZero();
  axisX(0) = 1; axisY(1) = 1; axisZ(2) = 1;

  std::vector<Eigen::Matrix3d> rotationMat;
  rotationMat.push_back(QuaternionToRotMat(gc.segment(3, 4)));
  rotationMat.push_back(rotation_X(gc(7)));
  rotationMat.push_back(rotation_Y(gc(8)));
  rotationMat.push_back(rotation_Y(gc(9)));

  footAngularVel += gv.segment(3,3);
  footAngularVel += gv(6) * rotationMat[0] * axisX;
  footAngularVel += gv(7) * rotationMat[0] * rotationMat[1] * axisY;
  footAngularVel += gv(8) * rotationMat[0] * rotationMat[1] * rotationMat[2] * axisY;

  return footAngularVel; /// replace this
}

Eigen::Matrix3d rotation_X(const double angle) {
  Eigen::Matrix3d R;
  R << 1, 0, 0, 0, cos(angle), -sin(angle), 0, sin(angle), cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Y(const double angle) {
  Eigen::Matrix3d R;
  R << cos(angle), 0, sin(angle), 0, 1, 0, -sin(angle), 0, cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Z(const double angle) {
  Eigen::Matrix3d R;
  R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
  return R;
}

Eigen::Matrix3d FixedFrameRPY(const Eigen::Vector3d &rpy) {
  return rotation_Z(rpy(2)) * rotation_Y(rpy(1)) * rotation_X(rpy(0));
}

Eigen::Matrix3d QuaternionToRotMat(const Eigen::Vector4d &q) {
  Eigen::Matrix3d R;
    R << 1 - 2 * std::pow(q(2), 2) - 2 * std::pow(q(3), 2),
        2 * q(1) * q(2) - 2 * q(0) * q(3),
        2 * q(1) * q(3) + 2 * q(0) * q(2),
        2 * q(1) * q(2) + 2 * q(0) * q(3),
        1 - 2 * std::pow(q(1), 2) - 2 * std::pow(q(3), 2),
        2 * q(2) * q(3) - 2 * q(0) * q(1),
        2 * q(1) * q(3) - 2 * q(0) * q(2),
        2 * q(2) * q(3) + 2 * q(0) * q(1),
        1 - 2 * std::pow(q(1), 2) - 2 * std::pow(q(2), 2);

  return R;
}
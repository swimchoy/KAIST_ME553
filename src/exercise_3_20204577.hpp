#pragma once

void kinovaConfig (Eigen::MatrixXd &rpySet, Eigen::MatrixXd &xyzSet);
void rotationMat (std::vector<Eigen::Matrix3d> &rotMat, const Eigen::VectorXd& gc, const Eigen::MatrixXd &rpySet);
Eigen::MatrixXd relativePosVecSet (const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> rotMat);
Eigen::MatrixXd endEffectorPosVecSet (Eigen::MatrixXd &relativePosVecSet);
Eigen::MatrixXd fixedBaseJ (Eigen::MatrixXd endEffectorPosVecSet, std::vector<Eigen::Matrix3d> rotMat);
Eigen::Matrix3d QuaternionToRotMat (const Eigen::Vector4d &q);
Eigen::MatrixXd pseudoInv (Eigen::MatrixXd &A);
Eigen::Vector4d RotMatToQuaternion (const Eigen::Matrix3d &M);
Eigen::Vector3d QuaternionToRotVec (const Eigen::Vector4d &q);


/// do not change the name of the method
inline Eigen::VectorXd getVelocityCommand (const Eigen::VectorXd& gc, const Eigen::Vector3d& pos, const Eigen::Vector4d& quat) {
  Eigen::MatrixXd rpySet, xyzSet;
  rpySet.setZero(7,3); xyzSet.setZero(7,3);
  kinovaConfig(rpySet, xyzSet);

  std::vector<Eigen::Matrix3d> rotMat;
  rotationMat(rotMat, gc, rpySet);

  Eigen::MatrixXd relativePosVecSet_, endEffectorPosVecSet_, J_;
  relativePosVecSet_ = relativePosVecSet(xyzSet, rotMat);
  endEffectorPosVecSet_ = endEffectorPosVecSet(relativePosVecSet_);
  J_ = fixedBaseJ(endEffectorPosVecSet_, rotMat);

  Eigen::VectorXd Err(6); Err.setZero();
  Err.head(3) = pos - endEffectorPosVecSet_.row(0).transpose();
  Err.tail(3) = rotMat[6] * QuaternionToRotVec(RotMatToQuaternion(rotMat[6].transpose() * QuaternionToRotMat(quat)));

  return pseudoInv(J_) * Err;
}


void kinovaConfig (Eigen::MatrixXd &rpySet, Eigen::MatrixXd &xyzSet) {
  xyzSet << 0.0, 0.0, 0.15675,
      0.0, 0.0016, -0.11875,
      0.0, -0.410, 0.0,
      0.0, 0.2073, -0.0114,
      0.0, 0.0, -0.10375,
      0.0, 0.10375, 0.0,
      0.0, 0.0, -0.1600;
  rpySet << 0.0, 3.14159265359, 0.0,
      -1.57079632679, 0.0, 3.14159265359,
      0.0, 3.14159265359, 0.0,
      -1.57079632679, 0.0, 3.14159265359,
      1.57079632679, 0.0, 3.14159265359,
      -1.57079632679, 0.0, 3.14159265359,
      3.14159265359, 0.0, 0.0;
}


Eigen::Matrix3d rotation_X (const double angle) {
  Eigen::Matrix3d R;
  R << 1, 0, 0, 0, cos(angle), -sin(angle), 0, sin(angle), cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Y (const double angle) {
  Eigen::Matrix3d R;
  R << cos(angle), 0, sin(angle), 0, 1, 0, -sin(angle), 0, cos(angle);
  return R;
}

Eigen::Matrix3d rotation_Z (const double angle) {
  Eigen::Matrix3d R;
  R << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
  return R;
}

Eigen::Matrix3d QuaternionToRotMat (const Eigen::Vector4d &q) {
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

Eigen::Vector4d RotMatToQuaternion (const Eigen::Matrix3d &M) {
  double t, r;
  Eigen::Vector4d q;
  t = M(0,0) + M(1,1) + M(2,2);
  r = std::sqrt(1+t);
  q << 0.5 * r,
      std::copysign(0.5*std::sqrt(1+M(0,0)-M(1,1)-M(2,2)), M(2,1) - M(1,2)),
      std::copysign(0.5*std::sqrt(1-M(0,0)+M(1,1)-M(2,2)), M(0,2) - M(2,0)),
      std::copysign(0.5*std::sqrt(1-M(0,0)-M(1,1)+M(2,2)), M(1,0) - M(0,1));
  return q;
}

Eigen::Vector3d QuaternionToRotVec (const Eigen::Vector4d &q){
  if (q.tail(3).norm() > 1e-10) {
    Eigen::Vector3d a((1/q.tail(3).norm())*q.tail(3));
    return 2 * atan2(q.tail(3).norm(), q(0)) * a;
  } else {
    return Eigen::Vector3d::Zero();
  }
}

Eigen::Vector3d RotMatToRotVec (const Eigen::Matrix3d &R) {
  double theta, r, t;
  Eigen::Vector3d u;
  u << R(2,1) - R(1,2),
      R(0,2) - R(2,0),
      R(1,0) - R(0,1);
  r = u.norm();
  t = R(0,0) + R(1,1) + R(2,2);
  theta = std::atan2(r, t-1);
  return theta * u;
}

Eigen::MatrixXd pseudoInv (Eigen::MatrixXd &A) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::MatrixXd U = svd.matrixU();
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd Diag = svd.singularValues();
  Eigen::MatrixXd S; S.setZero(V.cols(), U.cols());
  for (int idx = 0; idx < Diag.size(); ++idx) {
    if (Diag(idx, 0) > 1e-8) {
      S(idx, idx) = 1 / Diag(idx, 0);
    } else {
      S(idx, idx) = 0;
    }
  }
  return V * S * U.transpose();
}

Eigen::Matrix3d FixedFrameRPY (const Eigen::Vector3d &rpy) {
  return rotation_Z(rpy(2)) * rotation_Y(rpy(1)) * rotation_X(rpy(0));
}

void rotationMat (std::vector<Eigen::Matrix3d> &rotMat, const Eigen::VectorXd& gc, const Eigen::MatrixXd &rpySet) {
  Eigen::Matrix3d rotMat_; rotMat_.setIdentity();
  for (int idx = 0; idx < gc.size(); ++idx) {
    rotMat_ = rotMat_ * FixedFrameRPY(rpySet.row(idx)) * rotation_Z(gc(idx));
    rotMat.push_back(rotMat_);
  }
  rotMat_ = rotMat_ * FixedFrameRPY(rpySet.row(gc.size()));
  rotMat.push_back(rotMat_);
}

Eigen::MatrixXd relativePosVecSet (const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> rotMat) {
  Eigen::MatrixXd rPosVecSet_; rPosVecSet_.setZero(7,3);
  rPosVecSet_.row(0) = xyzSet.row(0);
  for (int idx = 0; idx < 6; ++idx) {
    rPosVecSet_.row(idx+1) = (rotMat[idx] * xyzSet.row(idx+1).transpose()).transpose();
  }
  return rPosVecSet_;
}

Eigen::MatrixXd endEffectorPosVecSet (Eigen::MatrixXd &relativePosVecSet) {
  Eigen::MatrixXd endEffectorPosVecSet_; endEffectorPosVecSet_.setZero(7,3);
  Eigen::Vector3d Zero3d; Zero3d.setZero();
  for (int idx = 0; idx < 7; ++idx) {
    endEffectorPosVecSet_.row(idx) = relativePosVecSet.colwise().sum();
    relativePosVecSet.row(idx) = Zero3d;
  }
  return endEffectorPosVecSet_;
}

Eigen::MatrixXd fixedBasePosJ (Eigen::MatrixXd endEffectorPosVecSet, std::vector<Eigen::Matrix3d> rotMat) {
  Eigen::MatrixXd fixedBasePosJ_; fixedBasePosJ_.setZero(3,6);
  Eigen::Vector3d Vec1, Vec2;
  for (int idx = 0; idx < 6; ++idx) {
    Vec1 = endEffectorPosVecSet.row(idx).transpose(); Vec2 = rotMat[idx].col(2);
    fixedBasePosJ_.col(idx) = - Vec1.cross(Vec2);
  }
  return fixedBasePosJ_;
}

Eigen::MatrixXd fixedBaseOriJ (std::vector<Eigen::Matrix3d> rotMat) {
  Eigen::MatrixXd fixedBaseOriJ_; fixedBaseOriJ_.setZero(3,6);
  for (int idx = 0; idx < 6; ++idx)
    fixedBaseOriJ_.col(idx) = rotMat[idx].col(2);
  return fixedBaseOriJ_;
}

Eigen::MatrixXd fixedBaseJ (Eigen::MatrixXd endEffectorPosVecSet, std::vector<Eigen::Matrix3d> rotMat) {
  Eigen::MatrixXd J; J.setZero(6,6);
  J.topRows(3) = fixedBasePosJ(endEffectorPosVecSet, rotMat);
  J.bottomRows(3) = fixedBaseOriJ(rotMat);
  return J;
}
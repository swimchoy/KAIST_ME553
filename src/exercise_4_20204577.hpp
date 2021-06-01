#pragma once

void kinovaConfig (Eigen::MatrixXd &rpySet, Eigen::MatrixXd &xyzSet, Eigen::MatrixXd &comSet, Eigen::VectorXd &massSet, std::vector<Eigen::Matrix3d> &InertiaVec);
void changeInertiaFrame (const std::vector<Eigen::Matrix3d> &InertiaVec_B, const std::vector<Eigen::Matrix3d> &R_, std::vector<Eigen::Matrix3d> &InertiaVec_W);
std::vector<Eigen::Matrix3d> VectorizedRotMat (const Eigen::VectorXd &gc, const Eigen::MatrixXd &rpySet);
std::vector<std::vector<Eigen::Vector3d>> effectivePosVec (const Eigen::MatrixXd &xyzSet, const Eigen::MatrixXd &comSet, const std::vector<Eigen::Matrix3d> &rotMatSet);
void VectorizedJ (std::vector<Eigen::MatrixXd> &VecJ, const std::vector<std::vector<Eigen::Vector3d>> &effPosVec, const std::vector<Eigen::Matrix3d> &rotMatVec);
void spatialInertiaMat (std::vector<Eigen::MatrixXd> &Mi, const Eigen::VectorXd &massSet, const std::vector<Eigen::Matrix3d> &InertiaVec);
void VectorizedRdot(std::vector<Eigen::Matrix3d> &dR_, const std::vector<Eigen::Matrix3d> &R_, const Eigen::VectorXd &gv, const std::vector<Eigen::MatrixXd> &J_);
std::vector<std::vector<Eigen::Vector3d>> effectivePositionDotVec (const Eigen::MatrixXd &xyzSet, const Eigen::MatrixXd &comSet,const std::vector<Eigen::Matrix3d> &dR_);
void PositionJdot (std::vector<Eigen::MatrixXd> &dPosJ_, const std::vector<std::vector<Eigen::Vector3d>> &effectivePosVec, \
                  const std::vector<std::vector<Eigen::Vector3d>> effectivePosDotVec, const std::vector<Eigen::Matrix3d> &R_, const std::vector<Eigen::Matrix3d> &dR_);
void AngularJdot (std::vector<Eigen::MatrixXd> &dAngJ_, const std::vector<Eigen::Matrix3d> &dR_);
Eigen::Matrix3d skew (const Eigen::Vector3d &w);
std::vector<std::vector<Eigen::Vector3d>> effectivePositionDotVec2 (const std::vector<Eigen::MatrixXd> &J_, const std::vector<Eigen::MatrixXd> &Jframe_, \
                                                                    const std::vector<Eigen::Matrix3d> &R_, const Eigen::VectorXd &gv, const std::vector<std::vector<Eigen::Vector3d>> &effectivePosVec);
void getFrameJ (std::vector<Eigen::MatrixXd> &Jframe_, const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> &R_);


/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  Eigen::MatrixXd rpySet, xyzSet, comSet, massMat;
  Eigen::VectorXd massSet;
  std::vector<Eigen::Matrix3d> inertiaVec_B, inertiaVec_W, R_;
  std::vector<Eigen::MatrixXd> M_, J_;
  std::vector<std::vector<Eigen::Vector3d>> effPosVec;

  rpySet.setZero(7, 3); xyzSet.setZero(7,3); comSet.setZero(8,3);
  massMat.setZero(6,6); massSet.setZero(8);
  kinovaConfig(rpySet, xyzSet, comSet, massSet, inertiaVec_B);

  R_ = VectorizedRotMat(gc, rpySet);
  changeInertiaFrame(inertiaVec_B, R_, inertiaVec_W);
  spatialInertiaMat(M_, massSet, inertiaVec_W);

  effPosVec = effectivePosVec(xyzSet, comSet, R_);
  VectorizedJ(J_, effPosVec, R_);

  for (int i = 0; i <= 7; ++i) 
    massMat += J_[i].transpose() * M_[i] * J_[i];

  return massMat;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  Eigen::MatrixXd rpySet, xyzSet, comSet;
  Eigen::VectorXd massSet, NonlinearVec;
  std::vector<Eigen::Matrix3d> inertiaVec_B, inertiaVec_W, R_, dR_;
  std::vector<Eigen::MatrixXd> M_, J_, dPosJ_, dAngJ_, Jframe_;
  std::vector<std::vector<Eigen::Vector3d>> effPosVec, effPosDotVec, effPosDotVec2;

  rpySet.setZero(7, 3); xyzSet.setZero(7,3); comSet.setZero(8,3);
  massSet.setZero(8); NonlinearVec.setZero(6);
  kinovaConfig(rpySet, xyzSet, comSet, massSet, inertiaVec_B);

  R_ = VectorizedRotMat(gc, rpySet);
  changeInertiaFrame(inertiaVec_B, R_, inertiaVec_W);
  spatialInertiaMat(M_, massSet, inertiaVec_W);

  effPosVec = effectivePosVec(xyzSet, comSet, R_);
  VectorizedJ(J_, effPosVec, R_);

  VectorizedRdot(dR_, R_, gv, J_);
  effPosDotVec = effectivePositionDotVec(xyzSet, comSet, dR_);

  getFrameJ(Jframe_,xyzSet,R_);
  effPosDotVec2 = effectivePositionDotVec2(J_, Jframe_, R_, gv, effPosVec);

  PositionJdot(dPosJ_, effPosVec, effPosDotVec, R_, dR_);
  AngularJdot(dAngJ_, dR_);

  /// gravity
  std::vector<Eigen::Vector3d> g_;
  Eigen::Vector3d g;
  for (int i = 0; i <= 7; ++i) {
    g.setZero(); g[2] = massSet(i) * -9.81; g_.push_back(g);
  }

  for(int i = 0; i <= 7; ++i) {
    NonlinearVec += J_[i].topRows(3).transpose() * massSet(i) * dPosJ_[i] * gv + \
                    J_[i].bottomRows(3).transpose() * inertiaVec_W[i] * dAngJ_[i] * gv + \
                    J_[i].bottomRows(3).transpose() * skew(J_[i].bottomRows(3) * gv) * (inertiaVec_W[i] * J_[i].bottomRows(3) * gv) + \
                    -J_[i].topRows(3).transpose() * g_[i];
  }

  return NonlinearVec;
}

Eigen::Matrix3d inertiaMat (const Eigen::VectorXd i) {
  Eigen::Matrix3d I;
  I << i(0), i(1), i(2), i(1), i(3), i(4), i(2), i(4), i(5);
  return I;
}

void kinovaConfig (Eigen::MatrixXd &rpySet, Eigen::MatrixXd &xyzSet, Eigen::MatrixXd &comSet, Eigen::VectorXd &massSet, std::vector<Eigen::Matrix3d> &InertiaVec) {
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
  comSet << 0, 0, 0.1255,
      0, -0.002, -0.0605,
      0, -0.2065, -0.01,
      0, 0.081, -0.0086,
      0, 0.0028848942, -0.0541932613,
      0, 0.0497208855, -0.0028562765,
      0, 0, -0.06,
      0, 0, 0;
  massSet << 0.46784, 0.7477, 0.99, 0.6763, 0.463, 0.463, 1.327, 0.01;
  Eigen::MatrixXd inertiaConfig; inertiaConfig.setZero(8,6);
  inertiaConfig << 0.000951270861568, 0, 0, 0.000951270861568, 0, 0.000374272,
      0.00152031725204, 0, 0, 0.00152031725204, 0, 0.00059816,
      0.010502207991, 0, 0, 0.000792, 0, 0.010502207991,
      0.00142022431908, 0, 0, 0.000304335, 0, 0.00142022431908,
      0.0004321316048, 0, 0, 0.0004321316048, 0, 9.26e-05,
      0.0004321316048, 0, 0, 9.26e-05, 0, 0.0004321316048,
      0.0004403232387, 0, 0, 0.0004403232387, 0, 0.0007416,
      0.01, 0, 0, 0.01, 0, 0.01;
  for (int idx = 0; idx < 8; ++idx) {
    InertiaVec.push_back(inertiaMat(inertiaConfig.row(idx)));
  }
}

void changeInertiaFrame (const std::vector<Eigen::Matrix3d> &InertiaVec_B, const std::vector<Eigen::Matrix3d> &R_, std::vector<Eigen::Matrix3d> &InertiaVec_W) {
  InertiaVec_W.clear();
  for (int i = 0; i < InertiaVec_B.size(); i++)
    InertiaVec_W.push_back(R_[i] * InertiaVec_B[i] * R_[i].transpose());
}

void spatialInertiaMat (std::vector<Eigen::MatrixXd> &Mi, const Eigen::VectorXd &massSet, const std::vector<Eigen::Matrix3d> &InertiaVec) {
  Eigen::MatrixXd _M;
  Eigen::Matrix3d I3; I3.setIdentity();
  for (int idx = 0; idx < 8; ++idx) {
    _M.setZero(6,6);
    _M.topLeftCorner(3,3) = massSet(idx) * I3;
    _M.bottomRightCorner(3,3) = InertiaVec[idx];
    Mi.push_back(_M);
  }
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

Eigen::Matrix3d FixedFrameRPY (const Eigen::Vector3d &rpy) {
  return rotation_Z(rpy(2)) * rotation_Y(rpy(1)) * rotation_X(rpy(0));
}

std::vector<Eigen::Matrix3d> VectorizedRotMat (const Eigen::VectorXd &gc, const Eigen::MatrixXd &rpySet) {
  std::vector<Eigen::Matrix3d> rotMatSet;
  Eigen::Matrix3d rotMat;
  rotMat.setIdentity();
  rotMatSet.push_back(rotMat);
  for (int joint = 0; joint < gc.size(); ++joint) {
    rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * rotation_Z(gc(joint)) ;
    rotMatSet.push_back(rotMat);
  }
  rotMat = rotMat * FixedFrameRPY(rpySet.row(gc.size()));
  rotMatSet.push_back(rotMat);
  return rotMatSet;
}

Eigen::MatrixXd relativePositionVec (const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> &rotMatSet) {
  Eigen::MatrixXd relativePosMat;
  relativePosMat.setZero(7,3);
  for (int idx = 0; idx < relativePosMat.rows(); ++idx) {
    relativePosMat.row(idx) = rotMatSet[idx] * xyzSet.row(idx).transpose();
  }
  return relativePosMat;
}

Eigen::MatrixXd relativeMassPositionVec (const Eigen::MatrixXd &comSet, const std::vector<Eigen::Matrix3d> &rotMatSet) {
  Eigen::MatrixXd relMassPosMat;
  relMassPosMat.setZero(8,3);
  for (int idx = 0; idx < relMassPosMat.rows(); ++idx) {
    relMassPosMat.row(idx) = rotMatSet[idx] * comSet.row(idx).transpose();
  }
  return relMassPosMat;
}

std::vector<std::vector<Eigen::Vector3d>> effectivePosVec (const Eigen::MatrixXd &xyzSet, const Eigen::MatrixXd &comSet, const std::vector<Eigen::Matrix3d> &rotMatSet) {
  Eigen::MatrixXd relativePosMat;
  Eigen::MatrixXd relativeMassPosMat;
  Eigen::MatrixXd tmpPosMat;
  Eigen::Vector3d zero3d; zero3d.setZero();
  std::vector<Eigen::Vector3d> ePosVec;
  std::vector<std::vector<Eigen::Vector3d>> effectivePosVec;

  relativePosMat = relativePositionVec(xyzSet, rotMatSet);
  relativeMassPosMat = relativeMassPositionVec(comSet, rotMatSet);

  for (int target = 7; target >= 0 ; --target) {
    ePosVec.clear();
    tmpPosMat.setZero(target + 1,3);
    tmpPosMat.topRows(target) = relativePosMat.topRows(target);
    tmpPosMat.bottomRows(1) = relativeMassPosMat.row(target);

    for (int idx = 0; idx <= target; ++idx) {
      ePosVec.push_back(tmpPosMat.colwise().sum());
      tmpPosMat.row(idx) = zero3d.transpose();
    }
    effectivePosVec.insert(effectivePosVec.begin(),ePosVec);
  }

  return effectivePosVec;
}

void PositionJ (Eigen::MatrixXd &PosJ, const std::vector<Eigen::Vector3d> &ePosVec, const std::vector<Eigen::Matrix3d> &rotMatVec) {
  PosJ.setZero(3,6);
  for (int idx = 0; idx < ePosVec.size() - 1; ++idx) {
    if(idx < PosJ.cols()) { PosJ.col(idx) = -ePosVec[idx + 1].cross(rotMatVec[idx + 1].col(2)); }
  }
}

void AngularJ (Eigen::MatrixXd &AngJ, const std::vector<Eigen::Matrix3d> &rotMatVec, const int &target) {
  AngJ.setZero(3,6);
  for (int idx = 0; idx < target; ++idx) {
    if(idx < AngJ.cols()) { AngJ.col(idx) = rotMatVec[idx + 1].col(2); }
  }
}

Eigen::MatrixXd Jacobi (const std::vector<Eigen::Vector3d> &ePosVec, const std::vector<Eigen::Matrix3d> &rotMatVec) {
  Eigen::MatrixXd J, PosJ, AngJ;
  int target = ePosVec.size() - 1;
  J.setZero(6,6);
  PositionJ(PosJ, ePosVec, rotMatVec);
  AngularJ(AngJ, rotMatVec, target);
  J.topRows(3) = PosJ;
  J.bottomRows(3) = AngJ;
  return J;
}

void VectorizedJ (std::vector<Eigen::MatrixXd> &VecJ, const std::vector<std::vector<Eigen::Vector3d>> &effPosVec, const std::vector<Eigen::Matrix3d> &rotMatVec) {
  for (int target = 0; target <= 7; ++target) {
    VecJ.push_back(Jacobi(effPosVec[target], rotMatVec));
  }
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

Eigen::Matrix3d skew (const Eigen::Vector3d &w) {
  Eigen::Matrix3d S;
  S << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;
  return S;
}

void VectorizedRdot(std::vector<Eigen::Matrix3d> &dR_, const std::vector<Eigen::Matrix3d> &R_, const Eigen::VectorXd &gv, const std::vector<Eigen::MatrixXd> &J_) {
  dR_.clear();
  for (int i = 0; i < R_.size(); ++i) {
    dR_.push_back(skew(J_[i].bottomRows(3) * gv) * R_[i]);
  }
}

std::vector<std::vector<Eigen::Vector3d>> effectivePositionDotVec (const Eigen::MatrixXd &xyzSet, const Eigen::MatrixXd &comSet,const std::vector<Eigen::Matrix3d> &dR_) {
  Eigen::MatrixXd relativePosDot; relativePosDot.setZero(xyzSet.rows(),3);
  for (int i = 0; i < relativePosDot.rows(); ++i)
    relativePosDot.row(i) = dR_[i] * xyzSet.row(i).transpose();
  
  Eigen::MatrixXd relativeComPosDot; relativeComPosDot.setZero(comSet.rows(),3);
  for (int i = 0; i < relativeComPosDot.rows(); ++i)
    relativeComPosDot.row(i) = dR_[i] * comSet.row(i).transpose();

  Eigen::MatrixXd tmpPosDot;
  Eigen::Vector3d zero3d; zero3d.setZero();
  std::vector<Eigen::Vector3d> ePosDot;
  std::vector<std::vector<Eigen::Vector3d>> effectivePosDotVec;

  for (int target = 7; target >= 0; --target) {
    ePosDot.clear();
    tmpPosDot.setZero(target + 1, 3);
    tmpPosDot.topRows(target) = relativePosDot.topRows(target);
    tmpPosDot.bottomRows(1) = relativeComPosDot.row(target);

    for (int frame = 0; frame <= target; ++frame) {
      ePosDot.push_back(tmpPosDot.colwise().sum());
      tmpPosDot.row(frame) = zero3d.transpose();
    }  
    effectivePosDotVec.insert(effectivePosDotVec.begin(),ePosDot);
  }

  return effectivePosDotVec;
}

std::vector<std::vector<Eigen::Vector3d>> effectiveFramePosVec (const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> &rotMatSet) {
  Eigen::MatrixXd relativePosMat;
  Eigen::MatrixXd tmpPosMat;
  Eigen::Vector3d zero3d; zero3d.setZero();
  std::vector<Eigen::Vector3d> ePosVec;
  std::vector<std::vector<Eigen::Vector3d>> effectivePosVec;

  relativePosMat = relativePositionVec(xyzSet, rotMatSet);

  for (int target = 7; target > 0 ; --target) {
    ePosVec.clear();
    tmpPosMat.setZero(target,3);
    tmpPosMat.topRows(target) = relativePosMat.topRows(target);

    for (int idx = 0; idx < target; ++idx) {
      ePosVec.push_back(tmpPosMat.colwise().sum());
      tmpPosMat.row(idx) = zero3d.transpose();
    }
    effectivePosVec.insert(effectivePosVec.begin(),ePosVec);
  }

  ePosVec.clear();
  ePosVec.push_back(zero3d);
  effectivePosVec.insert(effectivePosVec.begin(), ePosVec);

  return effectivePosVec;
}

void getFrameJ (std::vector<Eigen::MatrixXd> &Jframe_, const Eigen::MatrixXd &xyzSet, const std::vector<Eigen::Matrix3d> &R_) {
  Jframe_.clear();
  std::vector<std::vector<Eigen::Vector3d>> effFramePosVec;
  effFramePosVec = effectiveFramePosVec(xyzSet, R_);
  for (int target = 0; target <= 7; ++target) {
    Jframe_.push_back(Jacobi(effFramePosVec[target], R_));
  }
}

std::vector<std::vector<Eigen::Vector3d>> effectivePositionDotVec2 (const std::vector<Eigen::MatrixXd> &J_, const std::vector<Eigen::MatrixXd> &Jframe_, \
const std::vector<Eigen::Matrix3d> &R_, const Eigen::VectorXd &gv, const std::vector<std::vector<Eigen::Vector3d>> &effectivePosVec) {
  std::vector<Eigen::Vector3d> effectivePosDot;
  std::vector<std::vector<Eigen::Vector3d>> effectivePosDotVec;
  for (int e = 7; e >= 0; --e) {
    effectivePosDot.clear();
    for (int i = 0; i <= e; ++i) {
      effectivePosDot.push_back((J_[e]*gv).head(3) - (Jframe_[i]*gv).head(3));
    }
    effectivePosDotVec.insert(effectivePosDotVec.begin(), effectivePosDot);
  }

  return effectivePosDotVec;
}

void PositionJdot (std::vector<Eigen::MatrixXd> &dPosJ_, const std::vector<std::vector<Eigen::Vector3d>> &effectivePosVec, \
const std::vector<std::vector<Eigen::Vector3d>> effectivePosDotVec, const std::vector<Eigen::Matrix3d> &R_, const std::vector<Eigen::Matrix3d> &dR_) {
  dPosJ_.clear();
  Eigen::MatrixXd dPosJ; 
  for (int e = 0; e <= 7; ++e) {
    dPosJ.setZero(3, 6);
    for (int i = 0; i < e; ++i) {
      if (i < dPosJ.cols())
        dPosJ.col(i) = - effectivePosDotVec[e][i+1].cross(R_[i+1].col(2)) - effectivePosVec[e][i+1].cross(dR_[i+1].col(2));
    }
    dPosJ_.push_back(dPosJ);
  }
}

void AngularJdot (std::vector<Eigen::MatrixXd> &dAngJ_, const std::vector<Eigen::Matrix3d> &dR_) {
  dAngJ_.clear();
  Eigen::MatrixXd dAngJ;
  for (int e = 0; e <= 7; ++e) {
    dAngJ.setZero(3, 6);
    for (int i = 0; i < e; ++i) {
      if (i < dAngJ.cols())
        dAngJ.col(i) = dR_[i+1].col(2);
    }
    dAngJ_.push_back(dAngJ);
  }
}

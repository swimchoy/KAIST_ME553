#pragma once

Eigen::MatrixXd pInv (const Eigen::MatrixXd &A) {
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

class Robot {

 public:

  Eigen::Matrix3d inertiaMat (const Eigen::VectorXd i) {
    Eigen::Matrix3d I;
    I << i(0), i(1), i(2), i(1), i(3), i(4), i(2), i(4), i(5);
    return I;
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

  Eigen::Matrix3d skew (const Eigen::Vector3d &w) {
    Eigen::Matrix3d S;
    S << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;
    return S;
  }

};

class KINOVA : public Robot {

 public:

  KINOVA () {
    rpySet.setZero(7,3);
    xyzSet.setZero(7,3);
    comSet.setZero(8,3);
    massSet.setZero(8);
    b_I_.clear();

    kinovaConfig();

    gc.setZero(6);
    gv.setZero(6);
  }

  void update (const Eigen::VectorXd &gc_, const Eigen::VectorXd &gv_) {
    gc = gc_;
    gv = gv_;

    updateKinematics();
    updateDynamics();
  }

  void kinovaConfig () {

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
      b_I_.push_back(inertiaMat(inertiaConfig.row(idx)));
    }
  }

  void vectorized_R () {

    R_.clear();

    Eigen::Matrix3d rotMat;
    rotMat.setIdentity();
    R_.push_back(rotMat);

    for (int joint = 0; joint < gc.size(); ++joint) {
      rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * rotation_Z(gc(joint)) ;
      R_.push_back(rotMat);
    }

    rotMat = rotMat * FixedFrameRPY(rpySet.row(gc.size()));
    R_.push_back(rotMat);
  }

  void relativeJointPositions () {

    relativeJointPos.setZero(7,3);

    for (int idx = 0; idx < relativeJointPos.rows(); ++idx) {
      relativeJointPos.row(idx) = R_[idx] * xyzSet.row(idx).transpose();
    }
  }

  void relativeComPositions () {

    relativeComPos.setZero(8,3);

    for (int idx = 0; idx < relativeComPos.rows(); ++idx) {
      relativeComPos.row(idx) = R_[idx] * comSet.row(idx).transpose();
    }
  }

  void framePositions () {

    framePos.setZero(8,3);
    zero3d.setZero();

    Eigen::MatrixXd tmp_relativePos;
    tmp_relativePos = relativeJointPos;

    for (int idx = 0; idx < tmp_relativePos.rows(); ++idx) {
      framePos.row(framePos.rows() - (idx+1)) = tmp_relativePos.colwise().sum();
      tmp_relativePos.row(tmp_relativePos.rows() - (idx+1)) = zero3d;
    }
  }

  void getBodyJ (const int &root_body, const Eigen::Vector3d &pos, Eigen::MatrixXd &J) {

    J.setZero(6,6);

    for (int i = 0; i < root_body; ++i) {
      if (i < 6) {
        J.topRows(3).col(i) = -skew(framePos.row(root_body) + pos.transpose() - framePos.row(i+1))*(R_[i + 1].col(2));
        J.bottomRows(3).col(i) = R_[i + 1].col(2);
      }
    }
  }

  void getBodyJdot (const int &root_body, const Eigen::Vector3d &pos, Eigen::MatrixXd &dJ) {

    Eigen::MatrixXd J;
    J.setZero(6,6);
    dJ.setZero(6,6);

    for (int i = 0; i < root_body; ++i) {
      getBodyJ(root_body, pos, J);
      if (i < 6) {
        dJ.topRows(3).col(i) = -skew((J * gv).head(3) - (frameJ_[i+1] * gv).head(3)) * (R_[i+1].col(2)) - \
                              skew(framePos.row(root_body) + pos.transpose() - framePos.row(i+1)) * (dR_[i+1].col(2));
        dJ.bottomRows(3).col(i) = dR_[i+1].col(2);
      }
    }
  }

  void vectorized_comJ () {

    comJ_.clear();
    Eigen::MatrixXd J;

    for (int e = 0; e <= 7; ++e) {
      getBodyJ(e, relativeComPos.row(e), J);
      comJ_.push_back(J);
    }
  }

  void vectorized_com_dJ () {

    com_dJ_.clear();
    Eigen::MatrixXd dJ;

    for (int e = 0; e <= 7; ++e) {
      getBodyJdot(e, relativeComPos.row(e), dJ);
      com_dJ_.push_back(dJ);
    }
  }

  void vectorized_frameJ () {

    frameJ_.clear();
    Eigen::MatrixXd J;

    J.setZero(6,6);
    frameJ_.push_back(J);

    for (int e = 1; e < 7; ++e) {
      getBodyJ(e-1, relativeJointPos.row(e-1), J);
      frameJ_.push_back(J);
    }

  }

  void vectorized_dR () {
    dR_.clear();

    for (int i = 0; i < R_.size(); ++i) {
      dR_.push_back(skew(comJ_[i].bottomRows(3) * gv) * R_[i]);
    }
  }

  void vectorized_GravityForce () {
    mg_.clear();
    Eigen::Vector3d mg;
    for (int i = 0; i <= 7; ++i) {
      mg.setZero();
      mg[2] = massSet(i) * -9.81;
      mg_.push_back(mg);
    }
  }

  void updateDynamics () {
    /// mass matrix
    MassMatrix.setZero(6,6);
    vectorized_comJ();

    for (int i = 0; i <= 7; ++i)
      MassMatrix += comJ_[i].transpose() * M_[i] * comJ_[i];

    /// nonlinear term
    Nonlinearities.setZero(6);

    vectorized_frameJ();
    vectorized_dR();
    vectorized_com_dJ();
    vectorized_GravityForce();

    for (int i = 0; i <= 7; ++i) {
      Nonlinearities += comJ_[i].topRows(3).transpose() * massSet(i) * com_dJ_[i].topRows(3) * gv + \
                        comJ_[i].bottomRows(3).transpose() * w_I_[i] * com_dJ_[i].bottomRows(3) * gv + \
                        comJ_[i].bottomRows(3).transpose() * skew(comJ_[i].bottomRows(3) * gv) * (w_I_[i] * comJ_[i].bottomRows(3) * gv) + \
                        -comJ_[i].topRows(3).transpose() * mg_[i];
    }

  }

  void updateKinematics () {
    /// calculate R
    vectorized_R();

    /// b_I to w_I
    w_I_.clear();
    changeInertiaFrame(b_I_, R_, w_I_);

    /// set spatialInertiaMatrix
    M_.clear();
    spatialInertiaMat(M_, massSet, w_I_);

    /// relative positions e.i. world
    relativeJointPositions();
    relativeComPositions();

    /// frame positions e.i. world
    framePositions();
  }

  Eigen::MatrixXd getMassMatrix () { return MassMatrix; }
  Eigen::MatrixXd getNonlinearities () { return Nonlinearities; }
  Eigen::Matrix3d getFrameOrientation (const int & e) { return R_[e]; }
  Eigen::Vector3d getFramePosition (const int & e) { return framePos.row(e); }
  Eigen::Vector3d getRelativeJointPos (const int & e) { return relativeJointPos.row(e); }
  Eigen::Vector3d getRelativeComPos (const int & e) { return relativeComPos.row(e); }
  Eigen::MatrixXd getComJ (const int & e) {return comJ_[e];}


 private:

  Eigen::Vector3d zero3d;
  Eigen::VectorXd gc, gv, massSet, Nonlinearities;
  Eigen::MatrixXd rpySet, xyzSet, comSet, MassMatrix, relativeJointPos, relativeComPos, framePos;
  std::vector<Eigen::Vector3d> mg_;
  std::vector<Eigen::Matrix3d> R_, dR_, b_I_, w_I_;
  std::vector<Eigen::MatrixXd> M_, comJ_, frameJ_, com_dJ_;

};

/// do not change the name of the method
inline Eigen::VectorXd getGeneralizedForce(const raisim::ArticulatedSystem *robot) {

  /// !!!!!!!!!! You can use RaiSim function this time. Or you can use your results from exercise 4

  KINOVA kinova;
  double r_ball = 0.05;
  Eigen::Vector3d contactPos, b_contactPos, contactF;
  Eigen::VectorXd gc, gv;
  Eigen::MatrixXd contactJ, contactJdot;
  double Kd = 0.5;

  /// update
  robot->getState(gc, gv);
  kinova.update(gc,gv);

  /// contact position
  contactPos.setZero();
  contactPos[2] = -r_ball;
  b_contactPos = kinova.getFrameOrientation(7).transpose() * contactPos;

  /// contact force
  contactF.setZero();
  contactF[2] = 10;

  /// contact J and Jdot
  kinova.getBodyJ(7, contactPos, contactJ);
  kinova.getBodyJdot(7, contactPos, contactJdot);

  /// return the tau command
  return  pInv(contactJ.topRows(3) * kinova.getMassMatrix().inverse()) * \
          (contactJ.topRows(3) * kinova.getMassMatrix().inverse() * kinova.getNonlinearities() - \
          contactJ.topRows(3) * kinova.getMassMatrix().inverse() * contactJ.topRows(3).transpose() * contactF - \
          contactJdot.topRows(3) * gv) + \
          Kd * (-gv);
}



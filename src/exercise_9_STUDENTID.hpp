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
    ga.setZero(6);
    gf.setZero(6);
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

    for (int e = 1; e <= 7; ++e) {
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
      // TODO: Change to -9.81
      mg[2] = massSet(i) * -0.00;
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

  /// Exercise 9 ///
  Eigen::VectorXd S (const int & idx) {
    Eigen::VectorXd s;
    s.setZero(6);
    if (idx < 7 && idx > 0)
      s.tail(3) = R_[idx].col(2);
    return s;
  }

  Eigen::VectorXd dS (const int & idx) {
    Eigen::VectorXd ds = S(idx);
    ds.tail(3) = skew(comJ_[idx].bottomRows(3) * gv) * ds.tail(3);
    return ds;
  }

  Eigen::MatrixXd X (const int & B, const int & P) {
    Eigen::MatrixXd x(6,6);
    x.setIdentity();
    x.bottomLeftCorner(3,3) = skew(framePos.row(B) - framePos.row(P));
    return x;
  }

  Eigen::MatrixXd dX (const int & B, const int & P) {
    Eigen::MatrixXd dx;
    dx.setZero(6,6);
    dx.bottomLeftCorner(3,3) = skew(frameJ_[B].topRows(3) * gv - frameJ_[P].topRows(3) * gv);
    return dx;
  }

  void updateSingleBodyDynamics_toJoint (const raisim::ArticulatedSystem *robot) {
    M_j.clear();
    b_j.clear();

    Eigen::MatrixXd M, X_com(6, 6);
    Eigen::VectorXd b;

    //TODO: MAKE Composite body at the end effector
    for (int i = 0; i <= 6; ++i) {
      M.setZero(6,6);
      b.setZero(6);

      if (i == 6) {
        double c_m = massSet(i) + massSet(i+1);
        Eigen::Vector3d r_com = (massSet(i) * relativeComPos.row(i) + massSet(i+1) * (relativeComPos.row(i+1) + (framePos.row(i+1) - framePos.row(i)))) / c_m;
        Eigen::Vector3d r1 = relativeComPos.row(i).transpose() - r_com;
        Eigen::Vector3d r2 = relativeComPos.row(i+1).transpose() + (framePos.row(i+1) - framePos.row(i)).transpose() - r_com;
        Eigen::Matrix3d c_I = w_I_[i] + w_I_[i+1] - massSet(i) * skew(r1) * skew(r1) - massSet(i+1) * skew(r2) * skew(r2);

        M.topLeftCorner(3,3) = c_m * Eigen::Matrix3d::Identity();
        M.topRightCorner(3,3) = -c_m * skew(-r_com);
        M.bottomLeftCorner(3,3) = c_m * skew(-r_com);
        M.bottomRightCorner(3,3) = c_I - c_m * skew(-r_com) * skew(-r_com);

        b.head(3) = c_m * skew(comJ_[i].bottomRows(3) * gv) * skew(comJ_[i].bottomRows(3) * gv) * (-r_com);
        b.tail(3) = skew(comJ_[i].bottomRows(3) * gv) * (c_I - c_m * skew(-r_com) * skew(-r_com)) * \
                  (comJ_[i].bottomRows(3) * gv);

      } else {
        M.topLeftCorner(3,3) = massSet(i) * Eigen::Matrix3d::Identity();
        M.topRightCorner(3,3) = -massSet(i) * skew(-relativeComPos.row(i));
        M.bottomLeftCorner(3,3) = massSet(i) * skew(-relativeComPos.row(i));
        M.bottomRightCorner(3,3) = w_I_[i] - massSet(i) * skew(-relativeComPos.row(i)) * skew(-relativeComPos.row(i));

        //TODO: Check should I use frameJ or comJ
        b.head(3) = massSet(i) * skew(comJ_[i].bottomRows(3) * gv) * skew(comJ_[i].bottomRows(3) * gv) * (-relativeComPos.row(i).transpose());
        b.tail(3) = skew(comJ_[i].bottomRows(3) * gv) * (w_I_[i] - massSet(i) * skew(-relativeComPos.row(i)) * skew(-relativeComPos.row(i))) * \
                  (comJ_[i].bottomRows(3) * gv);

        raisim::Vec<3> w;
        if (i!=0)
          robot->getFrameAngularVelocity("kinova_joint_" + std::to_string(i), w);
        else
          robot->getFrameAngularVelocity("kinova_joint_base", w);
        std::cout<<"GT: frame angVel"<<i<<"\n"<<w.e()<<std::endl;
        std::cout<<"body AngVel"<<i<<"\n"<<comJ_[i].bottomRows(3) * gv<<std::endl;
      }

      //TODO: gravity
//      X_com.setIdentity();
//      X_com.bottomLeftCorner(3,3) = skew(-relativeComPos.row(i));
//      b += - X_com.leftCols(3) * mg_[i];

      M_j.push_back(M);
      b_j.push_back(b);
    }
  }

  void setGeneralizedForce(const Eigen::VectorXd &gf_) {
    gf = gf_;
  }

  void updateArticulatedDynamics (const raisim::ArticulatedSystem *robot) {
    ArtMassMat.clear();
    ArtBiasedForce.clear();
    ArtMassMat.resize(7);
    ArtBiasedForce.resize(7);

    Eigen::MatrixXd ArtM;
    Eigen::VectorXd ArtB;
    ArtM.setZero(6,6);
    ArtB.setZero(6);

    Eigen::VectorXd ga_true(6);
    ga_true = MassMatrix.inverse() * (gf - Nonlinearities);
    std::cout<<"ga_true\n"<< ga_true<<std::endl;

    std::vector<Eigen::VectorXd> a_raisim;
    raisim::Vec<3> pos1, pos2;
    a_raisim.resize(7);
    Eigen::VectorXd a;
    a.setZero(6);
    a_raisim[0] = a;

    for (int i = 0; i < 6; ++i) {
      if (i != 0) {
        robot->getFrameVelocity("kinova_joint_" + std::to_string(i), pos1);
        robot->getFrameAngularVelocity("kinova_joint_" + std::to_string(i), pos2);
      } else {
        robot->getFrameVelocity("kinova_joint_base", pos1);
        robot->getFrameAngularVelocity("kinova_joint_base", pos2);
      }
      Eigen::VectorXd vp(6);
      vp.head(3) = pos1.e();
      vp.tail(3) = pos2.e();

      a = dS(i+1) * gv(i) + S(i+1) * ga_true(i) + dX(i+1, i).transpose() * vp + X(i+1, i).transpose() * a;
      a_raisim[i+1] = a;
      std::cout<<"a_raisim_"<<i+1<<a_raisim[i+1]<<std::endl;
    }

    for (int i = 0; i < M_j.size(); ++i) {

      std::cout<<"M_j_"<<i<<"\n"<<M_j[i]<<std::endl;
      std::cout<<"b_j_"<<i<<"\n"<<b_j[i]<<std::endl;
    }

    for (int i = 0; i < 6; ++i) {
      if (i != 0) {
        robot->getFramePosition("kinova_joint_" + std::to_string(i), pos1);
        robot->getFramePosition("kinova_joint_" + std::to_string(i+1), pos2);
      } else {
        robot->getFramePosition("kinova_joint_base", pos1);
        robot->getFramePosition("kinova_joint_" + std::to_string(i+1), pos2);
      }
      std::cout<<"skew(rPB)"<<i<<i+1<<"\n"<<skew(pos2.e()-pos1.e())<<std::endl;
    }

    for (int i = 6; i >= 0 ; --i) {
      std::cout<<"\nnow processing... : "<<i<<std::endl;

      std::cout<<"check Eq.5"<<std::endl;
      if (i != 6 && i != 0) {
        std::cout << ga_true(i) - (1 / (S(i+1).transpose() * ArtM * S(i+1))) * (-S(i+1).transpose() * (ArtM * (dS(i+1) * gv(i) + \
        dX(i+1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i - 1)) + X(i+1, i) * a_raisim[i]) + ArtB) + gf(i)) << std::endl;
      }

      if (i == 6) {
        ArtM = M_j[i];
        ArtB = b_j[i];
//      } else if (i == 6) {
//        std::cout<<S(i+1).transpose() * ArtM * S(i+1)<<std::endl;
//        ArtM = M_j[i] + X(i+1, i) * ArtM * X(i+1, i).transpose();
//        ArtB = b_j[i] + X(i+1, i) * (ArtM * (dX(i+1, i).transpose() * (frameJ_[i] * gv)) + ArtB);
      } else if (i == 0) {
        ArtM = M_j[i] + X(i+1, i) * ArtM * X(i+1, i).transpose() + X(i+1, i) * ArtM * S(i+1) * \
               (S(i+1).transpose() * ArtM * S(i+1)).inverse() * (-S(i+1).transpose() * ArtM * X(i+1, i).transpose());
        ArtB = b_j[i] + X(i+1, i) * (ArtM * (S(i+1) * (S(i+1).transpose() * ArtM * S(i+1)).inverse() * \
               (-S(i+1).transpose() * (ArtM * (dS(i+1) * gv(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv)) + ArtB) + gf(i)) + \
               dS(i+1) * gv(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv)) + ArtB);
      } else {
        ArtM = M_j[i] + X(i+1, i) * ArtM * X(i+1, i).transpose() + X(i+1, i) * ArtM * S(i+1) * \
               (S(i+1).transpose() * ArtM * S(i+1)).inverse() * (-S(i+1).transpose() * ArtM * X(i+1, i).transpose());
        ArtB = b_j[i] + X(i+1, i) * (ArtM * (S(i+1) * (S(i+1).transpose() * ArtM * S(i+1)).inverse() * \
               (-S(i+1).transpose() * (ArtM * (dS(i+1) * gv(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i-1))) + ArtB) + gf(i)) + \
               dS(i+1) * gv(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i-1))) + ArtB);

//        std::cout<<"SMS\n"<<S(i+1).transpose() * ArtM * S(i+1)<<std::endl;
//        std::cout<<"inverse\n"<<(S(i+1).transpose() * ArtM * S(i+1)).inverse()<<std::endl;

        // v_p check --> OK
//        raisim::Vec<3> pos1, pos2;
//        robot->getFrameVelocity("kinova_joint_" + std::to_string(i), pos1);
//        robot->getFrameAngularVelocity("kinova_joint_" + std::to_string(i), pos2);
//        Eigen::VectorXd vp(6);
//        vp.head(3) = pos1.e();
//        vp.tail(3) = pos2.e();
//        std::cout << "v_parent\n" << (frameJ_[i] * gv + S(i) * gv(i - 1) - vp).norm() << std::endl;


        // kinematics costraints check --> OK
//        std::cout<<"kinematics constraints"<<std::endl;
//        if (((frameJ_[i+1] * gv + S(i+1) * gv(i) - S(i+1) * gv(i)) - X(i+1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i-1))).norm() < 1e-8)
//          std::cout<<"passed!!"<<std::endl;

      }

      //TODO: problem
      std::cout<<"check primitive"<<std::endl;
      if (i != 0 && i != 6)
        std::cout<<gf(i-1) - S(i).transpose() * (M_j[i] * a_raisim[i] + b_j[i] + X(i+1, i) * \
                  (ArtMassMat[i+1] * a_raisim[i+1] + ArtBiasedForce[i+1]))<<std::endl;

      std::cout<<"check Articulated Dynamics"<<std::endl;
      if (i != 0)
        std::cout<<gf(i-1) - S(i).transpose() * (ArtM * a_raisim[i] + ArtB)<<std::endl;

//      std::cout<<"a_raisim\n"<<a_raisim[i]<<std::endl;
      std::cout<<"XBP\n"<<X(i+1, i)<<std::endl;
      std::cout<<"ST\n"<<S(i).transpose()<<std::endl;

      std::cout<<"M\n"<<ArtM<<std::endl;
      std::cout<<"b\n"<<ArtB<<std::endl;

      ArtMassMat[i] = ArtM;
      ArtBiasedForce[i] = ArtB;
    }

  }

  void ABA (Eigen::VectorXd & ga_, const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf, const raisim::ArticulatedSystem *robot) {
    update(gc, gv);
    setGeneralizedForce(gf);

    updateSingleBodyDynamics_toJoint(robot);
    updateArticulatedDynamics(robot);

    Eigen::VectorXd a_p;
    a_p.setZero(6);

    for (int i = 0; i < ga.size(); ++i) {
//      std::cout<<dS(i+1) * gv(i)<<std::endl;
      if (i > 0) {
        ga(i) = (1 / (S(i + 1).transpose() * ArtMassMat[i + 1] * S(i + 1))) * (-S(i + 1).transpose() * \
              (ArtMassMat[i + 1] * (dS(i + 1) * gv(i) + dX(i + 1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i - 1)) + \
              X(i + 1, i).transpose() * a_p) + ArtBiasedForce[i + 1]) + gf(i));
        a_p = dS(i+1) * gv(i) + S(i+1) * ga(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv + S(i) * gv(i-1)) + X(i+1, i) * a_p;
      } else {
        ga(i) = (1 / (S(i + 1).transpose() * ArtMassMat[i + 1] * S(i + 1))) * (-S(i + 1).transpose() * \
              (ArtMassMat[i + 1] * (dS(i + 1) * gv(i) + dX(i + 1, i).transpose() * (frameJ_[i] * gv) + \
              X(i + 1, i).transpose() * a_p) + ArtBiasedForce[i + 1]) + gf(i));
        a_p = dS(i+1) * gv(i) + S(i+1) * ga(i) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + X(i+1, i) * a_p;
      }
      // TODO: test eq 4
//      std::cout<<"eq4 test\n"<<gf(i) - S(i+1).transpose() * (ArtMassMat[i+1] * a_p + ArtBiasedForce[i+1])<<std::endl;


    }
    std::cout<<"R_[5]\n"<<R_[5]<<std::endl;
    ga_ = ga;
  }

 private:

  Eigen::Vector3d zero3d;
  Eigen::VectorXd gc, gv, massSet, Nonlinearities;
  Eigen::MatrixXd rpySet, xyzSet, comSet, MassMatrix, relativeJointPos, relativeComPos, framePos;
  std::vector<Eigen::Vector3d> mg_;
  std::vector<Eigen::Matrix3d> R_, dR_, b_I_, w_I_;
  std::vector<Eigen::MatrixXd> M_, comJ_, frameJ_, com_dJ_;

  // for Articulated Body Algorithm
  Eigen::VectorXd ga, gf;
  std::vector<Eigen::MatrixXd> ArtMassMat, M_j;
  std::vector<Eigen::VectorXd> ArtBiasedForce, b_j;
};

/// do not change the name of the method
inline Eigen::MatrixXd getGaUsingABA (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf, const raisim::ArticulatedSystem *robot) {
  KINOVA kinova;
  Eigen::VectorXd ga;
  kinova.ABA(ga, gc, gv, gf, robot);
  std::cout<<ga<<std::endl;



  return ga;
}
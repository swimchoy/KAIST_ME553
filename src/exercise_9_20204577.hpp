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
  for (int idx = 0; idx < massSet.size(); ++idx) {
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
  double t, r, s;
  Eigen::Vector4d q;
  t = M(0,0) + M(1,1) + M(2,2);
  r = std::sqrt(1+t);
  s = 1 / (2*r);
  q << 0.5 * r,
      (M(2,1) - M(1,2)) * s,
      (M(0,2) - M(2,0)) * s,
      (M(1,0) - M(0,1)) * s;
  return q;
}
Eigen::Matrix3d RotVecToRotMat (Eigen::Vector3d r) {
  Eigen::Matrix3d R;
  r += 1e-20 * Eigen::Vector3d::Ones();
  Eigen::Vector3d u = (1/r.norm()) * r;
  R << u[0]*u[0]*(1-std::cos(r.norm()))+std::cos(r.norm()),     u[0]*u[1]*(1-std::cos(r.norm()))-u[2]*std::sin(r.norm()), u[0]*u[2]*(1-std::cos(r.norm()))+u[1]*std::sin(r.norm()),
      u[1]*u[0]*(1-std::cos(r.norm()))+u[2]*std::sin(r.norm()), u[1]*u[1]*(1-std::cos(r.norm()))+std::cos(r.norm()),      u[1]*u[2]*(1-std::cos(r.norm()))-u[0]*std::sin(r.norm()),
      u[2]*u[0]*(1-std::cos(r.norm()))-u[1]*std::sin(r.norm()), u[2]*u[1]*(1-std::cos(r.norm()))+u[0]*std::sin(r.norm()), u[2]*u[2]*(1-std::cos(r.norm()))+std::cos(r.norm());
  return R;
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


class Robot {

 public:

  void update (const Eigen::VectorXd &gc_, const Eigen::VectorXd &gv_) {
    gc = gc_;
    gv = gv_;

    updateKinematics();
    updateDynamics();
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

  void updateDynamics () {
    /// update Jacobians
    vectorized_comJ();
    vectorized_frameJ();

    if (algorithm == "CRBA+RNE") {
      /// CRBA to get MassMatrix
      CompositeRigidBodyAlgorithm();

      /// CRBA to get NonlinearTerm
      runRNE_Nonlinearities();

    } else {
      /// default: PNE to get MassMatrix and NonlinearTerm
      projectedNewtonEuler();
    }
  }

  void projectedNewtonEuler () {
    /// mass matrix
    MassMatrix.setZero(dof,dof);

    for (int i = 0; i < massSet.size(); ++i)
      MassMatrix += comJ_[i].transpose() * M_[i] * comJ_[i];

    /// nonlinear term
    Nonlinearities.setZero(dof);

    vectorized_dR();
    vectorized_com_dJ();
    vectorized_GravityForce();

    for (int i = 0; i < massSet.size(); ++i) {
      Nonlinearities += comJ_[i].topRows(3).transpose() * massSet(i) * com_dJ_[i].topRows(3) * gv + \
                        comJ_[i].bottomRows(3).transpose() * w_I_[i] * com_dJ_[i].bottomRows(3) * gv + \
                        comJ_[i].bottomRows(3).transpose() * skew(comJ_[i].bottomRows(3) * gv) * (w_I_[i] * comJ_[i].bottomRows(3) * gv) + \
                        -comJ_[i].topRows(3).transpose() * mg_[i];
    }
  }

  void setAlgorithm (const std::string &algo) { algorithm = algo; }

  void vectorized_R () {
    R_.clear();
    Eigen::Matrix3d rotMat;
    int idx = 0;
    int axis_idx = 0;

    for (int joint = 0; joint < jointSet.size(); ++joint) {
      if (joint == 0) {
        /// at the root.
        if (floating) {
          rotMat = QuaternionToRotMat(gc.segment(3, 4));
          a_gc = gc.tail(gc.size() - 7);
        } else {
          rotMat.setIdentity();
          a_gc = gc;
        }
      } else {
        /// at the leaves.
        if (jointSet[joint] == "revolute") {
          if ((axisSet.row(axis_idx) - Eigen::Matrix3d::Identity().row(0)).norm() < 1e-4) {
            rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * rotation_X(a_gc(idx));
          } else if ((axisSet.row(axis_idx) - Eigen::Matrix3d::Identity().row(1)).norm() < 1e-4) {
            rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * rotation_Y(a_gc(idx));
          } else if ((axisSet.row(axis_idx) - Eigen::Matrix3d::Identity().row(2)).norm() < 1e-4) {
            rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * rotation_Z(a_gc(idx));
          } else {
            std::cout<<"the axis of rotation is not aligned with the joint frame."<<std::endl;
          }
          ++idx; ++axis_idx;
        } else if (jointSet[joint] == "prismatic") {
          rotMat = rotMat * FixedFrameRPY(rpySet.row(joint));
          ++idx; ++axis_idx;
        } else if (jointSet[joint] == "spherical") {
          rotMat = rotMat * FixedFrameRPY(rpySet.row(joint)) * QuaternionToRotMat(a_gc.segment(idx, 4));
          idx += 4; ++axis_idx;
        } else if (jointSet[joint] == "fixed") {
          rotMat = rotMat * FixedFrameRPY(rpySet.row(joint));
        } else {
          std::cout<<"this kind of joint is not yet defined"<<std::endl;
        }
      }
      R_.push_back(rotMat);
    }
  }

  void relativeJointPositions () {
    relativeJointPos.setZero(jointSet.size() - 1,3);
    int a_idx = 0;

    for (int idx = 0; idx < relativeJointPos.rows(); ++idx) {
      if (jointSet[idx+1] == "prismatic") {
        relativeJointPos.row(idx) = R_[idx] * xyzSet.row(idx+1).transpose() + R_[idx+1] * a_gc(a_idx+1) * axisSet.row(a_idx+1).transpose();
      } else {
        relativeJointPos.row(idx) = R_[idx] * xyzSet.row(idx+1).transpose();
      }
      if (jointSet[idx] != "fixed") { ++a_idx; }
    }
  }

  void relativeComPositions () {
    relativeComPos.setZero(massSet.size(),3);

    for (int idx = 0; idx < relativeComPos.rows(); ++idx) {
      relativeComPos.row(idx) = R_[idx] * comSet.row(idx).transpose();
    }
  }

  void framePositions () {
    framePos.setZero(jointSet.size(),3);

    Eigen::MatrixXd tmp_relativePos;
    tmp_relativePos = relativeJointPos;

    for (int idx = 0; idx < framePos.rows(); ++idx) {
      if (floating) {
        framePos.row(framePos.rows() - (idx + 1)) = tmp_relativePos.colwise().sum() + gc.head(3).transpose();
      } else {
        framePos.row(framePos.rows() - (idx + 1)) = tmp_relativePos.colwise().sum();
      }
      if (idx >= framePos.rows() - 1) { break; }
      tmp_relativePos.row(tmp_relativePos.rows() - (idx+1)) = Eigen::Vector3d::Zero();
    }
  }

  void getBodyJ (const int &target_body, const Eigen::Vector3d &pos, Eigen::MatrixXd &J) {
    J.setZero(6,gv.size());
    int a_idx = 0, i;

    for (int joint = 0; joint <= target_body; ++joint) {
      if (joint == 0) {
        if (floating) {
          J.topLeftCorner(3, 3) = Eigen::Matrix3d::Identity();
          J.block(0, 3, 3, 3) = -skew(framePos.row(target_body) + pos.transpose() - gc.head(3).transpose());
          J.bottomLeftCorner(3, 3) = Eigen::Matrix3d::Zero();
          J.block(3, 3, 3, 3) = Eigen::Matrix3d::Identity();
          i = 6;
        } else {
          i = 0;
        }
      } else {
        if (jointSet[joint] == "revolute") {
          J.topRows(3).col(i) = -skew(framePos.row(target_body) + pos.transpose() - framePos.row(joint)) * (R_[joint] * axisSet.row(a_idx).transpose());
          J.bottomRows(3).col(i) = R_[joint] * axisSet.row(a_idx).transpose();
          ++a_idx; ++i;
        } else if (jointSet[joint] == "prismatic") {
          J.topRows(3).col(i) = R_[joint] * axisSet.row(a_idx).transpose();
          J.bottomRows(3).col(i) = Eigen::Vector3d::Zero();
          ++a_idx; ++i;
        } else if (jointSet[joint] == "spherical") {
          J.topRows(3).middleCols(i, 3) = -skew(framePos.row(target_body) + pos.transpose() - framePos.row(joint)) * R_[joint];
          J.bottomRows(3).middleCols(i, 3) = R_[joint];
          ++a_idx; i += 3;
        }
      }

    }
  }

  void getBodyJdot (const int &target_body, const Eigen::Vector3d &pos, Eigen::MatrixXd &dJ) {
    Eigen::MatrixXd J;
    getBodyJ(target_body, pos, J);
    dJ.setZero(6,gv.size());
    int a_idx = 0, i;

    for (int joint = 0; joint <= target_body; ++joint) {
      if (joint == 0) {
        if (floating) {
          dJ.topLeftCorner(3, 3) = Eigen::Matrix3d::Zero();
          dJ.block(0, 3, 3, 3) = -skew(J.topRows(3) * gv - gv.head(3));
          dJ.bottomLeftCorner(3, 3) = Eigen::Matrix3d::Zero();
          dJ.block(3, 3, 3, 3) = Eigen::Matrix3d::Zero();
          i = 6;
        } else {
          i = 0;
        }
      } else {
        if (jointSet[joint] == "revolute") {
          dJ.topRows(3).col(i) = -skew(J.topRows(3) * gv - frameJ_[joint].topRows(3) * gv) * (R_[joint] * axisSet.row(a_idx).transpose()) -
              skew(framePos.row(target_body) + pos.transpose() - framePos.row(joint)) * (dR_[joint] * axisSet.row(a_idx).transpose());
          dJ.bottomRows(3).col(i) = dR_[joint] * axisSet.row(a_idx).transpose();
          ++a_idx; ++i;
        } else if (jointSet[joint] == "prismatic") {
          dJ.topRows(3).col(i) = dR_[joint] * axisSet.row(a_idx).transpose();
          dJ.bottomRows(3).col(i) = Eigen::Vector3d::Zero();
          ++a_idx; ++i;
        } else if (jointSet[joint] == "spherical") {
          dJ.topRows(3).middleCols(i, 3) = -skew(J.topRows(3) * gv - frameJ_[joint].topRows(3) * gv) * R_[joint] -
              skew(framePos.row(target_body) + pos.transpose() - framePos.row(joint)) * dR_[joint];
          dJ.bottomRows(3).middleCols(i, 3) = dR_[joint];
          ++a_idx; i += 3;
        }
      }
    }
  }

  void vectorized_comJ () {
    comJ_.clear();
    Eigen::MatrixXd J;

    for (int e = 0; e < massSet.size(); ++e) {
      getBodyJ(e, relativeComPos.row(e), J);
      comJ_.push_back(J);
    }
  }

  void vectorized_com_dJ () {
    com_dJ_.clear();
    Eigen::MatrixXd dJ;

    for (int e = 0; e < massSet.size(); ++e) {
      getBodyJdot(e, relativeComPos.row(e), dJ);
      com_dJ_.push_back(dJ);
    }
  }

  void vectorized_frameJ () {
    frameJ_.clear();
    Eigen::MatrixXd J;

    for (int i = 0; i < jointSet.size(); ++i) {
      getBodyJ(i, Eigen::Vector3d::Zero(), J);
      frameJ_.push_back(J);
    }
  }

  void vectorized_dR () {
    dR_.clear();

    for (int i = 0; i < R_.size(); ++i) {
      dR_.push_back(skew(frameJ_[i].bottomRows(3) * gv) * R_[i]);
    }
  }

  void vectorized_GravityForce () {
    mg_.clear();
    Eigen::Vector3d mg;
    for (int i = 0; i < massSet.size(); ++i) {
      mg.setZero();
      mg[2] = massSet(i) * gravity;
      mg_.push_back(mg);
    }
  }

  /// <Exercise 8> composite rigid body algorithm ///

  void CompositeRigidBodyAlgorithm () {
    MassMatrix.setZero(dof, dof);

    Eigen::MatrixXd M, a_j;
    Eigen::VectorXd b;
    Eigen::Vector3d r_ij;
    int idx_i = 0, idx_j = 0;

    for (int joint_j = 0; joint_j < jointSet.size(); ++joint_j) {
      compositeBodyDynamics_toJoint(joint_j, massSet.size() - 1, M, b);
      idx_i = 0;
      for (int joint_i = 0; joint_i <= joint_j; ++joint_i) {
        if((a_joint(joint_j) != -1) && (a_joint(joint_i) != -1)) {
          r_ij = relativeJointPos.middleRows(joint_i, joint_j - joint_i).colwise().sum();
          a_j = X(joint_j, joint_i).transpose() * S(joint_i);
          if (jointSet[joint_i] == "spherical" && jointSet[joint_j] == "spherical") {
            MassMatrix.bottomRightCorner(a_dof, a_dof).block(idx_j, idx_i, 3, 3) = S(joint_j).transpose() * (M * a_j);
            idx_i += 3;
          } else if (jointSet[joint_j] == "spherical") {
            MassMatrix.bottomRightCorner(a_dof, a_dof).block(idx_j, idx_i, 3, 1) = S(joint_j).transpose() * (M * a_j);
            ++idx_i;
          } else if (jointSet[joint_i] == "spherical") {
            MassMatrix.bottomRightCorner(a_dof, a_dof).block(idx_j, idx_i, 1, 3) = S(joint_j).transpose() * (M * a_j);
            idx_i += 3;
          } else {
            MassMatrix.bottomRightCorner(a_dof, a_dof)(idx_j, idx_i) = (S(joint_j).transpose() * (M * a_j)).value();
            ++idx_i;
          }
        }
      }
      if (jointSet[joint_j] == "spherical")
        idx_j += 3;
      else if (jointSet[joint_j] != "fixed")
        ++idx_j;
    }

    /// symmetric
    for (int j = 0; j < a_dof; ++j) {
      for (int i = 0; i < j; ++i) {
        MassMatrix.bottomRightCorner(a_dof, a_dof)(i, j) = MassMatrix.bottomRightCorner(a_dof, a_dof)(j, i);
      }
    }

    if (floating) {
      compositeBodyDynamics_toJoint(0, massSet.size()-1, M, b);
      MassMatrix.topLeftCorner(6,6) = M;
      for (int joint = 0; joint < jointSet.size(); ++joint) {
        compositeBodyDynamics_toJoint(joint, massSet.size() - 1, M, b);
        if (a_joint(joint) != -1)
          MassMatrix.bottomLeftCorner(a_dof, 6).row(a_joint(joint)) = S(joint).transpose() * M * X(joint, 0).transpose();
        /// symmetric
        MassMatrix.topRightCorner(6, a_dof) = MassMatrix.bottomLeftCorner(a_dof, 6).transpose();
      }
    }

  }

  void detectActingJoint () {
    a_joint.setZero(jointSet.size());
    int active = 0;
    for (int joint = 0; joint < jointSet.size(); ++joint) {
      if (jointSet[joint] == "fixed") {
        a_joint[joint] = -1;
      } else {
        a_joint[joint] = active;
        ++active;
      }
    }
  }

  void compositeBodyDynamics_toJoint(const int &start, const int &end, Eigen::MatrixXd &M, Eigen::VectorXd &b) {
    double m_c = 0.;
    Eigen::Vector3d r_com;
    Eigen::Matrix3d I_c;

    r_com.setZero();
    I_c.setZero();
    M.setZero(6,6);
    b.setZero(6);

    for (int i = start; i <= end; ++i) {
      r_com += massSet(i) * (framePos.row(i) + relativeComPos.row(i));
      m_c += massSet(i);
    }
    r_com /= m_c;

    for (int i = start; i <= end; ++i)
      I_c += w_I_[i] - massSet(i) * skew(framePos.row(i) + relativeComPos.row(i) - r_com.transpose()) * skew(framePos.row(i) + relativeComPos.row(i) - r_com.transpose());

    M.topLeftCorner(3,3) = m_c * Eigen::Matrix3d::Identity();
    M.topRightCorner(3,3) = -m_c * skew(r_com.transpose() - framePos.row(start));
    M.bottomLeftCorner(3,3) = m_c * skew(r_com.transpose() - framePos.row(start));
    M.bottomRightCorner(3,3) = I_c - m_c * skew(r_com.transpose() - framePos.row(start)) * skew(r_com.transpose() - framePos.row(start));

    b.head(3) = m_c * skew(frameJ_[start].bottomRows(3) * gv) * skew(frameJ_[start].bottomRows(3) * gv) * (r_com - framePos.row(start).transpose());
    b.tail(3) = skew(frameJ_[start].bottomRows(3) * gv) * M.bottomRightCorner(3,3) * (frameJ_[start].bottomRows(3) * gv);
  }

  /// <Exercise 8> recursive Newton Euler ///

  void runRNE_Nonlinearities () {
    Eigen::VectorXd ga_ = ga;
    Eigen::VectorXd gf_ = gf;
    ga.setZero(dof);
    gf.setZero(dof);

    recursiveNewtonEuler();
    Nonlinearities = gf;

    ga = ga_;
    gf = gf_;
  }

  void recursiveNewtonEuler () {
    Eigen::VectorXd accel(6);
    if (floating) {
      accel = ga.head(6);
      accel[2] += -gravity;
    } else {
      accel.setZero();
      accel[2] += -gravity;
    }

    /// first pass
    int gIdx = 0;
    std::vector<Eigen::VectorXd> acceleration;
    acceleration.resize(jointSet.size());

    if (floating)
      gIdx += 6;

    for (int joint = 0; joint < jointSet.size(); ++joint) {
      acceleration[joint] = accel;
      if (joint == jointSet.size() - 1) { break; }

      /// parent to child
      accel.head(3) += skew(accel.tail(3)) * relativeJointPos.row(joint).transpose() +
          skew(frameJ_[joint].bottomRows(3) * gv) * (frameJ_[joint+1].topRows(3) * gv - frameJ_[joint].topRows(3) * gv);

      /// i to i'
      if (jointSet[joint+1] != "fixed") {
        if (jointSet[joint+1] == "spherical") {
          accel += dS(joint+1) * gv.segment(gIdx, 3) + S(joint+1) * ga.segment(gIdx, 3);
          gIdx += 3;
        } else {
          accel += dS(joint + 1) * gv(gIdx) + S(joint + 1) * ga(gIdx);
          ++gIdx;
        }
      }
    }

    /// second pass
    int gfIdx;
    if (floating) {
      gfIdx = a_dof + 6;
    } else {
      gfIdx = a_dof;
    }
    Eigen::VectorXd f(6);
    f.setZero();

    updateSingleBodyDynamics_toJoint();

    for (int joint = M_j.size() - 1; joint >= 0; --joint) {
      if (joint == M_j.size() - 1) {
        f = M_j[joint] * acceleration[joint] + b_j[joint];
      } else {
        f = M_j[joint] * acceleration[joint] + b_j[joint] + X(joint+1, joint) * f;
      }
      --gfIdx;
      if (joint == 0) {
        if (floating) {
          gf.head(6) = f;
        }
        break;
      }
      if (jointSet[joint] == "spherical") {
        gfIdx -= 2;
        gf.segment(gfIdx, 3) = S(joint).transpose() * f;
      } else {
        gf(gfIdx) = (S(joint).transpose() * f).value();
      }
    }
  }

  void setGeneralizedAcceleration (const Eigen::VectorXd &ga_) {
    ga = ga_;
  }

  /// <Exercise 9> articulated body algorithm ///

  void ABA (Eigen::VectorXd & ga_, const Eigen::VectorXd& gc_, const Eigen::VectorXd& gv_, const Eigen::VectorXd& gf_) {
    update(gc_, gv_);
    setGeneralizedForce(gf_);

    /// to root: update M^a, b^a
    updateSingleBodyDynamics_toJoint();
    updateArticulatedDynamics();

    /// to leaves: calculate ga, a_p
    ForwardDynamics_ABA();
    ga_ = ga;
  }

  void setGeneralizedForce(const Eigen::VectorXd &gf_) {
    gf = gf_;
  }

  Eigen::MatrixXd S (const int & joint) {
    Eigen::MatrixXd s;
    if (jointSet[joint] == "revolute") {
      s.setZero(6, 1);
      s.bottomRows(3) = R_[joint] * axisSet.row(a_joint(joint)).transpose();
    } else if (jointSet[joint] == "prismatic") {
      s.setZero(6, 1);
      s.topRows(3) = R_[joint] * axisSet.row(a_joint(joint)).transpose();
    } else if (jointSet[joint] == "spherical") {
      s.setZero(6, 3);
      s.topRows(3) = Eigen::Matrix3d::Zero();
      s.bottomRows(3) = R_[joint];
    } else if (jointSet[joint] == "fixed") {
      s.setZero(6, 1);
      s = 1e-15 * Eigen::MatrixXd::Ones(6, 1);
    } else
      std::cout<<"this kind of joint is not yet provided"<<std::endl;
    return s;
  }

  Eigen::MatrixXd dS (const int & joint) {
    Eigen::MatrixXd ds = S(joint);
    if (jointSet[joint] == "revolute")
      ds.bottomRows(3) = skew(frameJ_[joint].bottomRows(3) * gv) * ds.bottomRows(3);
    else if (jointSet[joint] == "prismatic")
      ds.topRows(3) = skew(frameJ_[joint].bottomRows(3) * gv) * ds.topRows(3);
    else if (jointSet[joint] == "spherical")
      ds.bottomRows(3) = skew(frameJ_[joint].bottomRows(3) * gv) * ds.bottomRows(3);
    else if (jointSet[joint] == "fixed")
      ds = 1e-15 * Eigen::MatrixXd::Ones(6, 1);
    else
      std::cout<<"this kind of joint is not yet provided"<<std::endl;
    return ds;
  }

  Eigen::MatrixXd X (const int & B, const int & P) {
    Eigen::MatrixXd x(6,6);
    x.setIdentity();
    x.bottomLeftCorner(3, 3) = skew(framePos.row(B) - framePos.row(P));
    return x;
  }

  Eigen::MatrixXd dX (const int & B, const int & P) {
    Eigen::MatrixXd dx;
    dx.setZero(6,6);
    dx.bottomLeftCorner(3,3) = skew(frameJ_[B].topRows(3) * gv - frameJ_[P].topRows(3) * gv);
    return dx;
  }

  void updateSingleBodyDynamics_toJoint () {
    M_j.clear();
    b_j.clear();

    Eigen::MatrixXd M;
    Eigen::VectorXd b;
    int start, end;
    bool compositeChain = false;

    for (int joint = 0; joint < jointSet.size(); ++joint) {

      M.setZero(6,6);
      b.setZero(6);

      /// fixed joints
      if (joint < jointSet.size() - 1) {
        if (a_joint(joint + 1) == -1) {
          if (!compositeChain) { start = joint; }
          compositeChain = true;
          continue;
        } else {
          if (compositeChain) {
            end = joint;
            compositeBodyDynamics_toJoint(start, end, M, b);
            compositeChain = false;
            M_j.push_back(M);
            b_j.push_back(b);
            continue;
          }
        }
      } else {
        if (compositeChain) {
          end = joint;
          compositeBodyDynamics_toJoint(start, end, M, b);
          compositeChain = false;
          M_j.push_back(M);
          b_j.push_back(b);
          break;
        }
      }

      /// movable joints
      assert(!compositeChain);
      M.topLeftCorner(3,3) = massSet(joint) * Eigen::Matrix3d::Identity();
      M.topRightCorner(3,3) = -massSet(joint) * skew(relativeComPos.row(joint));
      M.bottomLeftCorner(3,3) = massSet(joint) * skew(relativeComPos.row(joint));
      M.bottomRightCorner(3,3) = w_I_[joint] - massSet(joint) * skew(relativeComPos.row(joint)) * skew(relativeComPos.row(joint));

      b.head(3) = massSet(joint) * skew(comJ_[joint].bottomRows(3) * gv) * skew(comJ_[joint].bottomRows(3) * gv) * (relativeComPos.row(joint).transpose());
      b.tail(3) = skew(comJ_[joint].bottomRows(3) * gv) * M.bottomRightCorner(3,3) * (comJ_[joint].bottomRows(3) * gv);

      M_j.push_back(M);
      b_j.push_back(b);
    }
  }

  void updateArticulatedDynamics () {
    ArtMassMat.resize(M_j.size());
    ArtBiasedForce.resize(b_j.size());

    Eigen::MatrixXd ArtM(6,6);
    Eigen::VectorXd ArtB(6);
    int j, B, gIdx;

    j = M_j.size();
    gIdx = gv.size();

    for (int i = jointSet.size()-1; i >= 0 ; --i) {
      if (jointSet[i] == "fixed" && i != 0) { continue; }

      /// Assume first joint is "fixed"
//      j = a_joint(i) + 1;
//      if (floating)
//        B = a_joint(i+1) + 6;
//      else
//        B = a_joint(i+1);
      --j;

      if (j == M_j.size() - 1) {
        ArtM = M_j[j];
        ArtB = b_j[j];
      } else {
        ArtM = M_j[j] + X(i + 1, i) * ArtMassMat[j + 1] * X(i + 1, i).transpose() + X(i + 1, i) * ArtMassMat[j + 1] * S(i + 1) * \
             (S(i + 1).transpose() * ArtMassMat[j + 1] * S(i + 1)).inverse() * (-S(i + 1).transpose() * ArtMassMat[j + 1] * X(i + 1, i).transpose());
        if (jointSet[i + 1] == "spherical") {
          gIdx -= 3;
          ArtB = b_j[j] + X(i + 1, i) * (ArtMassMat[j + 1] * (S(i + 1) * (S(i + 1).transpose() * ArtMassMat[j + 1] * S(i + 1)).inverse() * \
                 (-S(i + 1).transpose() * (ArtMassMat[j + 1] * (dS(i + 1) * gv.segment(gIdx, 3) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtBiasedForce[j + 1]) + \
                 gf.segment(gIdx, 3)) + dS(i + 1) * gv.segment(gIdx, 3) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtBiasedForce[j + 1]);
        } else {
          gIdx -= 1;
          ArtB = b_j[j] + X(i + 1, i) * (ArtMassMat[j + 1] * (S(i + 1) * (S(i + 1).transpose() * ArtMassMat[j + 1] * S(i + 1)).inverse() * \
                 ((-S(i + 1).transpose() * (ArtMassMat[j + 1] * (dS(i + 1) * gv(gIdx) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtBiasedForce[j + 1])).value() + \
                 gf(gIdx)) + dS(i + 1) * gv(gIdx) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtBiasedForce[j + 1]);
        }

//        ArtB = b_j[j] + X(i + 1, i) * (ArtMassMat[j + 1] * (S(i + 1) * (S(i + 1).transpose() * ArtMassMat[j + 1] * S(i + 1)).inverse() * \
//             (-S(i + 1).transpose() * (ArtMassMat[j + 1] * (dS(i + 1) * gv(B) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtBiasedForce[j + 1]) + \
//             gf(B)) + dS(i + 1) * gv(B) + dX(i + 1, i).transpose() * (frameJ_[i] * gv)) + ArtB);
      }

      ArtMassMat[j] = ArtM;
      ArtBiasedForce[j] = ArtB;
    }

    for (int i = 0; i < ArtMassMat.size(); ++i) {
      std::cout<<"M^a_"<<i<<"\n"<<ArtMassMat[i]<<std::endl;
    }
  }

  void ForwardDynamics_ABA () {
    Eigen::VectorXd a_p;
    int j, B;

    /// set accel at root frame
    if (floating) {
      a_p = ArtMassMat[0].inverse() * (gf.head(6) - ArtBiasedForce[0]);
    } else {
      a_p.setZero(6);
      a_p[2] += -gravity;
    }

    /// set a_p as ga.head(6) if it's floating
    if (floating) {
      a_p[2] += gravity;
      ga.head(6) = a_p;
      a_p[2] -= gravity;
    }

    /// to the leaves
    B = 0;
    if (floating)
      B = 6;

    for (int i = 0; i < jointSet.size(); ++i) {
      if (i == jointSet.size()-1) { break; }
      if (jointSet[i+1] == "fixed") { continue; }
      j = a_joint(i) + 1;
//      if (floating)
//        B = a_joint(i+1) + 6;
//      else
//        B = a_joint(i+1);

      if (jointSet[i+1] == "spherical") {
        ga.segment(B, 3) = (S(i+1).transpose() * ArtMassMat[j+1] * S(i+1)).inverse() * (-S(i+1).transpose() * \
                       (ArtMassMat[j+1] * (dS(i+1) * gv.segment(B, 3) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + \
                       X(i+1, i).transpose() * a_p) + ArtBiasedForce[j+1]) + gf.segment(B,3));
        a_p = dS(i+1) * gv.segment(B, 3) + S(i+1) * ga.segment(B, 3) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + X(i+1, i).transpose() * a_p;
        B += 3;
      } else {
        ga(B) = (S(i+1).transpose() * ArtMassMat[j+1] * S(i+1)).inverse()(0, 0) * ((-S(i+1).transpose() * \
                 (ArtMassMat[j+1] * (dS(i+1) * gv(B) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + \
                 X(i+1, i).transpose() * a_p) + ArtBiasedForce[j+1])).value() + gf(B));
        a_p = dS(i+1) * gv(B) + S(i+1) * ga(B) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + X(i+1, i).transpose() * a_p;
        ++B;
      }
//      ga(B) = (1 / (S(i+1).transpose() * ArtMassMat[j+1] * S(i+1))) * (-S(i+1).transpose() * \
//                       (ArtMassMat[j+1] * (dS(i+1) * gv(B) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + \
//                       X(i+1, i).transpose() * a_p) + ArtBiasedForce[j+1]) + gf(B));
//      a_p = dS(i+1) * gv(B) + S(i+1) * ga(B) + dX(i+1, i).transpose() * (frameJ_[i] * gv) + X(i+1, i).transpose() * a_p;
    }
  }

  int getDof() { return dof; }
  Eigen::MatrixXd getMassMatrix () { return MassMatrix; }
  Eigen::MatrixXd getNonlinearities () { return Nonlinearities; }
  Eigen::Matrix3d getFrameOrientation (const int & e) { return R_[e]; }
  Eigen::Vector3d getFramePosition (const int & e) { return framePos.row(e); }
  Eigen::MatrixXd getFrameJacobian (const int & j) { return frameJ_[j]; }
  Eigen::Vector3d getRelativeJointPos (const int & e) { return relativeJointPos.row(e); }
  Eigen::Vector3d getRelativeComPos (const int & e) { return relativeComPos.row(e); }
  Eigen::MatrixXd getComJ (const int & e) {return comJ_[e];}


 protected:

  int dof=0, a_dof=0;
  double gravity=0;
  bool floating=false;
  std::string algorithm = "PNE";
  Eigen::VectorXd gc, gv, massSet, Nonlinearities, a_gc, a_joint;
  Eigen::MatrixXd rpySet, xyzSet, comSet, axisSet, MassMatrix, relativeJointPos, relativeComPos, framePos;
  std::vector<std::string> jointSet;
  std::vector<Eigen::Vector3d> mg_;
  std::vector<Eigen::Matrix3d> R_, dR_, b_I_, w_I_;
  std::vector<Eigen::MatrixXd> M_, comJ_, frameJ_, com_dJ_;

  /// Articulated Body Algorithm
  Eigen::VectorXd ga, gf;
  std::vector<Eigen::MatrixXd> ArtMassMat, M_j;
  std::vector<Eigen::VectorXd> ArtBiasedForce, b_j;
};


class ANYMAL_ONELEG : public Robot {

 public:

  ANYMAL_ONELEG () {
    a_dof = 3;
    dof = 6 + a_dof;
    floating = true;

    gravity = -9.81;

    xyzSet.setZero(6,3);
    rpySet.setZero(6,3);
    axisSet.setZero(3,3);
    jointSet.resize(6);

    comSet.setZero(6,3);
    massSet.setZero(6);
    b_I_.clear();

    anymalOneLegConfig();
    detectActingJoint();

    gc.setZero(dof + 1);
    gv.setZero(dof);
    ga.setZero(dof);
    gf.setZero(dof);
  }

  void anymalOneLegConfig () {
    /// joint configuration
    xyzSet << 0, 0, 0,
        0.277, 0.116, 0.0,
        0.0635, 0.041, 0.0,
        0.0, 0.109, -0.25,
        0.1, -0.02, 0.0,
        0.0, 0.0, -0.32125;

    rpySet << 0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0,
        0, 0, 0;

    axisSet << 1, 0, 0,
        0, 1, 0,
        0, 1, 0;

    jointSet[0] = "fixed";
    jointSet[1] = "revolute";
    jointSet[2] = "revolute";
    jointSet[3] = "revolute";
    jointSet[4] = "fixed";
    jointSet[5] = "fixed";

    /// link configuration
    comSet << -0.001960558279, -0.001413217745, 0.050207125344,
        0.064516258147, -0.003787101702, -0.000152184388,
        -0.003897968082, 0.054226618537, -0.214583373795,
        0.030816858139, -0.004617229294, 0.000893125713,
        -8.66e-10, -1.472e-09, -0.244345749188,
        0., 0., 0.;

    massSet << 16.793507758, 1.42462064, 1.634976467, 0.207204302, 0.140170767, 0.;

    Eigen::MatrixXd inertiaConfig; inertiaConfig.setZero(massSet.size(),6);
    inertiaConfig <<
                  0.217391101503, -0.00132873239126, -0.00228200226173, 0.639432546734, -0.00138078263145, 0.62414077654,
        0.00243023349564, -1.53023971e-05, -2.1819095354e-05, 0.00230257239103, 2.6473021273e-05, 0.0019806759227,
        0.0120367944369, 6.762065206e-05, 0.000287806340448, 0.0120643637939, -0.00140610131218, 0.00249422574881,
        0.0002104880248, -5.6750980345e-05, 1.0127699391e-05, 0.000676270210023, -8.22869024e-07, 0.000545032674924,
        0.00159938741862, -9.32e-13, 1.039e-11, 0.00159938741932, 1.7563e-11, 5.4423177329e-05,
        0., 0., 0., 0., 0., 0.;

    for (int idx = 0; idx < inertiaConfig.rows(); ++idx) {
      b_I_.push_back(inertiaMat(inertiaConfig.row(idx)));
    }
  }

};

class KINOVA : public Robot {

 public:

  KINOVA () {
    a_dof = 6;
    dof = 6;
    floating = false;

    gravity = -9.81;

    xyzSet.setZero(8,3);
    rpySet.setZero(8,3);
    axisSet.setZero(6,3);
    jointSet.resize(8);

    comSet.setZero(8,3);
    massSet.setZero(8);
    b_I_.clear();

    kinovaConfig();
    detectActingJoint();

    int gcDim = a_dof;
    if (floating)
      gcDim += 1;
    gcDim += std::count(jointSet.begin(), jointSet.end(), "spherical");

    gc.setZero(gcDim);
    gv.setZero(a_dof);
    ga.setZero(a_dof);
    gf.setZero(a_dof);
  }

  void kinovaConfig () {

    xyzSet << 0.0, 0.0, 0.0,
        0.0, 0.0, 0.15675,
        0.0, 0.0016, -0.11875,
        0.0, -0.410, 0.0,
        0.0, 0.2073, -0.0114,
        0.0, 0.0, -0.10375,
        0.0, 0.10375, 0.0,
        0.0, 0.0, -0.1600;

    rpySet << 0.0, 0.0, 0.0,
         0.0, 3.14159265359, 0.0,
        -1.57079632679, 0.0, 3.14159265359,
        0.0, 3.14159265359, 0.0,
        -1.57079632679, 0.0, 3.14159265359,
        1.57079632679, 0.0, 3.14159265359,
        -1.57079632679, 0.0, 3.14159265359,
        3.14159265359, 0.0, 0.0;

    axisSet << 0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1;

    jointSet[0] = "fixed";
    jointSet[1] = "revolute";
    jointSet[2] = "revolute";
    jointSet[3] = "revolute";
    jointSet[4] = "revolute";
    jointSet[5] = "revolute";
    jointSet[6] = "revolute";
    jointSet[7] = "fixed";

    comSet << 0, 0, 0.1255,
        0, -0.002, -0.0605,
        0, -0.2065, -0.01,
        0, 0.081, -0.0086,
        0, 0.0028848942, -0.0541932613,
        0, 0.0497208855, -0.0028562765,
        0, 0, -0.06,
        0, 0, 0;

    massSet << 0.46784, 0.7477, 0.99, 0.6763, 0.463, 0.463, 1.327, 0.01;

    Eigen::MatrixXd inertiaConfig; inertiaConfig.setZero(massSet.size(),6);
    inertiaConfig << 0.000951270861568, 0, 0, 0.000951270861568, 0, 0.000374272,
        0.00152031725204, 0, 0, 0.00152031725204, 0, 0.00059816,
        0.010502207991, 0, 0, 0.000792, 0, 0.010502207991,
        0.00142022431908, 0, 0, 0.000304335, 0, 0.00142022431908,
        0.0004321316048, 0, 0, 0.0004321316048, 0, 9.26e-05,
        0.0004321316048, 0, 0, 9.26e-05, 0, 0.0004321316048,
        0.0004403232387, 0, 0, 0.0004403232387, 0, 0.0007416,
        0.01, 0, 0, 0.01, 0, 0.01;

    for (int idx = 0; idx < inertiaConfig.rows(); ++idx) {
      b_I_.push_back(inertiaMat(inertiaConfig.row(idx)));
    }
  }

};

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrixUsingCRBA (const Eigen::VectorXd& gc, raisim::ArticulatedSystem * robot) {

  KINOVA kinova;

  kinova.setAlgorithm("CRBA+RNE");
  kinova.update(gc, Eigen::VectorXd::Zero(kinova.getDof()));

  return kinova.getMassMatrix();
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearitiesUsingRNE (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, raisim::ArticulatedSystem * robot) {

  KINOVA kinova;

  kinova.setAlgorithm("CRBA+RNE");
  kinova.update(gc, gv);

  return kinova.getNonlinearities();
}

inline Eigen::MatrixXd getGaUsingABA (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {

  KINOVA kinova;
  Eigen::VectorXd ga;
  kinova.ABA(ga, gc, gv, gf);

  return ga;
}


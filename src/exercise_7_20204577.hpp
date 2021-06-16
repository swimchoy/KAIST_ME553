#pragma once
#include <set>

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

//      std::copysign(0.5*std::sqrt(1+M(0,0)-M(1,1)-M(2,2)), M(2,1) - M(1,2)),
//      std::copysign(0.5*std::sqrt(1-M(0,0)+M(1,1)-M(2,2)), M(0,2) - M(2,0)),
//      std::copysign(0.5*std::sqrt(1-M(0,0)-M(1,1)+M(2,2)), M(1,0) - M(0,1));

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

struct Ball {
  double radius, mass;
  Eigen::Matrix3d b_I, w_I;
  Eigen::VectorXd gc, gv, ga, gf, Nonlinear;
  Eigen::MatrixXd MassMatrix;
};

class Contact {
 public:

  Contact(const int & obj1, const int & obj2) : obj1(obj1), obj2(obj2) {
    assert(obj1 != -1);
    if (obj2 == -1) {
      groundContact=true;
      dynamicContact=false;
    } else {
      dynamicContact=true;
      groundContact=false;
    }
  }

  void update(const std::vector<Ball> &ball) {
    if (groundContact) {
      distance = ball[obj1].gc[2] - ball[obj1].radius;

      /// contact positions
      contactPos1 = ball[obj1].gc.head(3);
      contactPos1[2] -= ball[obj1].radius;
      contactPos2.setZero();
      contactPos2.head(2) = ball[obj1].gc.head(2);

      /// contact frame & Jacobi
      Rc.setIdentity();
      updateContactJ(ball);

      /// apparent Mass Matrix
      MappInv1.setZero(3,3);
      MappInv1 = Jc1 * ball[obj1].MassMatrix.inverse() * Jc1.transpose();
      MappInv2.setZero(3,3);

    } else if (dynamicContact) {
      distance = (ball[obj1].gc.head(3) - ball[obj2].gc.head(3)).norm() - ball[obj1].radius - ball[obj2].radius;

      /// contact positions w.r.t. obj1 and obj2
      contactPos1 = ball[obj1].gc.head(3) + (ball[obj2].gc.head(3) - ball[obj1].gc.head(3)) * \
                    (ball[obj1].radius / (ball[obj1].radius + ball[obj2].radius));
      contactPos2 = ball[obj2].gc.head(3) + (ball[obj1].gc.head(3) - ball[obj2].gc.head(3)) * \
                    (ball[obj2].radius / (ball[obj1].radius + ball[obj2].radius));

      /// contact frame: Gram-Schmidt Process.
      Rc.col(2) = (ball[obj1].gc.head(3) - ball[obj2].gc.head(3));
      Rc.col(0) = QuaternionToRotMat(ball[obj1].gc.tail(4)).col(0) - (QuaternionToRotMat(ball[obj1].gc.tail(4)).col(0).dot(Rc.col(2)) / Rc.col(2).dot(Rc.col(2))) * Rc.col(2);
      if ((Rc.col(0) - Eigen::Vector3d::Zero()).norm() < 1e-8) {
        Rc.col(0) = (QuaternionToRotMat(ball[obj1].gc.tail(4)).col(0) + 1e-6*Eigen::Vector3d::Ones()) - \
            ((QuaternionToRotMat(ball[obj1].gc.tail(4)).col(0) + 1e-6*Eigen::Vector3d::Ones()).dot(Rc.col(2)) / Rc.col(2).dot(Rc.col(2))) * Rc.col(2);
      }
      Rc.col(2) /= Rc.col(2).norm();
      Rc.col(0) /= Rc.col(0).norm();
      Rc.col(1) = Rc.col(2).cross(Rc.col(0));

      assert(Rc.col(0).norm() - 1. < 1e-8);
      assert(Rc.col(1).norm() - 1. < 1e-8);
      assert(Rc.col(2).norm() - 1. < 1e-8);
      assert(Rc.col(0).dot(Rc.col(1)) < 1e-8);
      assert(Rc.col(1).dot(Rc.col(2)) < 1e-8);
      assert(Rc.col(2).dot(Rc.col(0)) < 1e-8);

      /// Jacobi
      updateContactJ(ball);

      /// apparent Mass matrix
      MappInv1.setZero(3,3);
      MappInv1 = Jc1 * ball[obj1].MassMatrix.inverse() * Jc1.transpose();
      MappInv2.setZero(3,3);
      MappInv2 = Jc2 * ball[obj2].MassMatrix.inverse() * Jc2.transpose();

    } else {
      std::cout<<"WARNING: both contact booleans are false"<<std::endl;
    }
  }

  void updateContactJ(const std::vector<Ball> &ball) {

    Jc1.setZero(3,6);
    Jc1.block(0,0,3,3).setIdentity();
    Jc1.block(0,3,3,3) = -skew(contactPos1 - ball[obj1].gc.head(3));

    Jc2.setZero(3, 6);
    if (dynamicContact) {
      Jc2.block(0, 0, 3, 3).setIdentity();
      Jc2.block(0, 3, 3, 3) = -skew(contactPos2 - ball[obj2].gc.head(3));
    }
  }

  bool isInCollision(const std::vector<Ball> &ball) {
    if (distance <= 0.) {
      return true;
    } else {
      return false;
    }
  }

  void reset(const std::vector<Ball> &ball) {
    Imp.setZero();
    if (groundContact) {
      prevVimp = Rc.transpose() * (Jc1 * ball[obj1].gv);
      prevVimp[2] = (1+c_r)*prevVimp[2];
    } else if (dynamicContact) {
      prevVimp = Rc.transpose() * (Jc1 * ball[obj1].gv - Jc2 * ball[obj2].gv);
      prevVimp[2] = (1+c_r)*prevVimp[2];
    }
  }

  double mu, c_r;
  int obj1, obj2;
  double distance;
  bool groundContact, dynamicContact;
  Eigen::Vector3d contactPos1, contactPos2, Imp, prevVimp, Vimp;
  Eigen::Matrix3d Rc;
  Eigen::MatrixXd Jc1, Jc2, MappInv1, MappInv2;
};


class SimulationClass {
 public:
  SimulationClass() {
    dt = 0.001;
    g = 9.81;
    zero3d.setZero();
    e_z << 0,0,1;

    ball.resize(2);
    for (int bodyIdx = 0; bodyIdx < ball.size(); ++bodyIdx) {
      ball[bodyIdx].gc.setZero(7);
      ball[bodyIdx].gv.setZero(6);
      ball[bodyIdx].ga.setZero(6);
      ball[bodyIdx].gf.setZero(6);
      ball[bodyIdx].MassMatrix.setZero(6,6);
      ball[bodyIdx].Nonlinear.setZero(6);
    }

    ball[0].radius = 0.5;
    ball[1].radius = 0.7;

    ball[0].gc << 0, 0, 0.5, 1, 0, 0, 0;
    ball[1].gc << 0.1, 0.1, 3, 1, 0, 0, 0;

    double density = 1000;
    ball[0].mass = density * 4/3 * M_PI * std::pow(0.5,3);
    ball[1].mass = density * 4/3 * M_PI * std::pow(0.7,3);

    ball[0].b_I = 0.4 * ball[0].mass * std::pow(ball[0].radius, 2) * Eigen::Matrix3d::Identity();
    ball[1].b_I = 0.4 * ball[1].mass * std::pow(ball[1].radius, 2) * Eigen::Matrix3d::Identity();

    Contacts.push_back(Contact(0,-1));
    Contacts.push_back(Contact(1,-1));
    Contacts.push_back(Contact(0,1));

    Contacts[0].mu = 0.8;
    Contacts[1].mu = 0.8;
    Contacts[2].mu = 0.8;

    Contacts[0].c_r = 0.0;
    Contacts[1].c_r = 0.0;
    Contacts[2].c_r = 0.0;
  }

  void integrate() {
    updateDynamics();
    updateContacts();

    for (int bodyIdx = 0; bodyIdx < ball.size(); ++bodyIdx) {

      /// u dot
      ball[bodyIdx].ga = ball[bodyIdx].MassMatrix.inverse() * (ball[bodyIdx].gf - ball[bodyIdx].Nonlinear);

      /// contact effects
      auto it = objsInContact.find(bodyIdx);
      if (it != objsInContact.end()) {
        for (int Idx = 0; Idx < InContact.size(); ++Idx) {
          int Id = InContact[Idx];

          if (bodyIdx == Contacts[Id].obj1) {
            ball[bodyIdx].ga += ball[bodyIdx].MassMatrix.inverse() * (Contacts[Id].Jc1.transpose() * Contacts[Id].Rc * (Contacts[Id].Imp / dt));
          } else if (bodyIdx == Contacts[Id].obj2) {
            ball[bodyIdx].ga -= ball[bodyIdx].MassMatrix.inverse() * (Contacts[Id].Jc2.transpose() * Contacts[Id].Rc * (Contacts[Id].Imp / dt));
          }
        }
      }

      /// integrate
      ball[bodyIdx].gv += dt * ball[bodyIdx].ga;
      ball[bodyIdx].gc.head(3) += dt * ball[bodyIdx].gv.head(3);
      Eigen::Vector4d q = RotMatToQuaternion(QuaternionToRotMat(ball[bodyIdx].gc.tail(4)) * RotVecToRotMat(dt * ball[bodyIdx].gv.tail(3)));
      ball[bodyIdx].gc.tail(4) = q / q.norm();

    }
  }

  void setPosition(raisim::Visuals * sphere1, raisim::Visuals * sphere2) {
    sphere1->setPosition(ball[0].gc.head(3));
    sphere1->setOrientation(ball[0].gc.tail(4));
    sphere2->setPosition(ball[1].gc.head(3));
    sphere2->setOrientation(ball[1].gc.tail(4));
  }

  void getJ(const int &bodyIdx, const Eigen::Vector3d &pos, Eigen::MatrixXd &J) {
    J.setZero(6,6);
    J.block(0,0,3,3).setIdentity();
    J.block(3,3,3,3).setIdentity();
    J.block(0,3,3,3) = -skew(pos);
  }

  void getJdot(const int &bodyIdx, const Eigen::Vector3d &pos, Eigen::MatrixXd &dJ) {
    dJ.setZero(6,6);
    Eigen::MatrixXd J_target, J_center;
    getJ(bodyIdx, pos, J_target);
    getJ(bodyIdx, zero3d, J_center);
    dJ.block(0,3,3,3) = -skew(J_target.topRows(3) * ball[bodyIdx].gv - J_center.topRows(3) * ball[bodyIdx].gv);
  }

  void updateDynamics() {
    Eigen::MatrixXd J, dJ;
    for (int bodyIdx = 0; bodyIdx < ball.size(); ++bodyIdx) {

      /// Inertia in world frame
      ball[bodyIdx].w_I = QuaternionToRotMat(ball[bodyIdx].gc.tail(4)) * ball[bodyIdx].b_I * QuaternionToRotMat(ball[bodyIdx].gc.tail(4)).transpose();

      /// MassMatrix
      getJ(bodyIdx, zero3d, J);
//      ball[bodyIdx].MassMatrix = J.transpose() * ball[bodyIdx].mass * J;
      ball[bodyIdx].MassMatrix.topLeftCorner(3,3) = ball[bodyIdx].mass * Eigen::Matrix3d::Identity();
      ball[bodyIdx].MassMatrix.bottomRightCorner(3,3) = ball[bodyIdx].w_I;

      /// Nonlinear Term
//      getJdot(bodyIdx, zero3d, dJ);
//      ball[bodyIdx].Nonlinear = J.topRows(3).transpose() * ball[bodyIdx].mass * dJ.topRows(3) * ball[bodyIdx].gv + \
//                                  J.bottomRows(3).transpose() * ball[bodyIdx].w_I * dJ.bottomRows(3) * ball[bodyIdx].gv + \
//                                  J.bottomRows(3).transpose() * skew(J.bottomRows(3) * ball[bodyIdx].gv) * (ball[bodyIdx].w_I * J.bottomRows(3) * ball[bodyIdx].gv) + \
//                                  -J.topRows(3).transpose() * (ball[bodyIdx].mass * -g * e_z);
      ball[bodyIdx].Nonlinear.head(3).setZero();
      ball[bodyIdx].Nonlinear.tail(3) = skew(J.bottomRows(3) * ball[bodyIdx].gv) * (ball[bodyIdx].w_I * J.bottomRows(3) * ball[bodyIdx].gv);
      ball[bodyIdx].Nonlinear -= J.topRows(3).transpose() * (ball[bodyIdx].mass * -g * e_z);
    }
  }

  void updateContacts() {

    InContact.clear(); // InContact: vector having contact IDs that currently in contact.
    objsInContact.clear(); // objsInContact: set having all currently contacting objects except to the ground.

    for (int Id = 0; Id < Contacts.size(); ++Id) {

      /// contact kinematics
      Contacts[Id].update(ball);

      /// set Vimp- & reset Imp to 0
      Contacts[Id].reset(ball);

      if (Contacts[Id].isInCollision(ball)) {

        InContact.push_back(Id);
        objsInContact.insert(Contacts[Id].obj1);
        objsInContact.insert(Contacts[Id].obj2);
      }
    }
    objsInContact.erase(-1);

    PGSsolver();
  }

  void PGSsolver () {
    /// PGS solver
    ///  - calculate impulses at each contact points.
    ///  - All impulses are assumed to pointing inward to Obj1

    double c_z, c_t, alpha = 0.5, totalErr = 1.0;

    if (!InContact.empty()) {
      while (totalErr > 1e-8) {

        totalErr = 0.;
        for (int Idx = 0; Idx < InContact.size(); ++Idx) {
          int Id = InContact[Idx];

          /// compute Vimp+
          if (Contacts[Id].dynamicContact) {
            Contacts[Id].Vimp = Contacts[Id].prevVimp + Contacts[Id].Rc.transpose() * (Contacts[Id].MappInv1 + Contacts[Id].MappInv2) * Contacts[Id].Rc * Contacts[Id].Imp + \
                                dt * Contacts[Id].Rc.transpose() * (Contacts[Id].Jc1 * ball[Contacts[Id].obj1].MassMatrix.inverse() * (ball[Contacts[Id].obj1].gf - ball[Contacts[Id].obj1].Nonlinear) -
                Contacts[Id].Jc2 * ball[Contacts[Id].obj2].MassMatrix.inverse() * (ball[Contacts[Id].obj2].gf - ball[Contacts[Id].obj2].Nonlinear));
          } else {
            Contacts[Id].Vimp = Contacts[Id].prevVimp + Contacts[Id].Rc.transpose() * (Contacts[Id].MappInv1 + Contacts[Id].MappInv2) * Contacts[Id].Rc * Contacts[Id].Imp;
            Contacts[Id].Vimp += dt * Contacts[Id].Rc.transpose() * (Contacts[Id].Jc1 * ball[Contacts[Id].obj1].MassMatrix.inverse() * (ball[Contacts[Id].obj1].gf - ball[Contacts[Id].obj1].Nonlinear));
          }

          /// G terms
          for (int jit = 0; jit < InContact.size(); ++jit) {
            int j = InContact[jit];
            if (j == Id) { continue; }
            if (Contacts[j].obj1 == Contacts[Id].obj1) {
              Contacts[Id].Vimp += Contacts[Id].Rc.transpose() * Contacts[Id].Jc1 * ball[Contacts[Id].obj1].MassMatrix.inverse() * Contacts[j].Jc1.transpose() * Contacts[j].Rc * Contacts[j].Imp;
            } else if (Contacts[j].obj1 == Contacts[Id].obj2) {
              Contacts[Id].Vimp -= Contacts[Id].Rc.transpose() * Contacts[Id].Jc2 * ball[Contacts[Id].obj2].MassMatrix.inverse() * Contacts[j].Jc1.transpose() * Contacts[j].Rc * Contacts[j].Imp;
            } else if (Contacts[j].obj2 == Contacts[Id].obj1) {
              Contacts[Id].Vimp -= Contacts[Id].Rc.transpose() * Contacts[Id].Jc1 * ball[Contacts[Id].obj1].MassMatrix.inverse() * Contacts[j].Jc2.transpose() * Contacts[j].Rc * Contacts[j].Imp;
            } else if (Contacts[j].obj2 == Contacts[Id].obj2) {
              if(Contacts[Id].obj2 == -1) { continue; }
              Contacts[Id].Vimp += Contacts[Id].Rc.transpose() * Contacts[Id].Jc2 * ball[Contacts[Id].obj2].MassMatrix.inverse() * Contacts[j].Jc2.transpose() * Contacts[j].Rc * Contacts[j].Imp;
            }
          }

          /// update coefficients
          c_z = alpha / (Contacts[Id].Rc.transpose() * (Contacts[Id].MappInv1 + Contacts[Id].MappInv2) * Contacts[Id].Rc)(2,2);
          c_t = alpha / std::max((Contacts[Id].Rc.transpose() * (Contacts[Id].MappInv1 + Contacts[Id].MappInv2) * Contacts[Id].Rc)(0,0),
                                 (Contacts[Id].Rc.transpose() * (Contacts[Id].MappInv1 + Contacts[Id].MappInv2) * Contacts[Id].Rc)(1,1));

          /// calculate err
          totalErr += std::abs(Contacts[Id].Imp[2] - prox_z(Contacts[Id].Imp[2] - Contacts[Id].Vimp[2] * c_z));
          totalErr += (Contacts[Id].Imp.head(2) - prox_t(Contacts[Id].Imp.head(2) - Contacts[Id].Vimp.head(2) * c_t, Contacts[Id].Imp[2], Contacts[Id].mu)).norm();

          /// update Imps
          Contacts[Id].Imp[2] = prox_z(Contacts[Id].Imp[2] - Contacts[Id].Vimp[2] * c_z);
          Contacts[Id].Imp.head(2) = prox_t(Contacts[Id].Imp.head(2) - Contacts[Id].Vimp.head(2) * c_t, Contacts[Id].Imp[2], Contacts[Id].mu);

        }
      }
    }

  }

  double prox_z (const double Imp_z) {
    if (Imp_z < 0) {
      return 0;
    } else {
      return Imp_z;
    }
  }
  Eigen::Vector2d prox_t (const Eigen::Vector2d & in, const double & lam_z, const double & mu) {
    if (in.norm() <= mu * lam_z) {
      return in;
    } else {
      return ((mu * lam_z) / in.norm()) * in;
    }
  }

 private:
  double dt, g;
  Eigen::Vector3d zero3d, e_z;
  std::set<int> objsInContact;
  std::vector<int> InContact;
  std::vector<Ball> ball;
  std::vector<Contact> Contacts;
};
// This file is part of RaiSim. You must obtain a valid license from RaiSim Tech
// Inc. prior to usage.

#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>

int main(int argc, char* argv[]) {
  // fixed size matrices and vectors
  Eigen::Matrix3d A, B;
  Eigen::Vector3d x, y;

  // matrix initialization methods
  A.setZero();
  std::cout<<"A.setZero() \n"<<A<<"\n\n";

  A.setIdentity();
  std::cout<<"A.setIdentity() \n"<<A<<"\n\n";

  A.setConstant(2);
  std::cout<<"A.setConstant(2) \n"<<A<<"\n\n";

  A.setOnes();
  std::cout<<"A.setOnes() \n"<<A<<"\n\n";

  x.setLinSpaced(3,5,7);
  std::cout<<"x.setLinSpaced(3,5,7) \n"<<x<<"\n\n";

  A << x, x, x;
  std::cout<<"A << x, x, x \n"<<A<<"\n\n";

  A << 1, 2, 3,
       7, 4, 2,
       3, 5, 2;
  std::cout<<"A << 1, 2, 3,\n"
             "     7, 4, 2,\n"
             "     3, 5, 2 \n"<<A<<"\n\n";
  // eigen is column-major!

  // operations
  y = A * x;
  std::cout<<"y = A * x\n"<<y<<"\n\n";
  x = A.inverse() * y;
  std::cout<<"x = A.inverse() * y\n"<<x<<"\n\n";

  double xAx = x.transpose() * A * x;
  std::cout<<"x.transpose() * A * x\n"<<xAx<<"\n\n";

  // block operations
  A.col(0) = x;
  std::cout<<"A.col(0) = x\n"<<A<<"\n\n";

  // others
  std::cout<<"A is "<<A.rows()<<"X"<<A.cols()<<"\n\n";

}

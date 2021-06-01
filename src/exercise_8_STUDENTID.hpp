#pragma once

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrixUsingCRBA (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


    return Eigen::MatrixXd::Ones(9,9);
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearitiesUsingRNE (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return Eigen::VectorXd::Ones(9);
}


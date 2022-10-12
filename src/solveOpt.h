#ifndef _SOLVEOPT_H_
#define _SOLVEOPT_H_

#include <eigen/Eigen/Dense>

/// @brief Calculate optimized input vector
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam ec length of error vector
/// @tparam uc length of input vector
/// @param invH gradient matrix
/// @param F gradient matrix
/// @param xk states vector at current step
/// @param refg reference global vector for the next n steps
/// @return input vector
template <int n, int sc, int ec, int uc>
Eigen::Matrix<double, n * uc, 1> calculateOptInputs(Eigen::Matrix<double, n * uc, n * uc> &invH, Eigen::Matrix<double, n * uc, sc + n * ec> &F, Eigen::Matrix<double, sc, 1> &xk, Eigen::Matrix<double, ec * n, 1> &refg);

Eigen::MatrixXd calculateOptInputsdyn(Eigen::MatrixXd &invH, Eigen::MatrixXd &F, Eigen::MatrixXd &xk, Eigen::MatrixXd &refg, int n, int sc, int ec);

#endif
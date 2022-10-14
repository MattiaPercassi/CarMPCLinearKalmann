#ifndef _GRADIENTMATRIXES_H_
#define _GRADIENTMATRIXES_H_

#include <eigen/Eigen/Dense>

/// @brief Build H^-1 matrix
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam uc length of input vector
/// @param Bg Global Bg matrix
/// @param Qg Global Qg matrix
/// @param Rg Global Rg matrix
/// @return Matrix H^-1 of the cost gradient
template <int n, int sc, int uc>
Eigen::Matrix<double, n * uc, n * uc> buildinverseH(Eigen::Matrix<double, n * sc, n * uc> &Bg, Eigen::Matrix<double, n * sc, n * sc> &Qg, Eigen::Matrix<double, n * uc, n * uc> &Rg);

Eigen::MatrixXd buildinverseHdyn(Eigen::MatrixXd &Bg, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Rg);

/// @brief Build F matrix
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam ec length of error vector
/// @tparam uc length of input vector
/// @param Ag Global Ag matrix
/// @param Qg Global Qg matrix
/// @param Bg Global Bg matrix
/// @param Tg Global Tg matrix
/// @return Matrix F of the cost gradient
template <int n, int sc, int ec, int uc>
Eigen::Matrix<double, n * uc, sc + n * ec> buildF(Eigen::Matrix<double, n * sc, sc> &Ag, Eigen::Matrix<double, n * sc, n * sc> &Qg, Eigen::Matrix<double, n * sc, n * uc> &Bg, Eigen::Matrix<double, n * ec, n * sc> &Tg);

Eigen::MatrixXd buildFdyn(Eigen::MatrixXd &Ag, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Bg, Eigen::MatrixXd &Tg, int n, int sc, int ec, int uc);

#endif
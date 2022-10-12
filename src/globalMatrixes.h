#ifndef _GLOBALMATRIXES_H_
#define _GLOBALMATRIXES_H_

#include <eigen/Eigen/Dense>

// use a template function where the template parameter is the number of steps for the control

/// @brief Build global A matrix
/// @tparam n integrations steps
/// @tparam sc number of states
/// @param At Augmented matrix of the LTI discrete
/// @return Global Ag matrix
template <int n, int sc>
Eigen::Matrix<double, n * sc, sc> buildGlobalA(Eigen::Matrix<double, sc, sc> &At);

Eigen::MatrixXd buildGlobalAdyn(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &At, int n);

/// @brief Build global B matrix
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam uc number of inputs
/// @param At Augmented matrix of the LTI discrete
/// @param Bt Augmented matrix of the LTI discrete
/// @return Global Bg matrix
template <int n, int sc, int uc>
Eigen::Matrix<double, n * sc, n * uc> buildGlobalB(Eigen::Matrix<double, sc, sc> &At, Eigen::Matrix<double, sc, uc> &Bt);

Eigen::MatrixXd buildGlobalBdyn(Eigen::MatrixXd &At, Eigen::MatrixXd &Bt, int n, int sc, int uc);

/// @brief Build global Q matrix
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam ec length of error vector
/// @param Q Weight matrix for error vectors
/// @param S Weight matrix for last error vector
/// @param Ct Augmented matrix of the discrete LTI
/// @return Global Qg matrix
template <int n, int sc, int ec>
Eigen::Matrix<double, n * sc, n * sc> buildGlobalQ(Eigen::Matrix<double, ec, ec> &Q, Eigen::Matrix<double, ec, ec> &S, Eigen::Matrix<double, ec, sc> &Ct);

Eigen::MatrixXd buildGlobalQdyn(Eigen::MatrixXd &Q, Eigen::MatrixXd &S, Eigen::MatrixXd &Ct, int n, int sc, int ec);

/// @brief Build global T matrix
/// @tparam n integration steps
/// @tparam sc number of states
/// @tparam ec length of error vector
/// @param Q Weight matrix for error vectors
/// @param S Weight matrix for last error vector
/// @param Ct Augmented matrix of the discrete LTI
/// @return Global Tg matrix
template <int n, int sc, int ec>
Eigen::Matrix<double, n * ec, n * sc> buildGlobalT(Eigen::Matrix<double, ec, ec> &Q, Eigen::Matrix<double, ec, ec> &S, Eigen::Matrix<double, ec, sc> &Ct);

Eigen::MatrixXd buildGlobalTdyn(Eigen::MatrixXd &Q, Eigen::MatrixXd &S, Eigen::MatrixXd &Ct, int n, int sc, int ec);

/// @brief Build global R matrix
/// @tparam n integration steps
/// @tparam uc length of input vector
/// @param R Weight matrix for inputs
/// @return Global Rg matrix
template <int n, int uc>
Eigen::Matrix<double, n * uc, n * uc> buildGlobalR(Eigen::Matrix<double, uc, uc> &R);

Eigen::MatrixXd buildGlobalRdyn(Eigen::MatrixXd &R, int n, int uc);

#endif
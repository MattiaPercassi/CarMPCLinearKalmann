#ifndef _MPC_H_
#define _MPC_H_

#include <eigen/Eigen/Dense>

/// @brief Creates the At (tilde) matrix of the augmented system
/// @param Ad Matrix of the discrete LTI system
/// @param Bd Matrix of the discrete LTI system
/// @return Matrix At of the augmented discrete system
Eigen::MatrixXd augmentA(Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd);

/// @brief Creates the Bt matrix of the augmented system
/// @param Bd Matrix of the discrete LTI system
/// @return Matrix Bt of the augmented discrete system
Eigen::MatrixXd augmentB(Eigen::MatrixXd &Bd);

/// @brief Creates the Ct output matrix of the augmented system
/// @param Cd Matrix of the discrete LTI system
/// @return Matrix Ct of the augmented discrete system
Eigen::MatrixXd augmentC(Eigen::MatrixXd &Cd);

/// @brief Build the discretized Ad matrix from the system state space A matrix
/// @param A State space matrix of the continuous system (LTI)
/// @param h integration interval
/// @return Ad matrix of the disctete LTI
Eigen::MatrixXd discretizeA(Eigen::MatrixXd &A, double h);

/// @brief Build the discretized Bd matrix from the system state space B matrix
/// @param B State space matrix B of the continuous system (LTI)
/// @param h integration interval
/// @return Bd matrix of the discrete LTI
Eigen::MatrixXd discretizeB(Eigen::MatrixXd &B, double h);

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
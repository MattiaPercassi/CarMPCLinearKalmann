#include <eigen/Eigen/Dense>
#include <cmath>
#include "openLoop.h"

/// @brief computes new states of open loop system
/// @param xk current state vector
/// @param Ad discrete A matrix
/// @param Bd discrete A matrix
/// @param delta control input
/// @return new states vector
Eigen::MatrixXd updateOLS(Eigen::MatrixXd &xk, Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd, double delta)
{
    xk = Ad * xk + Bd * delta;
    return xk;
};
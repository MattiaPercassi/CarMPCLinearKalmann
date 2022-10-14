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

Eigen::MatrixXd updateOLSnonLinear(Eigen::MatrixXd &xk, double xd, double dt, double delta, Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd)
{
    double psidot = Ad(0, 0) * xk(0, 0) + Ad(0, 1) * xk(1, 0) + Bd(0, 0) * delta;
    double ydot = Ad(1, 0) * xk(0, 0) + Ad(1, 1) * xk(1, 0) + Bd(1, 0) * delta;
    double psi = xk(2, 0) + Ad(2, 0) * xk(0, 0);
    double Y = xk(3, 0) + (xd * sin(xk(2, 0)) + xk(1, 0) * cos(xk(2, 0))) * dt;
    xk(0, 0) = psidot;
    xk(1, 0) = ydot;
    xk(2, 0) = psi;
    xk(3, 0) = Y;
    return xk;
};
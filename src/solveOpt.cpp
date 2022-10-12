#include <eigen/Eigen/Dense>
#include "solveOpt.h"

template <int n, int sc, int ec, int uc>
Eigen::Matrix<double, n * uc, 1> calculateOptInputs(Eigen::Matrix<double, n * uc, n * uc> &invH, Eigen::Matrix<double, n * uc, sc + n * ec> &F, Eigen::Matrix<double, sc, 1> &xk, Eigen::Matrix<double, ec * n, 1> &refg)
{
    Eigen::Matrix<double, sc + ec * n, 1> xtemp = Eigen::Matrix<double, sc + ec * n, 1>::Zero();
    xtemp.block(0, 0, sc, 1) = xk;
    xtemp.block(sc, 0, n * ec, 1) = refg;
    Eigen::Matrix<double, n * uc, 1> ug = -invH * F * xtemp.transpose();
    return ug;
};

Eigen::MatrixXd calculateOptInputsdyn(Eigen::MatrixXd &invH, Eigen::MatrixXd &F, Eigen::MatrixXd &xk, Eigen::MatrixXd &refg, int n, int sc, int ec)
{
    Eigen::MatrixXd xtemp = Eigen::MatrixXd::Zero(sc + ec * n, 1);
    xtemp.block(0, 0, sc, 1) = xk;
    xtemp.block(sc, 0, n * ec, 1) = refg;
    auto ug = -invH * F * xtemp;
    return ug;
};
#include <eigen/Eigen/Dense>
#include "augmentLTI.h"

Eigen::MatrixXd augmentA(Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd)
{
    Eigen::MatrixXd At = Eigen::Matrix<double, 5, 5>::Zero();
    At.block(0, 0, 4, 4) = Ad;
    At.block(0, 4, 4, 1) = Bd;
    At.block(4, 0, 1, 4) = Eigen::Matrix<double, 1, 4>::Zero();
    At.block(4, 4, 1, 1) = Eigen::Matrix<double, 1, 1>::Identity(); // or identity matrix
    return At;
};

Eigen::MatrixXd augmentB(Eigen::MatrixXd &Bd)
{
    Eigen::MatrixXd Bt = Eigen::Matrix<double, 5, 1>::Zero();
    Bt.block(0, 0, 4, 1) = Bd;
    Bt.block(4, 0, 1, 1) = Eigen::Matrix<double, 1, 1>::Identity(); // or identity matrix
    return Bt;
};

Eigen::MatrixXd augmentC(Eigen::MatrixXd &Cd)
{
    Eigen::MatrixXd Ct = Eigen::Matrix<double, 2, 5>::Zero();
    Ct.block(0, 0, 2, 4) = Cd;
    Ct.block(0, 4, 2, 1) = Eigen::Matrix<double, 2, 1>::Zero();
    return Ct;
}
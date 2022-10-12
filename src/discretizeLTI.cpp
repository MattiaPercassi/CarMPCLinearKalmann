#include <eigen/Eigen/Dense>
#include "discretizeLTI.h"

Eigen::MatrixXd discretizeA(Eigen::MatrixXd &A, double h)
{
    Eigen::MatrixXd Ad = Eigen::Matrix<double, 4, 4>::Identity();
    Ad += A * h;
    return Ad;
};

Eigen::MatrixXd discretizeB(Eigen::MatrixXd &B, double h)
{
    Eigen::MatrixXd Bd = B * h;
    return Bd;
};
#include <eigen/Eigen/Dense>
#include "gradientMatrixes.h"

template <int n, int sc, int uc>
Eigen::Matrix<double, n * uc, n * uc> buildinverseH(Eigen::Matrix<double, n * sc, n * uc> &Bg, Eigen::Matrix<double, n * sc, n * sc> &Qg, Eigen::Matrix<double, n * uc, n * uc> &Rg)
{
    Eigen::Matrix<double, n * uc, n *uc> H = Bg.transpose() * Qg * Bg + Rg;
    return H.inverse();
};

Eigen::MatrixXd buildinverseHdyn(Eigen::MatrixXd &Bg, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Rg)
{
    auto H = Bg.transpose() * Qg * Bg + Rg;
    return H.inverse();
};

template <int n, int sc, int ec, int uc>
Eigen::Matrix<double, n * uc, sc + n * ec> buildF(Eigen::Matrix<double, n * sc, sc> &Ag, Eigen::Matrix<double, n * sc, n * sc> &Qg, Eigen::Matrix<double, n * sc, n * uc> &Bg, Eigen::Matrix<double, n * ec, n * sc> &Tg)
{
    Eigen::Matrix<double, sc + n * ec, n *uc> Ft = Eigen::Matrix<double, sc + n * ec, n * uc>::Zero();
    Ft.block(0, 0, sc, n * uc) = Ag.transpose() * Qg * Bg;
    Ft.block(sc, 0, n * ec, n * uc) = -Tg * Bg;
    return Ft.transpose();
};

Eigen::MatrixXd buildFdyn(Eigen::MatrixXd &Ag, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Bg, Eigen::MatrixXd &Tg, int n, int sc, int ec, int uc)
{
    Eigen::MatrixXd Ft = Eigen::MatrixXd::Zero(sc + n * ec, n * uc);
    Ft.block(0, 0, sc, n * uc) = Ag.transpose() * Qg * Bg;
    Ft.block(sc, 0, n * ec, n * uc) = -Tg * Bg;
    return Ft.transpose();
};
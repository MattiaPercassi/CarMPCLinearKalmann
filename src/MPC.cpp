#include "MPC.h"
#include <eigen/Eigen/Dense>
#include <iostream>

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
};

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

Eigen::MatrixXd buildGlobalAdyn(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &At, const int n)
{
    const int sc = At.rows();
    Eigen::MatrixXd Ag = Eigen::MatrixXd::Zero(n * sc, sc);
    for (int i{0}; i < n; ++i)
    {
        auto temp = At;
        for (int j{0}; j < i; ++j)
        {
            temp *= At;
        }
        Ag.block(i * sc, 0, sc, sc) = temp;
    }
    return Ag;
};

Eigen::MatrixXd buildGlobalBdyn(Eigen::MatrixXd &At, Eigen::MatrixXd &Bt, int n, int sc, int uc)
{
    Eigen::MatrixXd Bg = Eigen::MatrixXd::Zero(n * sc, n * uc);
    Eigen::MatrixXd Bgcol = Eigen::MatrixXd::Zero(n * sc, uc);
    for (int i{0}; i < n; ++i)
    {
        auto temp = Bt;
        for (int j{0}; j < i; ++j)
        {
            temp = At * temp;
        };
        Bgcol.block(i * sc, 0, sc, uc) = temp;
    };

    // build each column of the Bg matrix
    for (int i{0}; i < n; ++i)
    {
        Bg.block(i * sc, i * uc, (n - i) * sc, uc) = Bgcol.block(0, 0, (n - i) * sc, uc);
    };
    return Bg;
};

Eigen::MatrixXd buildGlobalQdyn(Eigen::MatrixXd &Q, Eigen::MatrixXd &S, Eigen::MatrixXd &Ct, int n, int sc, int ec)
{
    Eigen::MatrixXd Qg = Eigen::MatrixXd::Zero(n * sc, n * sc);
    auto temp = Ct.transpose() * Q * Ct;
    for (int i{0}; i < n; ++i)
    {
        if (i == n - 1)
            Qg.block(i * sc, i * sc, sc, sc) = Ct.transpose() * S * Ct;
        else
        {
            Qg.block(i * sc, i * sc, sc, sc) = temp;
        }
    }
    return Qg;
};

Eigen::MatrixXd buildGlobalTdyn(Eigen::MatrixXd &Q, Eigen::MatrixXd &S, Eigen::MatrixXd &Ct, int n, int sc, int ec)
{
    Eigen::MatrixXd Tg = Eigen::MatrixXd::Zero(n * ec, n * sc);
    auto temp = Q * Ct;
    for (int i{0}; i < n; ++i)
    {
        if (i == n - 1)
            Tg.block(i * ec, i * sc, ec, sc) = S * Ct;
        else
        {
            Tg.block(i * ec, i * sc, ec, sc) = temp;
        }
    }
    return Tg;
};

Eigen::MatrixXd buildGlobalRdyn(Eigen::MatrixXd &R, int n, int uc)
{
    Eigen::MatrixXd Rg = Eigen::MatrixXd::Zero(n * uc, n * uc);
    for (int i{0}; i < n; ++i)
    {
        Rg.block(i * uc, i * uc, uc, uc) = R;
    }
    return Rg;
};

Eigen::MatrixXd buildinverseHdyn(Eigen::MatrixXd &Bg, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Rg)
{
    auto H = Bg.transpose() * Qg * Bg + Rg;
    return H.inverse();
};

Eigen::MatrixXd buildFdyn(Eigen::MatrixXd &Ag, Eigen::MatrixXd &Qg, Eigen::MatrixXd &Bg, Eigen::MatrixXd &Tg, int n, int sc, int ec, int uc)
{
    Eigen::MatrixXd Ft = Eigen::MatrixXd::Zero(sc + n * ec, n * uc);
    Ft.block(0, 0, sc, n * uc) = Ag.transpose() * Qg * Bg;
    Ft.block(sc, 0, n * ec, n * uc) = -Tg * Bg;
    return Ft.transpose();
};

Eigen::MatrixXd calculateOptInputsdyn(Eigen::MatrixXd &invH, Eigen::MatrixXd &F, Eigen::MatrixXd &xk, Eigen::MatrixXd &refg, int n, int sc, int ec)
{
    Eigen::MatrixXd xtemp = Eigen::MatrixXd::Zero(sc + ec * n, 1);
    xtemp.block(0, 0, sc, 1) = xk;
    xtemp.block(sc, 0, n * ec, 1) = refg;
    auto ug = -invH * F * xtemp;
    return ug;
};
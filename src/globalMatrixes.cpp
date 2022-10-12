#include <eigen/Eigen/Dense>
#include "globalMatrixes.h"
#include <iostream>

template <int n, int sc>
Eigen::Matrix<double, n * sc, sc> buildGlobalA(Eigen::Matrix<double, sc, sc> &At)
{
    Eigen::Matrix<double, n * sc, sc> Ag = Eigen::Matrix<double, n * sc, sc>::Zero();
    for (int i{0}; i < n; ++i)
    {
        Eigen::Matrix<double, sc, sc> temp = At;
        for (int j{0}; j < i; ++j)
        {
            temp *= temp;
        };
        Ag.block(i * sc, 0, sc, sc) = temp;
    }
    return Ag;
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

template <int n, int sc, int uc>
Eigen::Matrix<double, n * sc, n * uc> buildGlobalB(Eigen::Matrix<double, sc, sc> &At, Eigen::Matrix<double, sc, uc> &Bt)
{
    Eigen::Matrix<double, n * sc, n *uc> Bg = Eigen::Matrix<double, n * sc, n * uc>::Zero();
    // build longest matrix of the "first column" of Bg
    Eigen::Matrix<double, n * sc, uc> Bgcol;
    for (int i{0}; i < n; ++i)
    {
        Eigen::Matrix<double, sc, uc> temp = Bt;
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

template <int n, int sc, int ec>
Eigen::Matrix<double, n * sc, n * sc> buildGlobalQ(Eigen::Matrix<double, ec, ec> &Q, Eigen::Matrix<double, ec, ec> &S, Eigen::Matrix<double, ec, sc> &Ct)
{
    Eigen::Matrix<double, n * sc, n *sc> Qg = Eigen::Matrix<double, n * sc, n * sc>::Zero();
    Eigen::Matrix<double, sc, sc> temp = Ct.transpose() * Q * Ct;
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

template <int n, int sc, int ec>
Eigen::Matrix<double, n * ec, n * sc> buildGlobalT(Eigen::Matrix<double, ec, ec> &Q, Eigen::Matrix<double, ec, ec> &S, Eigen::Matrix<double, ec, sc> &Ct)
{
    Eigen::Matrix<double, n * ec, n *sc> Tg = Eigen::Matrix<double, n * ec, n * sc>::Zero();
    Eigen::Matrix<double, ec, sc> temp = Q * Ct;
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

template <int n, int uc>
Eigen::Matrix<double, n * uc, n * uc> buildGlobalR(Eigen::Matrix<double, uc, uc> &R)
{
    Eigen::Matrix<double, n * uc, n *uc> Rg = Eigen::Matrix<double, n * uc, n * uc>::Zero();
    for (int i{0}; i < n; ++i)
    {
        Rg.block(i * uc, i * uc, uc, uc) = R;
    }
    return Rg;
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
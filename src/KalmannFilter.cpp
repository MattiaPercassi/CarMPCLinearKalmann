#include <eigen/Eigen/Dense>
#include "KalmannFilter.h"

KalmannFilter::KalmannFilter(Eigen::MatrixXd h, Eigen::MatrixXd r) : H{h}, R{r} {};

int KalmannFilter::calcInn(Eigen::MatrixXd &z, Eigen::MatrixXd &x)
{
    // TODO - add exception handling
    innovation = z - H * x;
    return 0;
};

int KalmannFilter::calcInnCov(Eigen::MatrixXd &P)
{
    // TODO - add exception handling
    S = H * P * H.transpose() + R;
    return 0;
};

int KalmannFilter::calcTuneMatrix(Eigen::MatrixXd &P)
{
    K = P * H.transpose() * S.inverse();
    return 0;
};

Eigen::MatrixXd KalmannFilter::updateEst(Eigen::MatrixXd &x)
{
    Eigen::MatrixXd xplus = x + K * innovation;
    return xplus;
};

Eigen::MatrixXd KalmannFilter::updateEstCov(Eigen::MatrixXd &P)
{
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(P.rows(), P.cols());
    Eigen::MatrixXd Pplus = (I - K * H) * P;
    return Pplus;
};

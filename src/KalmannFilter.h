#ifndef _KALMANNFILTER_H_
#define _KALMANNFILTER_H_

#include <eigen/Eigen/Dense>

class KalmannFilter
{
protected:
    Eigen::MatrixXd H; // matrix defining relationship between measurement and states
    Eigen::MatrixXd R; // covariance matrix (or scalar) of the measurement
    Eigen::MatrixXd innovation;
    Eigen::MatrixXd S;
    Eigen::MatrixXd K;

public:
    /// @brief initialize kalmann filter matrixes
    /// @param h states to measurement matrix
    /// @param r measurement covariance matrix
    KalmannFilter(Eigen::MatrixXd h, Eigen::MatrixXd r);
    /// @brief Calculate Innovation
    /// @param z measurement vector
    /// @param x a priori estimate state vector
    /// @return !=0 for error
    int calcInn(Eigen::MatrixXd &z, Eigen::MatrixXd &x);
    /// @brief Calculate innovation covariance
    /// @param P a priori state covariance matrix
    /// @return !=0 for error
    int calcInnCov(Eigen::MatrixXd &P);
    /// @brief Calculate tune matrix
    /// @param P a priori state covariance matrix
    /// @return !=0 for error
    int calcTuneMatrix(Eigen::MatrixXd &P);
    /// @brief update states estimate
    /// @param x a priori estimate state vector
    /// @return a posteriori estimate state vector
    Eigen::MatrixXd updateEst(Eigen::MatrixXd &x);
    /// @brief update states estimete covariance
    /// @param P a priori state covariance matrix
    /// @return a posteriori state covariance matrix
    Eigen::MatrixXd updateEstCov(Eigen::MatrixXd &P);
};

#endif
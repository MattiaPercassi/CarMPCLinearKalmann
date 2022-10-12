#ifndef _TRAJECTORYGEN_H_
#define _TRAJECTORYGEN_H_

#include <eigen/Eigen/Dense>

/// @brief generate tanh-like trajectory
/// @param steps number of simulation steps
/// @param ec number of references
/// @param xd longitudinal vehicle speed
/// @param h time interval for sample step
/// @param Ys initial value of Y
/// @param Ye final value of Y
/// @return reference trajectory for psi and Y
Eigen::MatrixXd tanhTraj(double steps, int ec, double xd, double h, double Ys, double Ye);

/// @brief generate sin-like trajectory
/// @param steps number of simulation steps
/// @param ec number of references
/// @param xd longitudinal vehicle speed
/// @param h time interval for sample step
/// @param a amplitude
/// @return reference trajectory for psi and Y
Eigen::MatrixXd sinTraj(double steps, int ec, double xd, double h, double a);

#endif
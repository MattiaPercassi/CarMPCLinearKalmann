#ifndef _OPENLOOP_H_
#define _OPENLOOP_H_

#include <eigen/Eigen/Dense>

Eigen::MatrixXd updateOLS(Eigen::MatrixXd &xk, double delta);

Eigen::MatrixXd updateOLSnonLinear(Eigen::MatrixXd &xk, double xd, double dt, double delta, Eigen::MatrixXd &Ad, Eigen::MatrixXd &Bd);

#endif
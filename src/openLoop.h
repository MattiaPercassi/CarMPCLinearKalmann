#ifndef _OPENLOOP_H_
#define _OPENLOOP_H_

#include <eigen/Eigen/Dense>

Eigen::MatrixXd updateOLS(Eigen::MatrixXd &xk, double delta);

#endif
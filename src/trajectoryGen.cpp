#include <cmath>
#include <eigen/Eigen/Dense>
#include "trajectoryGen.h"
Eigen::MatrixXd tanhTraj(double steps, int ec, double xd, double h, double Ys, double Ye)
{
    Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(steps * ec, 1);
    // y reference
    for (int i{0}; i < steps; ++i)
    {
        ref(i * ec + 1, 0) = std::tanh(h * (i - steps / 2)) * (Ye - Ys) / 2 + (Ye + Ys) / 2;
    };
    // psi reference
    double dx{h * xd}; // considered constant in approx of small psi
    double dy{0};
    for (int i{1}; i < steps; ++i)
    {
        dy = ref(i * ec + 1, 0) - ref(i * ec - 1, 0);
        ref(i * ec, 0) = std::atan(dy / dx);
    };
    return ref;
};

Eigen::MatrixXd sinTraj(double steps, int ec, double xd, double h, double a)
{
    Eigen::MatrixXd ref = Eigen::MatrixXd::Zero(steps * ec, 1);
    for (int i{0}; i < steps; ++i)
    {
        ref(i * ec + 1, 0) = std::sin(h * (i - steps / 2)) * a + a;
    };
    // psi reference
    double dx{h * xd}; // considered constant in approx of small psi
    double dy{0};
    for (int i{1}; i < steps; ++i)
    {
        dy = ref(i * ec + 1, 0) - ref(i * ec - 1, 0);
        ref(i * ec, 0) = std::atan(dy / dx);
    };
    return ref;
}
#include <eigen/Eigen/Dense>
#include <iostream>
#include <fstream>
#include "discretizeLTI.h"
#include "augmentLTI.h"
#include "globalMatrixes.h"
#include "gradientMatrixes.h"
#include "solveOpt.h"
#include "openLoop.h"
#include "trajectoryGen.h"
#include "KalmannFilter.h"
#include <ctime>
#include <random>
#include <cstdlib>

int main()
{
    // setup random generator for the creation of data
    srand(time(NULL));
    std::seed_seq sd{rand()};
    std::mt19937 rng{sd};

    // pi
    constexpr double pi = 3.1415;
    // For simplicity matrixes dimentions are defined in initialization
    constexpr int n = 20; // prediction steps
    constexpr int ec = 2; // error vector size
    constexpr int uc = 1; // number of inputs
    constexpr int sc = 5; // number of states
    // constexpr int frac = 20;  // number of fractions for fine status calculation at each step
    constexpr double h = 0.02; // [s] integration interval, 1/sample frequency
    // constexpr double smallh = h / frac; // [s] integration interval of the open loop system simulation
    constexpr int simTime = 10;   // [s] total time of simulation
    constexpr double m = 1500;    //[kg] vehicle mass
    constexpr double j = 3000;    //[kg*m^2] moment if inertia
    constexpr double caf = 19000; // [N] front wheel lateral force coefficient
    constexpr double car = 33000; // [N] rear wheel lateral force coefficient
    constexpr double lf = 2;      // [m] distance of center of mass and front wheel
    constexpr double lr = 3;      // [m] distance of center of mass and rear wheel
    constexpr double xd = 20;     // [m/s] longitudianl velocity
    // initial states vector : psidot, ydot, psi, Y, u(k-1)
    Eigen::MatrixXd xk = Eigen::MatrixXd::Zero(sc, 1);
    // vector of the system outputs
    Eigen::MatrixXd yk = Eigen::MatrixXd::Zero(ec, 1);
    // initialization input vector
    Eigen::MatrixXd ug = Eigen::MatrixXd::Zero(n * uc, 1);
    // initialization reference, reference vector is extended after hte ned of the simulation time samples for allowing definiton of refg until the last simulation step
    //[TANH]
    // Eigen::MatrixXd ref = tanhTraj(simTime / h + n * sc, ec, xd, h, 2, 10);
    //[SIN]
    Eigen::MatrixXd ref = sinTraj(simTime / h + n * sc, ec, xd, h, 3);
    Eigen::MatrixXd refg = Eigen::MatrixXd::Zero(n * ec, 1);

    // matrix constants
    double a1 = -2 / (j * xd) * (caf * lf * lf + car * lr * lr);
    double a2 = 2 / (j * xd) * (-caf * lf + car * lr);
    double a3 = -xd - (2 / (m * xd) * (caf * lf - car * lr));
    double a4 = -2 / (m * xd) * (caf + car);
    double b1 = 2 * caf * lf / j;
    double b2 = 2 * caf / m;
    // build matrix A
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 4);
    A << a1, a2, 0, 0, a3, a4, 0, 0, 1, 0, 0, 0, 0, 1, xd, 0;
    // build matrix B
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(4, 1);
    B << b1, b2, 0, 0;
    // build matrix C
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(2, 4);
    C << 0, 0, 1, 0, 0, 0, 0, 1;
    // Build matrix Q
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(2, 2);
    Q << 100, 0, 0, 1;
    // Build matrix S
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(2, 2);
    S << 300, 0, 0, 30;
    // Build matrix R
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(uc, uc);
    R << 200;

    // discretize
    auto Ad = discretizeA(A, h);
    auto Bd = discretizeB(B, h);
    auto Cd = C;

    // augment
    auto At = augmentA(Ad, Bd);
    auto Bt = augmentB(Bd);
    auto Ct = augmentC(Cd);

    // global matrixes
    auto Ag = buildGlobalAdyn(At, n);
    auto Bg = buildGlobalBdyn(At, Bt, n, sc, uc);
    auto Qg = buildGlobalQdyn(Q, S, Ct, n, sc, ec);
    auto Tg = buildGlobalTdyn(Q, S, Ct, n, sc, ec);
    auto Rg = buildGlobalRdyn(R, n, uc);

    // gradient matrixes
    auto invH = buildinverseHdyn(Bg, Qg, Rg);
    auto F = buildFdyn(Ag, Qg, Bg, Tg, n, sc, ec, uc);

    // Kalmann filter initialization, the filter will take in input the discretized matrix of the non augmented system as other computation is not required
    Eigen::MatrixXd Hgyro = Eigen::MatrixXd::Zero(1, 4);
    Hgyro << 0, 0, 1, 0;
    Eigen::MatrixXd Rgyro = Eigen::MatrixXd::Zero(1, 1);
    constexpr double gyrovar = pi / 60; // arbitrary for exercise purposes
    Rgyro << gyrovar;
    std::normal_distribution<> Errgyro{0, gyrovar};
    KalmannFilter gyroF{Hgyro, Rgyro};

    // initializing covariance matrix of the 4 states (excluded u(k-1))
    Eigen::MatrixXd stateestimate = Eigen::MatrixXd::Zero(4, 1);
    Eigen::MatrixXd Pstatecov = Eigen::MatrixXd::Zero(4, 4);
    Pstatecov << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pi / 60, 0, 0, 0, 0, 0; // arbitrary for exercises purposes

    // steering angle
    double delta{0};
    std::ofstream file{"simRes.txt"};
    if (!file)
    {
        std::cout << "Could not open file" << std::endl;
        return 1;
    }
    for (int i{0}; i < simTime / h; ++i)
    {
        refg = ref.block(i * ec, 0, n * ec, 1);
        Eigen::MatrixXd error = ref.block(i * ec, 0, ec, 1) - stateestimate.block(2, 0, 2, 1);
        ug = calculateOptInputsdyn(invH, F, xk, refg, n, sc, ec);
        delta += ug(0, 0);
        // delta must be constrained between the real available steering angles +- pi/6
        if (delta < -pi / 6)
            delta = -pi / 6;
        else if (delta > pi / 6)
            delta = pi / 6;

        xk(4, 0) = delta;

        // For test purposes update the states with same refinement as controller sample frequency, this represents both the system and the system state prediction
        xk.block(0, 0, 4, 1) = Ad * xk.block(0, 0, 4, 1) + Bd * delta;
        yk = Ct * xk;
        stateestimate = xk.block(0, 0, 4, 1);

        // create a measurement for the gyro
        Eigen::MatrixXd z = Eigen::MatrixXd::Zero(1, 1);
        z << xk(2, 0) + Errgyro(rng);
        gyroF.calcInn(z, stateestimate);
        gyroF.calcInnCov(Pstatecov);
        gyroF.calcTuneMatrix(Pstatecov);
        stateestimate = gyroF.updateEst(stateestimate);
        Pstatecov = gyroF.updateEstCov(Pstatecov);

        // write to file file
        file << yk(0, 0) << "," << refg(0, 0) << "," << stateestimate(2, 0) << "," << yk(1, 0) << "," << refg(1, 0) << "," << stateestimate(3, 0) << "," << i * h << std::endl;
    };
    file.close();
    return 0;
};
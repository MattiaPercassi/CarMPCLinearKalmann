Car lateral MPC + Linear Kalmann Filter

This project is meant to be an implementation of a simulation of an MPC controller on a car following a trajectory in the XY plane and the application of a Kalmann filter on simulated measurements.

-PHYSICAL MODEL-
The model used is a single rigid body bycicle model.
The system can be described in terms of 3 dof, the absolute position X, the absolute position Y and the yat angle PSI.
For simplification, the equations of motions are written with respect to the local referece xy that is attached to the center of mass of the moving vehicle.
The control action on the system is the angle delta, which is the steering angle of the car, through the variation of which we can control the modulus of the velocity vector along the y local axis.
The longitudinal velocity of the car is assumed to be constant xd, and due to this semplification no forces will be considered along the x direction.
Due to the previous assumption, the component of position X and Y are linked together from the velocity in the local y direction and the yaw angle PSI, therefore we will exclude the variable X in our state vector considering the following:
Yd = xd * sin(PSI) + yd * cos(PSI)
Xd = xd * cos(PSI) - yd * sin (PSI)
From which it is evident that PSI and yd are the only two variables

The model used is a LTI system fed into a basic MPC controller
Finally the discrete LTI model is augmented to take as input the Ddelta (variation of the steering angle) as it is assumed the minimization will not be done on the steering angle alone (this would produce low steering and potentially vibrating around 0), but instead on the variation of the steering.

-MATHEMATICAL DETAILS-
All integration of differential equations is done through the forward Euler method.
The transformation of the system equations in state space system is meant to rewrite the problem in term of ODE equations only.

-KALMANN FILTER-
A linear Kalmann filter is used.
The measurement of the gyroscope is obtained through adding a 0 mean error as a standard distribution on the current system state. Due to the use of the same LTI model for both the MPC and for updating the states of the simulation the error is small.
However as a proof of concept we can observe the error decreasing over time, as expected from linear Kalmann filter properies.
The filter is applied on the basic discrete system with 4 states while the augmented system is used for the MPC code.

-SIMULATION-
To solve linear algebra the library Eigen is used (https://eigen.tuxfamily.org/index.php?title=Main_Page).
The time of the simulation is set to a number of secods, based on which the trajectory vector will be generated and with it the main loop will proceed.
The simulation will proceed along the following logic

- Initialize all quantities, matrixes and trajectory vector

- Main simulation loop (one loop per each MPC update step, the update step is a variable to be defined)
    - Select proper reference vector elements for the time step
    - Calculate error between the reference vector and the state estimation vector
    - Calculate the optimized Ddelta vectors for the next N time steps
    - Calculate steering angle (real system input) with the first element of the Ddelta vector (steering angle will be limited within an interval to represent the physical limit of a real steering system)
    - Calculation of the a priori system state estimation through the LTI model
    - System status update loop (one loop per each fraction of update step, update step to be divided in M fractions, value to be defined) [CURRENTLY NOT YET IMPLEMENTED]
        - use the input Ddelta and current status to evaluate the system with the non-linear state space equations
    - Creation of the gyro measurement and evaluation of the a posteriori system estimate
    - Write the states of the system in a convenient location for data plotting

-DATA VISUALIZATION-
A simple python script using matplotlib will be used to plot the various states with respect to time and analyze the results

Here below results obtained with a tanh-like and sin-like trajectory

TANH

![XYt](XYplaneTanh.png)
![Yawt](yawAngleSinTanh.png)
![errort](errorsTanh.png)

SIN

![XYs](XYplaneSin.png)
![Yaws](yawAngleSin.png)
![errors](errorsSin.png)
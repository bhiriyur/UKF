#include <iostream>
#include "ukf.h"

UKF::UKF() {
    //TODO Auto-generated constructor stub
    Init();
}

UKF::~UKF() {
    //TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

    //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //create example sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    Xsig_aug <<
             5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
            1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
            2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
            0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
            0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
            0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
            0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

    double delta_t = 0.1; //time diff in sec
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

    for (int i = 0; i < 2*n_aug+1; ++i) {
        //predict sigma points
        VectorXd Xi = Xsig_aug.col(i);
        double px = Xi(0);
        double py = Xi(1);
        double v = Xi(2);
        double yaw = Xi(3);
        double dyawdt = Xi(4);
        double noise_v = Xi(5);
        double noise_yaw = Xi(6);

        VectorXd Xpred_i(5);
        if (std::fabs(dyawdt)<0.0001) {
            Xpred_i <<  px     + v*cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*noise_v,
                        py     + v*sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*noise_v,
                        v      + 0                  + delta_t*noise_v,
                        yaw    + dyawdt*delta_t     + 0.5*delta_t*delta_t*noise_yaw,
                        dyawdt + 0                  + delta_t*noise_yaw;
        } else {
            Xpred_i <<  px     + v/dyawdt*( sin(yaw+dyawdt*delta_t)-sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*noise_v,
                        py     + v/dyawdt*(-cos(yaw+dyawdt*delta_t)+cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*noise_v,
                        v      + 0                                            + delta_t*noise_v,
                        yaw    + dyawdt*delta_t                               + 0.5*delta_t*delta_t*noise_yaw,
                        dyawdt + 0                                            + delta_t*noise_yaw;

        }
        Xsig_pred.col(i) = Xpred_i;
    }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

    //print result
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

    //write result
    *Xsig_out = Xsig_pred;

}

#include <iostream>
#include "Eigen/Dense"
// #include <vector>
#include "ukf.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

int main() {

    //Create a UKF instance
    UKF ukf;

/*******************************************************************************
* Programming assignment calls
*******************************************************************************/

    MatrixXd Xsig;
    ukf.GenerateSigmaPoints(&Xsig);

    MatrixXd Xaug;
    ukf.AugmentedSigmaPoints(&Xaug);

    MatrixXd Xsig_pred;
    ukf.SigmaPointPrediction(&Xsig_pred);

    VectorXd x_pred;
    MatrixXd P_pred;
    ukf.PredictMeanAndCovariance(&x_pred, &P_pred);

    VectorXd z_meas;
    MatrixXd S_meas;
    ukf.PredictRadarMeasurement(&z_meas, &S_meas);

    VectorXd x_new;
    MatrixXd P_new;
    ukf.UpdateState(&x_new, &P_new);

    return 0;
}
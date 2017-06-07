#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  n_x_ = 5;
  n_aug_ = 7;
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  //cout << "PROCESS MEASUREMENT BEGIN" << endl;

  if (!is_initialized_) {

    P_ << 1,   0,    0,    0,    0,
          0,   1,    0,    0,    0,
          0,   0,    1,    0,    0,
          0,   0,    0,    1,    0,
          0,   0,    0,    0,    1;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho0 = meas_package.raw_measurements_[0];
      float phi0 = meas_package.raw_measurements_[1];
      float v0 =  meas_package.raw_measurements_[2];

      float x0 = rho0*cos(phi0);
      float y0 = rho0*sin(phi0);
      float vx0 = v0*cos(phi0);
      float vy0 = v0*sin(phi0);

      x_ << x0, y0, vx0, vy0, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1],
            1,
            1,
            0;
    }

    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  //cout << "PROCESS MEASUREMENT END" << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  GenerateSigmaPoints(&Xsig_pred_);
  AugmentedSigmaPoints(&Xsig_pred_);
  SigmaPointPrediction(delta_t, &Xsig_pred_);
  PredictMeanAndCovariance(&x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  PredictRadarMeasurement(&z, &S);
  UpdateStateRadar(Zsig, z_pred, z, S);
}


void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  lambda_ = 3 - n_x_;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  Xsig.col(0) = x_;
  A *= sqrt(lambda_ + n_x_);
  for (int i = 0; i < n_x_; ++i) {
    Xsig.col(1 + i      )  = x_ + A.col(i);
    Xsig.col(1 + i + n_x_) = x_ - A.col(i);
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig = " << std::endl << Xsig << std::endl;

  //write result
  *Xsig_out = Xsig;

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {


  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //create augmented mean state
  x_aug << x_, 0, 0;

  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_+0,n_x_+0) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //std::cout << P_aug << std::endl;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  A *= sqrt(lambda_ + n_aug_);

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i < n_aug_; ++i) {
    Xsig_aug.col(1+i      ) = x_aug + A.col(i);
    Xsig_aug.col(1+i+n_aug_) = x_aug - A.col(i);
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;

}



void UKF::SigmaPointPrediction(double delta_t, MatrixXd* Xsig_aug_) {


/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  MatrixXd Xsig_aug = *Xsig_aug_;

  for (int i = 0; i < 2*n_aug_+1; ++i) {
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
                  v        + 0                                            + delta_t*noise_v,
                  yaw    + dyawdt*delta_t                               + 0.5*delta_t*delta_t*noise_yaw,
                  dyawdt + 0                                            + delta_t*noise_yaw;

    }
    Xsig_pred_.col(i) = Xpred_i;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  VectorXd x;
  MatrixXd P;
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //set weights
  weights_(0) = lambda_/(lambda_ + n_aug_);
  int n_sigma = 2*n_aug_ + 1;
  for (int i=1; i<n_sigma; ++i) {
    weights_(i) = 1.0/(2*(lambda_ + n_aug_));
  }

  //predict state mean
  x *= 0.0;
  for (int i=0; i<n_sigma; ++i) {
    x += weights_(i)*Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P *= 0.0;
  VectorXd r;
  for (int i=0; i<n_sigma; ++i) {
    r = Xsig_pred_.col(i) - x;
    P += weights_(i)*r*r.transpose();
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}



void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug_;

  int n_sigma = 2*n_aug_ + 1;
//  //set vector for weights
//  VectorXd weights = VectorXd(2*n_aug+1);
//  double weight_0 = lambda/(lambda+n_aug);
//  weights(0) = weight_0;
//  for (int i=1; i<2*n_aug+1; i++) {
//    double weight = 0.5/(n_aug+lambda);
//    weights(i) = weight;
//  }
//
//  //radar measurement noise standard deviation radius in m
//  double std_radr = 0.3;
//
//  //radar measurement noise standard deviation angle in rad
//  double std_radphi = 0.0175;
//
//  //radar measurement noise standard deviation radius change in m/s
//  double std_radrd = 0.1;
//
//  //create example matrix with predicted sigma points
//  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
//  Xsig_pred <<
//            5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
//      1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
//      2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
//      0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
//      0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sigma);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space

  double px, py, v, yaw, yaw_rate, wi;
  double rho, phi, rhodot;

  z_pred *= 0;
  for (int i = 0; i < n_sigma; ++i) {
    px = Xsig_pred_(0, i);
    py = Xsig_pred_(1, i);
    v = Xsig_pred_(2, i);
    yaw = Xsig_pred_(3, i);
    yaw_rate = Xsig_pred_(4, i);
    wi = weights_(i);

    rho = std::sqrt(px*px + py*py);
    phi = std::atan2(py, px);
    rhodot = 0;
    if (std::fabs(rho) > 1.0e-4) {
      rhodot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
    }

    Zsig(0,i) = rho;
    Zsig(1,i) = phi;
    Zsig(2,i) = rhodot;

    //calculate mean predicted measurement
    z_pred += wi*Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd R(3,3);
  R <<
    std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;

  VectorXd diff(3);
  S *= 0;
  for (int i = 0; i < n_sigma; ++i) {
    diff = Zsig.col(i) - z_pred;
    S += weights_(i)*diff*diff.transpose();
  }

  S += R;

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}


void UKF::UpdateStateRadar(MatrixXd Zsig, VectorXd z_pred, VectorXd z, MatrixXd S) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

//  //define spreading parameter
//  double lambda = 3 - n_aug;
//
//  //set vector for weights
//  VectorXd weights = VectorXd(2*n_aug+1);
//  double weight_0 = lambda/(lambda+n_aug);
//  weights(0) = weight_0;
//  for (int i=1; i<2*n_aug+1; i++) {
//    double weight = 0.5/(n_aug+lambda);
//    weights(i) = weight;
//  }
//
//  //create example matrix with predicted sigma points in state space
//  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
//  Xsig_pred <<
//            5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
//      1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
//      2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
//      0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
//      0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
//
//  //create example vector for predicted state mean
//  VectorXd x = VectorXd(n_x);
//  x <<
//    5.93637,
//      1.49035,
//      2.20528,
//      0.536853,
//      0.353577;
//
//  //create example matrix for predicted state covariance
//  MatrixXd P = MatrixXd(n_x,n_x);
//  P <<
//    0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
//      -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
//      0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
//      -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
//      -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;
//
//  //create example matrix with sigma points in measurement space
//  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
//  Zsig <<
//       6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
//      0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
//      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;
//
//  //create example vector for mean predicted measurement
//  VectorXd z_pred = VectorXd(n_z);
//  z_pred <<
//         6.12155,
//      0.245993,
//      2.10313;
//
//  //create example matrix for predicted measurement covariance
//  MatrixXd S = MatrixXd(n_z,n_z);
//  S <<
//    0.0946171, -0.000139448,   0.00407016,
//      -0.000139448,  0.000617548, -0.000770652,
//      0.00407016, -0.000770652,    0.0180917;
//
//  //create example vector for incoming radar measurement
//  VectorXd z = VectorXd(n_z);
//  z <<
//    5.9214,   //rho in m
//      0.2187,   //phi in rad
//      2.0062;   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  int n_sigma = 2*n_aug_ + 1;
  VectorXd diff_state, diff_meas;
  Tc *= 0.0;
  for (int i = 0; i < n_sigma; ++i) {
    diff_state = Xsig_pred_.col(i) - x_;
    diff_meas = Zsig.col(i) - z_pred;
    Tc += weights_(i)*diff_state*diff_meas.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc*S.inverse();

  //update state mean and covariance matrix
  x_ += K*(z - z_pred);
  P_ -= K*S*K.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}
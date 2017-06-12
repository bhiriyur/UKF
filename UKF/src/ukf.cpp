#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


UKF::UKF() {

    /* =========================================================================
     * Initialize class variables
     * =========================================================================
     */

    use_laser_ = true;    // as name indicates
    use_radar_ = true;    // as name indicates
    n_x_ = 5;             // Number of original sigma points
    n_aug_ = 7;           // Number of augmented sigma points (accounting for process noise)
    std_a_ = 1.8;         // Process noise standard deviation longitudinal acceleration in m/s^2
    std_yawdd_ = 0.6;     // Process noise standard deviation yaw acceleration in rad/s^2
    std_laspx_ = 0.15;    // Laser measurement noise standard deviation position1 in m
    std_laspy_ = 0.15;    // Laser measurement noise standard deviation position2 in m
    std_radr_ = 0.3;      // Radar measurement noise standard deviation radius in m
    std_radphi_ = 0.03;   // Radar measurement noise standard deviation angle in rad
    std_radrd_ = 0.3;     // Radar measurement noise standard deviation radius change in m/s
    DEBUG = false;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    if (!is_initialized_) {
        /*
         * FIRST TIME ONLY
         */

        x_ = VectorXd(n_x_);
        P_ = MatrixXd(n_x_, n_x_);

        P_ <<   1,   0,    0,    0,    0,
                0,   1,    0,    0,    0,
                0,   0,    1,    0,    0,
                0,   0,    0,    1,    0,
                0,   0,    0,    0,    1;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            double rho0 = meas_package.raw_measurements_[0];
            double phi0 = meas_package.raw_measurements_[1];
            double v0 =  meas_package.raw_measurements_[2];

            double x0 = rho0*cos(phi0);
            double y0 = rho0*sin(phi0);
            double vx0 = v0*cos(phi0);
            double vy0 = v0*sin(phi0);

            x_ << x0, y0, vx0, vy0, 0;

        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

            //set the state with the initial location and zero velocity
            x_ << meas_package.raw_measurements_[0],
                    meas_package.raw_measurements_[1],
                    0,
                    0,
                    0;
        }

        // done initializing, no need to predict or update
        previous_timestamp_ = meas_package.timestamp_;
        is_initialized_ = true;
        return;
    }

    // Compute delta_t and store timestamp for future use
    delta_t = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    Prediction();

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        if (use_radar_) {
            UpdateRadar(meas_package);
        }
    } else {
        if (use_laser_) {
            UpdateLidar(meas_package);
        }
    }
}

// ================================================================================================
void UKF::Prediction() {

    AugmentedSigmaPoints();
    SigmaPointPrediction();
    PredictMeanAndCovariance();
}

// ================================================================================================
void UKF::UpdateLidar(MeasurementPackage meas_package) {

    PredictLidarMeasurement();
    VectorXd z(2);
    z <<    meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1);
    UpdateStateLidar(z);
}

// ================================================================================================
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    PredictRadarMeasurement();
    VectorXd z(3);
    z <<    meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1),
            meas_package.raw_measurements_(2);
    UpdateStateRadar(z);
}


// ================================================================================================
void UKF::AugmentedSigmaPoints() {

    //define spreading parameter
    lambda_ = 3 - n_aug_;

    n_sigma_ = 2*n_aug_ + 1;

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma_);

    //create augmented mean state
    x_aug << x_, 0, 0;

    //create augmented covariance matrix
    P_aug *= 0;
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_+0,n_x_+0) = std_a_*std_a_;
    P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

    //create square root matrix
    MatrixXd A = P_aug.llt().matrixL();
    A *= sqrt(lambda_ + n_aug_);

    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i=0; i < n_aug_; ++i) {
        Xsig_aug.col(1+i      ) = x_aug + A.col(i);
        Xsig_aug.col(1+i+n_aug_) = x_aug - A.col(i);
    }

    //write result
    Xsig_pred_= Xsig_aug;

    //print result
    if (DEBUG) std::cout << "***** AugmentedSigmaPoints *****" << std::endl;
}

// ================================================================================================
void UKF::SigmaPointPrediction() {
    
    MatrixXd Xsig_aug = Xsig_pred_;

    MatrixXd Xsig_out = MatrixXd(n_x_, n_sigma_);

    for (int i = 0; i < n_sigma_; ++i) {
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
        Xsig_out.col(i) = Xpred_i;
    }

    // write result
    Xsig_pred_ = Xsig_out;

    //print result
    if (DEBUG) std::cout << "***** SigmaPointPrediction *****" << std::endl;
}

// ================================================================================================
void UKF::PredictMeanAndCovariance() {


    VectorXd x = VectorXd(n_x_);
    MatrixXd P = MatrixXd(n_x_, n_x_);

    //set weights
    weights_ = VectorXd(n_sigma_);
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

    //write result
    x_ = x;
    P_ = P;

    //print result
    if (DEBUG) {
        std::cout << "***** PredictMeanAndCovariance *****" << std::endl;
        std::cout << "Predicted state" << std::endl;
        std::cout << x_ << std::endl;
        std::cout << "Predicted covariance matrix" << std::endl;
        std::cout << P_ << std::endl;
    }

}

// ================================================================================================
void UKF::PredictLidarMeasurement(){

    int n_z = 2;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    //transform sigma points into measurement space

    double px, py, v, yaw, yaw_rate, wi;

    z_pred *= 0;
    for (int i = 0; i < n_sigma_; ++i) {
        px = Xsig_pred_(0, i);
        py = Xsig_pred_(1, i);
        v = Xsig_pred_(2, i);
        yaw = Xsig_pred_(3, i);
        yaw_rate = Xsig_pred_(4, i);
        wi = weights_(i);

        Zsig(0,i) = px;
        Zsig(1,i) = py;

        //calculate mean predicted measurement
        z_pred += wi*Zsig.col(i);
    }

    //calculate measurement covariance matrix S
    MatrixXd R(n_z, n_z);
    R <<
      std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;

    VectorXd diff(n_z);
    S *= 0;
    for (int i = 0; i < n_sigma_; ++i) {
        diff = Zsig.col(i) - z_pred;
        S += weights_(i)*diff*diff.transpose();
    }

    S += R;

    //write result
    Zsig_ = Zsig;
    z_pred_ = z_pred;
    S_ = S;

    //print result
    if (DEBUG) {
        std::cout << "***** PredictLidarMeasurement *****" << std::endl;
        std::cout << "z_pred: " << std::endl << z_pred_ << std::endl;
        std::cout << "S: " << std::endl << S_ << std::endl;
    }
}

// ================================================================================================
void UKF::PredictRadarMeasurement() {

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, n_sigma_);

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    //transform sigma points into measurement space

    double px, py, v, yaw, yaw_rate, wi;
    double rho, phi, rhodot;

    z_pred *= 0;
    for (int i = 0; i < n_sigma_; ++i) {
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

        phi = normalize_angle(phi);

        Zsig(0,i) = rho;
        Zsig(1,i) = phi;
        Zsig(2,i) = rhodot;

        //calculate mean predicted measurement
        z_pred += wi*Zsig.col(i);
    }

    //calculate measurement covariance matrix S
    MatrixXd R(n_z,n_z);
    R <<
            std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;

    VectorXd diff(n_z);
    S *= 0;
    for (int i = 0; i < n_sigma_; ++i) {
        diff = Zsig.col(i) - z_pred;
        diff(1) = normalize_angle(diff(1));

        S += weights_(i)*diff*diff.transpose();
    }

    S += R;

    //write result
    Zsig_ = Zsig;
    z_pred_ = z_pred;
    S_ = S;

    //print result
    if (DEBUG) {
        std::cout << "***** PredictRadarMeasurement *****" << std::endl;
        std::cout << "z_pred: " << std::endl << z_pred_ << std::endl;
        std::cout << "S: " << std::endl << S_ << std::endl;
    }
}

// ================================================================================================
void UKF::UpdateStateLidar(VectorXd z) {

    //set measurement dimension, lidar can measure px, py
    int n_z = 2;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    int n_sigma = 2*n_aug_ + 1;
    VectorXd diff_state, diff_meas;
    Tc *= 0.0;
    for (int i = 0; i < n_sigma; ++i) {
        diff_state = Xsig_pred_.col(i) - x_;
        diff_meas = Zsig_.col(i) - z_pred_;
        Tc += weights_(i)*diff_state*diff_meas.transpose();
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc*S_.inverse();

    //update state mean and covariance matrix
    VectorXd zdiff = z - z_pred_;
    x_ += K*zdiff;
    P_ -= K*S_*K.transpose();

    //Update Lidar NIS
    NIS_laser_ = zdiff.transpose()*S_.inverse()*zdiff;

    //print result
    if (DEBUG) {
        std::cout << "***** UpdateStateLaser *****" << std::endl;
        std::cout << "Updated state x: " << std::endl << x_ << std::endl;
        std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
        std::cout << "NIS Lidar : " << NIS_laser_ << std::endl;
    }
}

// ================================================================================================
void UKF::UpdateStateRadar(VectorXd z) {

    //set measurement dimension, radar can measure r, phi, and r_dot
    int n_z = 3;

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    int n_sigma = 2*n_aug_ + 1;
    VectorXd diff_state, diff_meas;
    Tc *= 0.0;
    for (int i = 0; i < n_sigma; ++i) {
        diff_state = Xsig_pred_.col(i) - x_;
        diff_meas = Zsig_.col(i) - z_pred_;
        diff_meas(1) = normalize_angle(diff_meas(1));

        Tc += weights_(i)*diff_state*diff_meas.transpose();
    }

    //calculate Kalman gain K;
    MatrixXd K = Tc*S_.inverse();

    //update state mean and covariance matrix
    VectorXd zdiff = z - z_pred_;
    zdiff(1) = normalize_angle(zdiff(1));

    x_ += K*zdiff;
    P_ -= K*S_*K.transpose();

    //Update Radar NIS
    NIS_radar_ = zdiff.transpose()*S_.inverse()*zdiff;

    //print result
    if (DEBUG) {
        std::cout << "***** UpdateStateRadar *****" << std::endl;
        std::cout << "Updated state x: " << std::endl << x_ << std::endl;
        std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
        std::cout << "NIS Radar : " << NIS_radar_ << std::endl;
    }
}


// ================================================================================================
double UKF::normalize_angle(double val) {
    // Normalize angle to be in range -M_PI to +M_PI

    bool normalized = false;
    if (DEBUG and std::fabs(val)>M_PI) {
        std::cout << "NORMALIZING VAL " << val;
        normalized = true;
    }

    // normalize
    double out = val;
    while (out >  M_PI) out -= 2*M_PI;
    while (out < -M_PI) out += 2*M_PI;

    if (normalized) std::cout << " TO " << out << std::endl;

    return out;
}
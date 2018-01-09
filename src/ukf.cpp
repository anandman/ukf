#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::IOFormat;
using std::vector;

#define PI 3.14159265358979323
#define PI2 (2.0 * PI)

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //time when the state is true, in us
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_+n_aug_);

  // initial state vector
  // px, py, vel, yaw, yawd
  x_ = VectorXd::Zero(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() = default;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      x_ << cos(phi) * rho, sin(phi) * rho, sqrt(cos(phi) * rho_dot + sin(phi) * rho_dot), 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize the state with the initial location and zero velocity.
      */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // ignore LIDAR or RADAR measurements if specified
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) return;
  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) return;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  //compute the time elapsed between the current and previous measurements
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;  //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  Update(meas_package);

  // print the output
  IOFormat xFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "x_ = [", "]");
  IOFormat pFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "[", "]", "P_ = [", "]");
  cout << x_.format(xFmt) << endl;
  cout << P_.format(pFmt) << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // augmented state mean vector
  VectorXd x_aug_ = VectorXd::Zero(n_aug_);
  x_aug_.head(n_x_)=x_;

  // augmented state covariance matrix
  MatrixXd P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;
  MatrixXd sqrtP_aug_ = P_aug_.llt().matrixL();

  // augmented sigma points matrix
  MatrixXd Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.col(0) = x_aug_;
  double ofs = sqrt(lambda_ + n_aug_);
  for (int i=0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1) = x_aug_ + ofs * sqrtP_aug_.col(i);
    Xsig_aug_.col(n_aug_+i+1) = x_aug_ - ofs * sqrtP_aug_.col(i);
  }

  // predict sigma points
  Xsig_pred_.fill(0.0);
  for (int i=0; i < Xsig_aug_.cols(); i++) {
    //double px = Xsig_aug_(0,i);
    //double py = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yaw_dot = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yaw = Xsig_aug_(6,i);

    // calculate prediction integral
    VectorXd integral = VectorXd::Zero(n_x_);
    // check for yaw_dot = 0 which implies car is moving in a straight line
    if (fabs(yaw_dot) < 0.001) {
      integral << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              0, //yaw_dot * delta_t,
              0;
    }
    else {
      integral << v/yaw_dot * (sin(yaw + (yaw_dot*delta_t)) - sin(yaw)),
              v/yaw_dot * (cos(yaw) - cos(yaw + (yaw_dot*delta_t))),
              0,
              yaw_dot * delta_t,
              0;
    }

    // calculate noise influence on predicted state
    VectorXd noise = VectorXd::Zero(n_x_);
    noise << 0.5 * pow(delta_t,2) * cos(yaw) * nu_a,
            0.5 * pow(delta_t,2) * sin(yaw) * nu_a,
            delta_t * nu_a,
            0.5 * pow(delta_t,2) * nu_yaw,
            delta_t * nu_yaw;

    Xsig_pred_.col(i) = Xsig_aug_.col(i).head(5) + integral + noise;
  }

  // mean of predicted sigma points
  x_.fill(0.0);
  P_.fill(0.0);

  x_ = Xsig_pred_ * weights_;

  for (int i=0; i < Xsig_pred_.cols(); i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3) < -PI) x_diff(3) += PI2;
    while (x_diff(3) > PI) x_diff(3) -= PI2;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
* Updates the state and the state covariance matrix using a radar or laser measurement
* @param n_z measurement dimensions (RADAR[0:2] = rho, phi, rho_dot; LIDAR[0:1] = px, py)
* @param meas_package The measurement at k+1
*/
void UKF::Update(MeasurementPackage meas_package) {
  int n_z;
  bool is_radar;

  // figure out whether we're operating on laser or radar
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // RADAR updates
    n_z = 3;
    is_radar = true;
  } else {
    // LIDAR updates
    n_z = 2;
    is_radar = false;
  }

  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
    // extract values
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    if (is_radar) {
      double v  = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);
      double v1 = v * cos(yaw);
      double v2 = v * sin(yaw);

      // radar measurement model
      Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //rho
      Zsig(1,i) = atan2(p_y,p_x);                                 //phi
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //rho_dot
    } else {
      // laser measurement model
      Zsig(0,i) = p_x;
      Zsig(1,i) = p_y;
    }
  }

  // calculate predicted state mean in measurement space
  VectorXd z_ = Zsig * weights_;

  // normalize angles
  if (is_radar) {
    while (z_(1) > PI) z_(1) -= 2. * PI;
    while (z_(1) < -PI) z_(1) += 2. * PI;
  }

  // calculate predicted state covariance in measurement space
  MatrixXd S_ = MatrixXd::Zero(n_z, n_z);

  for (int i=0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_;

    //normalize angles
    if (is_radar) {
      while (z_diff(1) < -PI) z_diff(1) += PI2;
      while (z_diff(1) > PI) z_diff(1) -= PI2;
    }

    S_ += weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd::Zero(n_z, n_z);
  if (is_radar) {
    R_(0,0) = std_radr_ * std_radr_;
    R_(1,1) = std_radphi_ * std_radphi_;
    R_(2,2) = std_radrd_ * std_radrd_;
  } else {
    R_(0,0) = std_laspx_ * std_laspx_;
    R_(1,1) = std_laspy_ * std_laspy_;
  }
  S_ += R_;

  // calculate cross-correlation between state and measurement spaces
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);
  for (int i=0; i < Zsig.cols(); i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalize angles
    while (x_diff(3) < -PI) x_diff(3) += PI2;
    while (x_diff(3) > PI) x_diff(3) -= PI2;

    VectorXd z_diff = Zsig.col(i) - z_;
    if (is_radar) {
      //normalize angles
      while (z_diff(1) < -PI) z_diff(1) += 2 * PI;
      while (z_diff(1) > PI) z_diff(1) -= 2 * PI;
    }

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // calculate Kalman gain
  MatrixXd K = Tc * S_.inverse();

  // calculate error
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_diff = z - z_;
  if(is_radar) {
    while (z_diff(1) < -PI) z_diff(1) += PI2;
    while (z_diff(1) > PI) z_diff(1) -= PI2;
  }

  // update state
  x_ += K * z_diff;
  P_ -= K * S_ * K.transpose();

  // normalize angles
  while (x_(3)> PI) x_(3)-=PI2;
  while (x_(3)< -PI) x_(3)+=PI2;

  // calculate NIS
  double nis = z_diff.transpose() * S_.inverse() * z_diff;
  cout << "NIS = " << nis << endl;
}
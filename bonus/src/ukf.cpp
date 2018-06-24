#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2; // 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2; // 30;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  time_us_ = 0;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 + n_x_;

  radar_nis_ = 0;

  lidar_nis_ = 0;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);

  is_initialized_ = false;
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
  if (!is_initialized_) {
    time_us_ = meas_package.timestamp_;

    x_ << 1, 1, 1, 1, 0.1;

    P_ << 0.15,   0,  0, 0, 0,
             0, 0.15, 0, 0, 0,
             0,   0,  1, 0, 0,
             0,   0,  0, 1, 0,
             0,   0,  0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "Initialize with Radar" << endl;
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << "Initialize with Laser" << endl;
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }

    is_initialized_ = true;
    return;
  }

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
	time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }

  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
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

  // Generate Sigma Points
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  MatrixXd P_sqrt = P_.llt().matrixL();

  Xsig.col(0) = x_;
  for (int i = 0; i < n_x_; i++) {
    Xsig.col(i+1)      = x_ + sqrt(lambda_+n_x_) * P_sqrt.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * P_sqrt.col(i);
  }

  /// Augment Sigma Points
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  double lambda_aug = 3 - n_aug_;

  /// Create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  /// Create augmented covar matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  /// Create square root matrix
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  /// Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++) {
      Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_aug+n_aug_) * P_aug_sqrt.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug+n_aug_) * P_aug_sqrt.col(i);
  }

  // Sigma Point Prediction
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    /// Predicted values
    double px_p, py_p;
    double v_p, yaw_p, yawd_p;

    if (fabs(yawd) > 0.001) {
      px_p = px + (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + (v / yawd) * (cos(yaw) - cos(yaw + (yawd * delta_t)));
    }
    else {
      px_p = px + (v * delta_t * cos(yaw));
      py_p = py + (v * delta_t * sin(yaw));
    }

    v_p = v;
    yaw_p = yaw + (yawd * delta_t);
    yawd_p = yawd;

    /// Incorporate noise
    px_p = px_p + (0.5 * nu_a * delta_t*delta_t) * cos(yaw);
    py_p = py_p + (0.5 * nu_a * delta_t*delta_t) * sin(yaw);
    v_p = v_p + (nu_a * delta_t);

    yaw_p = yaw_p + (0.5 * nu_yawdd * delta_t*delta_t);
    yawd_p = yawd_p + (nu_yawdd * delta_t);

    /// Store sigma points
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predict
  weights_(0) = lambda_aug / (lambda_aug + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_aug);
  }

  /// Calculate predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  /// Calculate predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    /// state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    /// angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  cout << "Predicted state" << endl << x_ << endl;
  cout << "Predicted covariance matrix" << endl << P_ << endl << endl;
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
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);

  // Convert sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    /// extract values
    double px  = Xsig_pred_(0,i);
    double py  = Xsig_pred_(1,i);

    /// measurement model: px, py
    Zsig(0, i) = px;
    Zsig(1, i) = py;
  }

  // Calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //2n+1 simga points
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  lidar_nis_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
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

  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  // Convert sigma points into measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    /// extract values
    double px  = Xsig_pred_(0,i);
    double py  = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    /// measurement model: r, rho, rho_dot
    Zsig(0,i) = sqrt(px*px + py*py);
    Zsig(1,i) = atan2(py, px);
    Zsig(2,i) = ((px*v1) + (py*v2)) / sqrt(px*px + py*py);
  }

  // Calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    /// angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {  //2n+1 simga points
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1) > M_PI)  z_diff(1) -= 2.*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3) > M_PI)  x_diff(3) -= 2.*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1) > M_PI)   z_diff(1) -= 2.*M_PI;
  while (z_diff(1) < -M_PI)  z_diff(1) += 2.*M_PI;

  radar_nis_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

VectorXd UKF::FuturePrediction(double delta_t) {
  // Generate Sigma Points
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  MatrixXd P_sqrt = P_.llt().matrixL();

  Xsig.col(0) = x_;
  for (int i = 0; i < n_x_; i++) {
    Xsig.col(i+1)      = x_ + sqrt(lambda_+n_x_) * P_sqrt.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * P_sqrt.col(i);
  }

  /// Augment Sigma Points
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  double lambda_aug = 3 - n_aug_;

  /// Create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  /// Create augmented covar matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  /// Create square root matrix
  MatrixXd P_aug_sqrt = P_aug.llt().matrixL();

  /// Create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++) {
      Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_aug+n_aug_) * P_aug_sqrt.col(i);
      Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug+n_aug_) * P_aug_sqrt.col(i);
  }

  // Sigma Point Prediction
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double px       = Xsig_aug(0, i);
    double py       = Xsig_aug(1, i);
    double v        = Xsig_aug(2, i);
    double yaw      = Xsig_aug(3, i);
    double yawd     = Xsig_aug(4, i);
    double nu_a     = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    /// Predicted values
    double px_p, py_p;
    double v_p, yaw_p, yawd_p;

    if (fabs(yawd) > 0.001) {
      px_p = px + (v / yawd) * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = py + (v / yawd) * (cos(yaw) - cos(yaw + (yawd * delta_t)));
    }
    else {
      px_p = px + (v * delta_t * cos(yaw));
      py_p = py + (v * delta_t * sin(yaw));
    }

    v_p = v;
    yaw_p = yaw + (yawd * delta_t);
    yawd_p = yawd;

    /// Incorporate noise
    px_p = px_p + (0.5 * nu_a * delta_t*delta_t) * cos(yaw);
    py_p = py_p + (0.5 * nu_a * delta_t*delta_t) * sin(yaw);
    v_p = v_p + (nu_a * delta_t);

    yaw_p = yaw_p + (0.5 * nu_yawdd * delta_t*delta_t);
    yawd_p = yawd_p + (nu_yawdd * delta_t);

    /// Store sigma points
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predict
  weights_(0) = lambda_aug / (lambda_aug + n_aug_);
  for (int i = 1; i < 2*n_aug_+1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_aug);
  }

  /// Calculate predicted state mean
  VectorXd x_guess = VectorXd(n_x_);
  x_guess.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_guess = x_guess + weights_(i) * Xsig_pred_.col(i);
  }

  return x_guess;
}


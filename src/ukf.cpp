#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
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
  std_a_ = 1.1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  is_initialized_ = false;
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  
  weights_ = VectorXd(2*n_aug_+1);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i = 1; i < weights_.size(); i++) {
    weights_(i) = 1/(2*(lambda_+n_aug_));
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  // If this is our first time receiving a measurement, we must initialize the state vector and covariance matrix before we can make any predictions
  if(!is_initialized_) {
    is_initialized_ = true;
   	P_.setIdentity(n_x_, n_x_);
    
    VectorXd raw_vals = meas_package.raw_measurements_;
    
    double px, py;
    
    switch(meas_package.sensor_type_) {
      case MeasurementPackage::LASER:
        // Laser measurements directly measure Cartesian coordinates
        px = raw_vals[0];
        py = raw_vals[1];

        // Initialize yaw and yaw rate with values so we converge faster
        x_ << px, py, 0.0, 0.3, 0.1;

        // The standard deviation for laser measurements is smaller, so we need not use the identity matrix directly.
        P_(0,0) = 0.55;
        P_(1,1) = 0.55;
        break;
      case MeasurementPackage::RADAR:
        double rho, phi, rhod;
        rho = raw_vals[0];
        phi = raw_vals[1];
        rhod = raw_vals[2];

        // Radar measurements are in polar coordinates and must be converted to Cartesian before use
        px = rho*cos(phi);
        py = rho*sin(phi);

        x_ << px, py, 0.0, 0.3, 0.1;

        // The standard deviation for the radar measurements is larger than that of laser measurements, so we may have more variance, but still less than 1 for px and py
        P_(0,0) = 0.75;
        P_(1,1) = 0.75;
        break;
      default:
        std::cout << "Unsupported sensor type in UKF::ProcessMeasurement at ukf.cpp" << std::endl;
        return;
    }
    
    // Initialize the first timestamp for use with prediction function
    time_us_ = meas_package.timestamp_;
    
    return;
  }
  
  // We must divide by 1 million here because the time provided by the meas_package is in microseconds
  double time_change = (meas_package.timestamp_ - time_us_)/1000000.0;
  
  // UKF follows a prediction -> update cycle
  Prediction(time_change);
  
  // Instead of having one update function with many conditionals, we can simply divide the update step based on the sensor type
  switch(meas_package.sensor_type_) {
    case MeasurementPackage::LASER:
      UpdateLidar(meas_package);
      break;
    case MeasurementPackage::RADAR:
      UpdateRadar(meas_package);
      break;
    default:
      std::cout << "Unsupported sensor type in UKF::ProcessMeasurement at ukf.cpp" << std::endl;
      return;
  }
  
  // Update the timestamp for subsequent calls to this function
  time_us_ = meas_package.timestamp_;
  
  return;
}

void UKF::Prediction(double delta_t) {
  // Credit to the instructors & their solutions
  // This function is based on the programming exercises that came in previous lessons, and has been tweaked using the solutions publicly given
  
  // Create an augmented state vector, covariance matrix, and sigma point matrix, which all include acceleration noise and yaw acceleration noise
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  
  // The first elements of the x_(aug) and P_(aug) vectors are identical
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  // Follow the formulas for determining sigma points
  MatrixXd P_root = P_aug.llt().matrixL();
  
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * P_root.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * P_root.col(i);
  }
  
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_p;
    double py_p;

    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  // Update the state vector and covariance matrix using the predicted sigma points, and weight each
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) { 
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) { 

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Normalize angle to be in between -PI and +PI
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // This function was not written with assistance from the instructors, but its implementation was extrapolated from that of UKF::UpdateRadar, which was explored in previous lessons

  // Only update based on lidar if we allow it
  if(!use_laser_) return;
  
  // Initialize z vector and related matrices to transform the measurement space to the state space
  int n_z = 2;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for(int i = 0; i < 2 * n_aug_ + 1; i++)
  {
   	double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }
  
  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }
  
  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  // Add noise from sensor standard deviation
  S = S + R;
  
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd z = meas_package.raw_measurements_;
  
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;

  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // Update x_ and P_ with a combination of the predicted sigma points and actual measurements 
  // This updates P_ hopefully leading to more accurate predictions for the state vector in the future
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Credit to the instructors & their solutions
  // This function is based on the programming exercises that came in previous lessons, and has been tweaked using the solutions publicly given
  
  // Only update based on radar if we allow it
  if(!use_radar_) return;
  
  // Similarly to lidar, create vectors and matrices to help us update
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) { 
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y); 
    Zsig(1,i) = atan2(p_y,p_x);
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  }
  
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;
  // Add noise due to sensor standard deviation
  S = S + R;
  
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  VectorXd z = meas_package.raw_measurements_;
  
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) { 
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Normalize angle to between -PI and +PI
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  MatrixXd K = Tc * S.inverse();


  VectorXd z_diff = z - z_pred;

  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}
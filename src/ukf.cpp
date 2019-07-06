#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter lol
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
  std_a_ = 1.5;

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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  
  weights_ = VectorXd(2*n_aug_+1);
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i = 1; i < weights_.size(); i++) {
    weights_(i) = 1/2*(lambda_+n_aug_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if(!is_initialized_) {
    is_initialized_ = true;
   	P_.setIdentity(n_x_, n_x_);
    
    VectorXd raw_vals = meas_package.raw_measurements_;
    
    double px, py;
    
    switch(meas_package.sensor_type_) {
      case MeasurementPackage::LASER:
        px = raw_vals[0];
        py = raw_vals[1];

        x_ << px, py, 0, 0, 0;
        break;
      case MeasurementPackage::RADAR:
        double rho, phi, rhod;
        rho = raw_vals[0];
        phi = raw_vals[1];

        px = rho*cos(phi);
        py = rho*sin(phi);

        x_ << px, py, 0, 0, 0;
        break;
      default:
        return;
    }
    
    time_us_ = meas_package.timestamp_;
    
    return;
  }
  
  double delta_t = meas_package.timestamp_ - time_us_;
  Prediction(delta_t);
  
  switch(meas_package.sensor_type_) {
    case MeasurementPackage::LASER:
      UpdateLidar(meas_package);
      break;
    case MeasurementPackage::RADAR:
      UpdateRadar(meas_package);
      break;
    default:
      return;
  }
  
  time_us_ = meas_package.timestamp_;
  
  return;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}
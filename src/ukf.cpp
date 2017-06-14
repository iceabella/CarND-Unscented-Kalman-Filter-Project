#include "ukf.h"
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
  /// initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // time when the state is true, in us
  time_us_ = 0;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3-n_x_; // design parameter

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  //set weights
  /*weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i=1; i<2*n_aug_+1; i++)
    weights_(i) = 0.5/(lambda_ + n_aug_);*/

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;	// tunable parameter

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;	// tunable parameter

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

  1. Predict
  2. Update with new measurement
  */
  // if first measurement initialize only
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;
    //x_ = VectorXd(n_x_);

    // save time stamp of measurement
    time_us_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	double rho = meas_package.raw_measurements_(0);
	double phi = meas_package.raw_measurements_(1); 

	x_(0) = rho*cos(phi);
	x_(1) = rho*sin(phi);
	x_(2) = 0.0;
	x_(3) = 0.0; 
	x_(4) = 0.0; // we cannot calculate v, yaw or yaw rate based on only one radar measurement
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	x_(0) = meas_package.raw_measurements_(0);
	x_(1) = meas_package.raw_measurements_(1);
	x_(2) = 0.0;
	x_(3) = 0.0; 
	x_(4) = 0.0; // we cannot measure v, yaw or yaw rate for laser
    }
    // Initialize covariance matrix
    //P_ = MatrixXd(n_x_,n_x_);
    
    P_ << 0.8, 0, 0, 0, 0,
	  0, 0.8, 0, 0, 0,
	  0, 0, 10, 0, 0,
	  0, 0, 0, 10, 0,
	  0, 0, 0, 0, 10; // we are very uncertain of our starting values

    // Initialize weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for(int i=1; i<2*n_aug_+1; i++)
    	weights_(i) = 0.5/(lambda_ + n_aug_);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  // Define time interval between measurements
  double dt = (meas_package.timestamp_-time_us_)/1000000.0; // s WHAT HAPPENS IF WE DO NOT USE THE MEASUREMENT?
  time_us_ = meas_package.timestamp_;

  // 1. PREDICT
  Prediction(dt);

  // 2. UPDATE BASED ON NEW MEASUREMENT

  // RADAR MEASUREMENT
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
	UpdateRadar(meas_package);

  } // LIDAR MEASUREMENT
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
	UpdateLidar(meas_package);
  }
  else
	std::cout << "ERROR: Invalid sensor type!" << std::endl;

  return;
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

  1. Generate augmented sigma points
  2. Predict sigma points
  3. Predict mean and covariance
  */

  // 1. GENERATE AUGMENTED SIGMA POINTS
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented covariance matrix
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  MatrixXd tempMat = sqrt(lambda_+n_aug_)*A;
  
  for(int i=0; i<=n_aug_-1; i++){
    Xsig_aug.col(i+1) = x_aug + tempMat.col(i);
    Xsig_aug.col(n_aug_+1+i) = x_aug - tempMat.col(i);
  }
  
  // 2. PREDICT SIGMA POINTS
  
  /*double px, py, v, psi, psi_dot;
  double mu_a, mu_psi;
  double px_k1,py_k1, v_k1,psi_k1, psi_dot_k1;*/

  Xsig_pred_.fill(0.0);
  
  //predict sigma points
  for(int i=0;i<2*n_aug_+1;i++){
	// define variables for readability        
	double px = Xsig_aug(0,i);
        double py = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double psi = Xsig_aug(3,i);
        double psi_dot = Xsig_aug(4,i);
        double mu_a = Xsig_aug(5,i);
        double mu_psi = Xsig_aug(6,i);
	double px_k1,py_k1, v_k1,psi_k1, psi_dot_k1; // predicted sigma point
        
        //avoid division by zero
        if(fabs(psi_dot) < 0.001){
            px_k1 = px + v*cos(psi)*delta_t + 0.5*delta_t*delta_t*cos(psi)*mu_a;
            py_k1 = py + v*sin(psi)*delta_t + 0.5*delta_t*delta_t*sin(psi)*mu_a;            
        } else{
            px_k1 = px + v/psi_dot*(sin(psi+psi_dot*delta_t)-sin(psi)) + 0.5*delta_t*delta_t*cos(psi)*mu_a;
            py_k1 = py + v/psi_dot*(-cos(psi+psi_dot*delta_t)+cos(psi)) + 0.5*delta_t*delta_t*sin(psi)*mu_a;
        }

	v_k1 = v + delta_t*mu_a;
        psi_k1 = psi + psi_dot*delta_t + 0.5*delta_t*delta_t*mu_psi;
        psi_dot_k1 = psi_dot + delta_t*mu_psi;
        
        //write predicted sigma points into right column
        Xsig_pred_.col(i) << px_k1, py_k1, v_k1, psi_k1, psi_dot_k1;
  }

  //3. PREDICT MEAN AND COVARIANCE

  //predict state mean
  x_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++)
    x_ += weights_(i)*Xsig_pred_.col(i);
  //predicted state covariance matrix
  P_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i)*x_diff*x_diff.transpose();
  }
  return;
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

  1. Predict measurement
  2. Update state

  Measurement model is linear why we will use a standard kalman filter, since this is cheaper computationally (according to lecture 7: 12)
  */
  /**
    * update the state by using Kalman Filter equations i.e. lidar data
  */
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(2);
  z = meas_package.raw_measurements_;

  // define measurement function
  MatrixXd H = MatrixXd(2,5);
  H << 1.0, 0, 0, 0, 0,
       0, 1.0, 0, 0, 0;	// we only measure position in x and y

  // define measurement noise
  MatrixXd R = MatrixXd(2,2);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  VectorXd y = z-H*x_;
  MatrixXd HT = H.transpose();
  MatrixXd S = H*P_*HT + R;
  MatrixXd K = P_*HT*S.inverse();
  MatrixXd I(n_x_,n_x_); 
  I = MatrixXd::Identity(n_x_,n_x_);

  // new state and covariance
  x_ = x_ + K*y;
  P_ = (I-K*H)*P_;

  // Check NIS
  //VectorXd epsilon = (z-z_pred).transpose()*S_inv*(z-z_pred);
  //std::cout << "NIS: " << epsilon << std::endl;

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
  
  1. Predict measurement
  2. Update state
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  // 1. PREDICT MEASUREMENT
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  //transform sigma points into measurement space
  for(int i=0; i<2*n_aug_+1; i++){
      // define some variables for readability
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double psi = Xsig_pred_(3,i);
      // transformation made with assumption w = 0 (z=h(x)+w)
      Zsig(0,i) = sqrt(px*px+py*py);
      Zsig(1,i) = atan2(py,px);
      Zsig(2,i) = (px*cos(psi)*v + py*sin(psi)*v)/Zsig(0,i);
  }
  //calculate mean predicted measurement
  for(int i=0; i<2*n_aug_+1; i++)
    z_pred += weights_(i)*Zsig.col(i);
  
  MatrixXd R = MatrixXd(3,3);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;
  //calculate measurement covariance matrix S
  for(int i=0; i<2*n_aug_+1; i++){
      // calculate residual
      VectorXd z_diff = Zsig.col(i)-z_pred;

      //angle normalization
      while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
      while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

      S += weights_(i)*z_diff*z_diff.transpose();
  }
  // add measurement noise covariance
  S += R;
  
  // 2. UPDATE STATE

  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0; i<2*n_aug_+1; i++){
    // calculate residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // calculate state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc += weights_(i)*x_diff*z_diff.transpose(); 
  }

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc*S_inv;

  //calculate residual
  VectorXd z_diff = z-z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K*(z-z_pred);
  P_ = P_ - K*S*K.transpose();

  // Check NIS
  VectorXd epsilon = (z-z_pred).transpose()*S_inv*(z-z_pred);
  //std::cout << "NIS: " << epsilon << std::endl;

  return;

}

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

  is_initialized_  =  false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  n_x_ = x_.size();

  n_aug_ = 7;

  lambda_ = 3.0 - n_aug_;

  Xsig_pred_ = MatrixXd (n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 
  //cout << is_initialized_ << endl; 

  if (!is_initialized_) {



    time_us_ = meas_package.timestamp_;

    // Setting initial state vector to zero except sidot.
    x_ << 1,1,0,0,0;

    //Setting initial covariance matrix to diagonal identity matrix 
    P_.fill(0.0);

    P_(0,0) = 2;
    P_(1,1) = 4;
    P_(2,2) = 0.5;
    P_(3,3) = 0.5;
    P_(4,4) = 0.25;


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      x_(0) = ro * cos(theta);
      x_(1) = ro * sin(theta);
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    //cout << "This is the initialized state " << x_ << endl;
    //cout << P_ <<endl;

    // Initialized
    is_initialized_ = true;

    //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__);

    return;
  }

  //printf ("This is line %d of file %s (function %s)\n",__LINE__, __FILE__, __func__);
  //cout << "State has been initialized to " << is_initialized_ << endl;

  //cout << meas_package.timestamp_ << endl;

  //cout << meas_package.sensor_type_ << endl ;
  //cout << meas_package.raw_measurements_ << endl;

  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;

  //cout << delta_t << endl;

  /*while (delta_t > 0.1 ) {
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  }
  */
  if (delta_t >= 0.01) {
    Prediction(delta_t);
  }  

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && use_radar_) {
    UpdateRadar(meas_package);
  }
  else if((meas_package.sensor_type_ == MeasurementPackage::LASER) && use_laser_) {
    UpdateLidar(meas_package);
  }
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
  //cout << "Predicting" << endl;

  //Generating sigma points

  VectorXd x_aug_  = VectorXd::Zero(n_aug_); // Augmented state vector for Posterior state 
  MatrixXd P_aug_ = MatrixXd::Zero(n_aug_,n_aug_); // Augmented covariance matrix 
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1); // Augmented state sigma point matrix
  Xsig_aug_.fill(0.0);

  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  ///cout << "Posterior state x is \n" << x_aug_ << endl;
  //cout << "POsterior state P_ is \n" << P_aug_ << endl; 

  MatrixXd L = P_aug_.llt().matrixL(); // Square root matrix

  //cout << "POsterior state square root matrix is \n" << L << endl;

  Xsig_aug_.col(0) = x_aug_;

  for (int i = 1; i <= n_aug_ ; i++) {
    Xsig_aug_.col(i) = x_aug_ + sqrt(lambda_+ n_aug_) * L.col(i-1);
    Xsig_aug_.col(i + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i-1); 
  }

  //cout << " Posterior state sigma points are \n" << Xsig_aug_ << endl;

  // Generating A priori state sigma points 
  //cout << Xsig_pred_ << endl;

  for (int i = 0; i < 2*n_aug_ + 1 ; i ++) {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else 
    {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //cout << " A priori sigma predicted points are \n" << Xsig_pred_ << endl;

  // Converting sigma points to mean vector and covariance

  //Defining the weights vector 
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i =1; i < 2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  //Predicted state mean
  VectorXd x_pred_ = VectorXd(5);
  x_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    x_pred_ = x_pred_ + weights_(i) * Xsig_pred_.col(i);
  }

  //cout << "A priori x mean is \n" << x_pred_ << endl;

  //Predicted state covariance 
  MatrixXd P_pred_ =  MatrixXd(5,5);
  P_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {

    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    P_pred_ = P_pred_ + weights_(i) * x_diff_ * x_diff_.transpose();
  }

  P_ = P_pred_;

  //cout << " A priori P covariance is \n" << P_pred_ << endl;

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

  //cout << "Updating Lidar \n" ;

  int n_lid_ = 2;
  MatrixXd Zsig_ = MatrixXd(n_lid_, 2*n_aug_+1);

  //Converitng Zsig_pred_ to measurement space
  for (int i = 0; i < 2*n_aug_+1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    Zsig_(0,i) = p_x;
    Zsig_(1,i) = p_y;
  }

  //cout << " Predicted sigma points in Lidar measurement space are \n" << Zsig_;

  //converting sigma points to mean state vector
  VectorXd z_pred_ = VectorXd(n_lid_);
  z_pred_.fill(0.0);

  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ + weights_(i)*Zsig_.col(i);
  }

  //cout << " Predicted state in Lidar measurement space is \n" << z_pred_ << endl;

  //Need to define new measurement covariance 
  MatrixXd S = MatrixXd(n_lid_, n_lid_);
  S.fill(0.0);

  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

    S = S + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  // Add measurement noise to S
  MatrixXd R = MatrixXd(n_lid_,n_lid_);
  R.fill(0.0);
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  S = S + R;

  //cout << " Predicted covariance in Lidar measurement space is \n" << S << endl;

  //Updating state and covariance vector

  // Creating matrix for correlation
  MatrixXd Tc = MatrixXd(n_x_, n_lid_);

  //Calculate Tc
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    //angle normalization
    //while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    //while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff_ * z_diff_.transpose();
  }

  //Kalman Gain
  MatrixXd K = Tc * S.inverse();

  //cout << "Kalman Gain (Lidar) is \n" << K << endl;

  //residual 
  VectorXd z_diff_ = meas_package.raw_measurements_ - z_pred_;

  //cout << " raw measurement (Lidar) \n" << meas_package.raw_measurements_  << endl;


  //cout << "z_diff (Lidar) \n" << z_diff_ << endl;

  //Updating 
  x_ = x_ + K*z_diff_;
  P_ = P_ - K*S*K.transpose();

  //cout << "Updated (Lidar) x_ is \n" << x_ << endl;

  //cout << "UPdated P_(Lidar) is \n" << P_ << endl;

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

  //cout << "Updating Radar \n";

  int n_rad_ = 3;
  MatrixXd Zsig_ = MatrixXd(n_rad_,2*n_aug_+1);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //2n+1 simga points
    // extract values for better readibility

    //cout << "Printing Predicted sigma points" << Xsig_pred_ << endl;

    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                       //r
    Zsig_(1,i) = atan2(p_y,p_x);
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot                 
    
  }

  //cout << " Predicted sigma points in Radar measurement space are \n" << Zsig_ << endl;

  //converting these sigma points into state vector
  VectorXd z_pred_ = VectorXd(n_rad_);
  z_pred_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //cout << "Predicted state in Radar measurement space is \n" << z_pred_ << endl;

  //COvariance in Measurement space
  MatrixXd S = MatrixXd(n_rad_, n_rad_);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  //add measurement noise
  MatrixXd R = MatrixXd(n_rad_, n_rad_);
  R.fill(0.0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  S = S + R;

  //cout << "Predicted covariance in Radar measurement space is \n" << S << endl;

  //Updating mean and covariance 
  MatrixXd Tc = MatrixXd(n_x_,n_rad_);
  Tc.fill(0.0);
  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    
    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;
     

    Tc = Tc + weights_(i) * x_diff_ * z_diff_.transpose();
    

  }

  //Kalman Gain
  MatrixXd K = Tc * S.inverse();
  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 


  //cout << "Kalman Gain (Radar) is \n" << K << endl;

  //residual
  VectorXd z_diff_ = meas_package.raw_measurements_ - z_pred_;
  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 

  //cout << " raw measurement (Radar) \n" << meas_package.raw_measurements_  << endl;


  //cout << "z_diff (Radar) \n" << z_diff_ << endl;

  //angle normalization
  while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 
  while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;
  //printf ("This is line %d of file %s (function %s)\n",\
                      __LINE__, __FILE__, __func__); 

  //updating
  x_ = x_ + K * z_diff_;
  P_ = P_ - K*S*K.transpose();

  //cout << "Updated (Radar) x_ is \n" << x_ << endl;

  //cout << "UPdated P_(Radar) is \n" << P_ << endl;
}

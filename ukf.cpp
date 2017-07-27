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
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  // Setting it to diagonal matrix with values 1
  P_ = MatrixXd(5, 5);
  P_.fill (0.0);
  for (int i = 0; i < 5; i++)
  {
    P_(i,i) = 1;
  }

  /***************************************/
    /* Need to tune later */ 

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30; 

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

   /* Need to tune later */ 
  /***************************************/

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

  is_initialized_ = false;
  use_radar_ = true;
  use_laser_ = true;

  ///* predicted sigma points matrix
  //MatrixXd Xsig_pred_; 
  ///* will be initialized and defined later 

  ///* Weights of sigma points
  //VectorXd weights_;
  ///* will be initialized and defined later 

  ///* State dimension
  n_x_ = x_.size();

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  double lambda_ = 3 - n_aug_;

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

  /***************************************************/
  // Initialization of First Measurement 
  /***************************************************/

  x_ << 1,1,1,1,0.5; // Random choice 


  if (!is_initialized_)
  {
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      float ro = meas_package.raw_measurements_[0];  // The distance from origin to the object
      float theta = meas_package.raw_measurements_[1]; // angle made by ro with the y axis;
      x_(0) = ro*cos(theta);
      x_(1) = ro*sin(theta); 
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {

      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    return;
  }
  /***************************************************/
  // First Measurement Initializated
  /***************************************************/


  /***************************************************/
  // If not the first measurement, we update the timestamp
  /***************************************************/
  time_us_ = meas_package.timestamp_ - time_us_;
  time_us_ = meas_package.timestamp_;
  /***************************************************/
  // Timestamp updated
  /***************************************************/



  /***************************************************/
  // Generating Sigma Points
  /***************************************************/
  ///*Creating Sigma Point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  ///* Creating Squareroot matrix
  MatrixXd A = P_.llt().matrixL();

  ///* Setting first column to be same as X
  Xsig.col(0) = x_;

  ///* Using the given formula to generate sigma points for the Posterior state.
  ///* Xk∣k =[xk∣k     xk∣k + sqrt((λ+nx)Pk∣k)  xk∣k - sqrt((λ+nx)Pk∣k​​)]
  for (int i = 1; i < n_x_ ; i++)
  {
    Xsig.col(i) = x_ + sqrt(lambda_ + n_x_) * A.col(i - 1);
    Xsig.col(i + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i - 1);
  }
  /***************************************************/
  // Generated Sigma Points
  /***************************************************/



  /***************************************************/
  // Generating Sigma Points for the Augmented State (including process noises)
  /***************************************************/
  ///* Creating a augmented state vector
  VectorXd x_aug = VectorXd(n_aug_); // 5 state variables and two process noise variables 

  ///* Create Augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  ///*Creating Augmented Sigma Point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  ///* Updating augmented mean state vector matrix with noise mean values 
  x_aug.head(5) = x_;
  x_aug(5) = 0; // mean of process noise is zero
  x_aug(6) = 0; // mean of process noise is zero

  ///* Filling in P_aug matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;  // Covariance matrix corresponding to state mean position 
  P_aug(5,5) = std_a_ * std_a_;   // Covariance matrix corresponding to state process noise 
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  ///* Creating square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  ///* Filling in Xsig_aug matrix (Generating augmented sigma points)
  Xsig_aug.col(0) = x_aug;

  for (int i = 1; i< n_aug_; i++)
  {
    Xsig_aug.col(i)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i-1);
  }
  /***************************************************/
  // Generated Sigma Points for the Augmented State in Xsig_aug matrix (includes process noises)
  /***************************************************/

  Prediction(time_us_/1000000.0);

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) 
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /***************************************************/
  // Predicting Sigma Points through the defined process model
  /***************************************************/

  ///* Process model equation
  ///* https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/0949fca6-b379-42af-a919-ee50aa304e6a/lessons/daf3dee8-7117-48e8-a27a-fc4769d2b954/concepts/63674ce2-43ed-418c-bf8b-d16ae73dffc0

  ///* Creating predicted sigma points matrix
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //extract values for better readability
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

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
  /***************************************************/
  // Predicted Sigma Points through the defined process model
  /***************************************************/



  /***************************************************/
  // Converting sigma points into Predicted Mean and covariance 
  /***************************************************/

  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);

  ///* Defining weights 
  ///* Weights wi=λ/(λ+na) for i=1 //na = 7
  ///* Weights wi=0.5/(λ+na) for i=2...nσ  //nσ = 15 
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;  // 1st weight 

  for (int i=1; i<2*n_aug_+1; i++) // 2nd to 2n+1 weights
  {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  ///* Converting predicted sigma points to predicted state vector and covariance matrix
  ///* Predicted Mean : xk+1∣k = ∑i=1 to nσ​  (wi * Xk+1∣k,i)

  //predicted state mean
  VectorXd x_pred_ = VectorXd(5);
  x_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points
    x_pred_ = x_pred_ + weights_(i) * Xsig_pred_.col(i);
  }

  ///* Predicted Covariance
  ///*Pk+1∣k = ∑i=1 to ​nσ (wi * (Xk+1∣k,i − xk+1∣k) * (Xk+1∣k,i − xk+1∣k)T
  MatrixXd P_pred_ = MatrixXd(5,5);
  P_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred_ = P_pred_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  /***************************************************/
  // Converted sigma points into Predicted Mean (x_pred_) and covariance (P_pred_) 
  /***************************************************/

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package.raw_measurements_);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package.raw_measurements_);
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(VectorXd z) 
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /***************************************************/
  // Converting the predicted sigma points into the measurement state space and 
  // then again changing them into mean vector and covariance model (not sigma points format)
  /***************************************************/
  
  ///* we use the same sigma points as predicted before Xsig_pred;transform sigma points into measurement space
  ///* Measurement Model : zk+1∣k = h(xk+1) + wk+1  (h is the conversion matrix and w is the measurement noise)


  int n_lid = 2; // set measurement dimension, lidar can measure px andn py only

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_lid, 2 * n_aug + 1); // 1 column for each sigma point 

  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
  
    // measurement model
    Zsig(0,i) = p_x;     // px
    Zsig(1,i) = p_y;     //py
  }

  //mean predicted measurement in measurement space 
  // Predicted Mean : zk+1∣k = ∑i=1 to nσ​  (wi * Zk+1∣k,i)
  VectorXd z_pred = VectorXd(n_lid);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) 
  {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  ///* Sk+1∣k = ∑i=1 to ​nσ (wi * (Zk+1∣k,i − zk+1∣k) * (Zk+1∣k,i − zk+1∣k)T  +  R;
  /***** we dont use augmentation here because the noise R is just additive part of the funtion ******/
  MatrixXd S = MatrixXd(n_lid,n_lid); 
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_lid,n_lid);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  /***************************************************/
  // Updating 
  /***************************************************/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_lid);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();


    //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(VectorXd z) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */

  
  /***************************************************/
  // Converting the predicted sigma points into the measurement state space and 
  // then again changing them into mean vector and covariance model (not sigma points format)
  /***************************************************/
  
  ///* we use the same sigma points as predicted before Xsig_pred;transform sigma points into measurement space
  ///* Measurement Model : zk+1∣k = h(xk+1) + wk+1  (h is the conversion matrix and w is the measurement noise)


  int n_rad = 3; // set measurement dimension, radar can measure r, phi, and r_dot

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_rad, 2 * n_aug + 1); // 1 column for each sigma point 

  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement in measurement space 
  // Predicted Mean : zk+1∣k = ∑i=1 to nσ​  (wi * Zk+1∣k,i)
  VectorXd z_pred = VectorXd(n_rad);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) 
  {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  ///* Sk+1∣k = ∑i=1 to ​nσ (wi * (Zk+1∣k,i − zk+1∣k) * (Zk+1∣k,i − zk+1∣k)T  +  R;
  /***** we dont use augmentation here because the noise R is just additive part of the funtion ******/
  MatrixXd S = MatrixXd(n_rad,n_rad); 
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_rad,n_rad);
  R <<    std_radr*std_radr, 0, 0,
          0, std_radphi*std_radphi, 0,
          0, 0,std_radrd*std_radrd;
  S = S + R;

  /***************************************************/
  // Updating 
  /***************************************************/

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_rad);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights(i) * x_diff * z_diff.transpose();


    //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  }
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  H_laser_ << 1,0,0,0,
              0,1,0,0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  ekf_.F_ = MatrixXd::Identity(4,4);

  ekf_.P_ = MatrixXd::Identity(4,4);

  ekf_.Q_ = MatrixXd::Zero(4,4);




}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "Initializing EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    cout << "measurement_pack: " <<  measurement_pack.sensor_type_ << measurement_pack.raw_measurements_ << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        auto ro = measurement_pack.raw_measurements_(0);
        auto phi = measurement_pack.raw_measurements_(1);
        auto ro_dot = measurement_pack.raw_measurements_(2);
        auto sin_phi = sin(phi);
        auto cos_phi = cos(phi);

        ekf_.x_ << ro * cos_phi,  ro * sin_phi, 0, 0; //ro_dot*cos_phi, ro_dot*sin_phi;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
        ekf_.P_(2,2) = 100;
        ekf_.P_(3,3) = 100;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "KF init done!" << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  long long delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0,2) = delta_t;
  ekf_.F_(1,3) = delta_t;


  ekf_.Q_(0,0) = pow(delta_t,4)/4*ekf_.noise_ax;
  ekf_.Q_(0,2) = pow(delta_t,3)/2*ekf_.noise_ax;
  ekf_.Q_(1,1) = pow(delta_t,4)/4*ekf_.noise_ay;
  ekf_.Q_(1,3) = pow(delta_t,3)/2*ekf_.noise_ay;
  ekf_.Q_(2,0) = ekf_.Q_(0,2);
  ekf_.Q_(2,2) = pow(delta_t,2)*ekf_.noise_ax;
  ekf_.Q_(3,1) = ekf_.Q_(1,3);
  ekf_.Q_(3,3) = pow(delta_t,2)*ekf_.noise_ay;

  cout << "Predict step" << endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  cout << "Update step" << endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      Hj_ = tools.CalculateJacobian(ekf_.x_);

      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

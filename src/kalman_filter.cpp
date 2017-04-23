#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
    this->x_ = this->F_*this->x_;
    this->P_ = this->F_*this->P_*this->F_.transpose() + this->Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    MatrixXd y = z - this->H_ * this->x_;
    MatrixXd S = this->H_ * this->P_ * this->H_.transpose() + this->R_;
    MatrixXd K = this->P_*this->H_.transpose() * S.inverse();
    this->x_ = this->x_ + (K*y);
    auto rows = this->P_.rows();
    auto cols = this->P_.cols();
    MatrixXd I = MatrixXd::Identity(rows, cols);
    this->P_ = (I - K*this->H_) * this->P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}

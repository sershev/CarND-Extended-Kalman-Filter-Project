#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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

void KalmanFilter::UpdateValues(const MatrixXd &y){
    MatrixXd PH_transpose = P_ * H_.transpose();
    MatrixXd S = H_ * PH_transpose + R_;
    MatrixXd K = PH_transpose * S.inverse();
    x_ = x_ + (K*y);
    int rows = P_.rows();
    int cols = P_.cols();
    MatrixXd I = MatrixXd::Identity(rows, cols);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::Predict() {
    x_ = F_*x_;
    P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    MatrixXd y = z - H_ * x_;
    KalmanFilter::UpdateValues(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    VectorXd y = z - CartesianToPolar();
    if (y(1) > M_PI){
        y(1) = y(1) - pi2;
    }
    else if(y(1) < -M_PI){
        y(1) = y(1) + pi2;
    }
    KalmanFilter::UpdateValues(y);
}

VectorXd KalmanFilter::CartesianToPolar(){
    VectorXd polar(3);
    double phi;
    auto ro = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    if (x_(0) == 0){
        phi = phi_fallback;
    }else{
        phi = atan2(x_(1), x_(0));
    }
    auto ro_dot = (x_(0)*x_(2)+ x_(1)*x_(3)) / std::max(ro, ro_fallback);

    cout << "ro: " << ro << endl;
    cout << "phi: " << phi << endl;
    cout << "ro_dot: " << ro_dot << endl;
    polar << ro, phi, ro_dot;
    return polar;
}


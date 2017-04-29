#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd mean(4);
    mean << 0,0,0,0;
    int n = estimations.size();
    assert (n > 0);
    assert (n == ground_truth.size());

    for (int i=0; i < n; ++i){
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array().pow(2);
        mean += residual;
    }
    mean = (mean / n).array().sqrt();
    return mean;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    float pxpy_squared = pow(px,2) + pow(py,2);
    float pxpy_squared_square_root = sqrt(pxpy_squared);
    float vypx = vy*px;
    float vxpy = vx*py;

    //check division by zero
    if(px==0 || py==0){
        cerr << "Tools::CalculateJacobian() - Error - Devision by Zero" << endl;
    }

    auto px_by_pxpy_squared_square_root = px/pxpy_squared_square_root;
    auto py_by_pxpy_squared_square_root = py/pxpy_squared_square_root;

    auto pdot_nach_py = px*(vypx-vxpy)/pow(pxpy_squared_square_root, 3);
    auto pdot_nach_px = py*(vxpy-vypx)/pow(pxpy_squared_square_root, 3);

    Hj << px_by_pxpy_squared_square_root, py_by_pxpy_squared_square_root, 0,0,
            -py/pxpy_squared, px/pxpy_squared, 0,0,
            pdot_nach_px, pdot_nach_py, px_by_pxpy_squared_square_root, py_by_pxpy_squared_square_root;
    return Hj;
}


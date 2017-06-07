#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  /**
    * Calculate the RMSE here.
  */

  // Input validation
  assert(estimations.size()>0);
  assert(estimations.size()==ground_truth.size());

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  //std::cout << "RMSE: " << rmse << std::endl;

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd diff(4);
    diff = (estimations[i]-ground_truth[i]).array().pow(2);
    //std::cout << "GROUNDTRUTH: " << ground_truth[i] << std::endl;
    //std::cout << "ESTIMATIONS: " << estimations[i] << std::endl;
    //std::cout << "diff: " << diff << std::endl;
    rmse += diff;
    float check = rmse(1) / estimations.size();
    // std::cout << sqrt(check) << std::endl;
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;

}

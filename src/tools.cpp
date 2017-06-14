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
      * Calculate the RMSE
  */
  VectorXd rmse(4);
  rmse.fill(0);

  //  * the estimation vector size should not be zero
  if(estimations.size() == 0){
    std::cout << "Estimation vector size is zero" << std::endl;
    return rmse;
  }
    
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == ground_truth.size()){
	VectorXd temp(4);

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
    		temp = estimations[i]-ground_truth[i];
    		temp = temp.array()*temp.array();
    		rmse += temp;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();
  }
  else
    	std::cout << "Estimation vector and ground truth vector are of different size" << std::endl;

  //return the result
  return rmse;
}

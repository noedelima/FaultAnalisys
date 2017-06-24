/*
 * faultanalysis.hpp
 *
 *  Created on: 23 de jun de 2017
 *      Author: noedelima
 */

#ifndef SRC_FAULTANALYSIS_HPP_
#define SRC_FAULTANALYSIS_HPP_

#include <iostream>
#include <Eigen/Eigen>
using namespace std;
using namespace Eigen;

MatrixXcd Ifault(MatrixXcd Ybus, string genpath, int bus, complex<double> Zf);

#endif /* SRC_FAULTANALYSIS_HPP_ */

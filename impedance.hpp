/*
 * impedance.hpp
 *
 *  Created on: 17 de jun de 2017
 *      Author: noedelima
 */

#ifndef SRC_IMPEDANCE_HPP_
#define SRC_IMPEDANCE_HPP_
#include "DataStruct/database.hpp"
#include "network.hpp"
#include <Eigen/Eigen>
using namespace Eigen;

MatrixXcd Ythev(MatrixXcd Y, string genpath);
MatrixXcd kron(MatrixXcd A);
VectorXcd Zthev(MatrixXcd Ythevenin, int n);

#endif /* SRC_IMPEDANCE_HPP_ */

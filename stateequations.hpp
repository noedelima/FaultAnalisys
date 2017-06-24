/*
 * stateequations.hpp
 *
 *  Created on: 21 de jun de 2017
 *      Author: noedelima
 */

#ifndef SRC_STATEEQUATIONS_HPP_
#define SRC_STATEEQUATIONS_HPP_
#include <iostream>
#include <Eigen/Eigen>
#include <complex>
using namespace std;
using namespace Eigen;

VectorXcd Voltage(MatrixXcd Ybus);
MatrixXcd State(MatrixXcd Ybus, VectorXcd V);
MatrixXcd Current(MatrixXcd Ybus, VectorXcd V);

Vector3cd Symmetrical(Vector3cd X);
Vector3cd backSymmetrical(Vector3cd X);

complex<double> dPijdTk(MatrixXcd Y, VectorXcd V, int i, int j, int k);
complex<double> dQijdTk(MatrixXcd Y, VectorXcd V, int i, int j, int k);
complex<double> dPijdVk(MatrixXcd Y, VectorXcd V, int i, int j, int k);
complex<double> dQijdVk(MatrixXcd Y, VectorXcd V, int i, int j, int k);
complex<double> dVidTj(MatrixXcd Y, VectorXcd V, int i, int j);
complex<double> dVidVj(MatrixXcd Y, VectorXcd V, int i, int j);

#endif /* SRC_STATEEQUATIONS_HPP_ */

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

void printmatrix(MatrixXcd M);
VectorXcd Vfaultsym(MatrixXcd Ybus, string genpath, int bus, complex<double> Zf);
Matrix<Matrix3cd, Dynamic, Dynamic> ProdMatrixOverload(Matrix<Matrix3cd, Dynamic, Dynamic> A, Matrix<Matrix3cd, Dynamic, Dynamic> B);

#endif /* SRC_FAULTANALYSIS_HPP_ */

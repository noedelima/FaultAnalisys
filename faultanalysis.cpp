/*
 * faultanalysis.cpp
 *
 *  Created on: 23 de jun de 2017
 *      Author: noedelima
 */


#include "faultanalysis.hpp"
#include "stateequations.hpp"
#include "impedance.hpp"

#define pi 3.14159265358979;

void printmatrix(MatrixXcd M){
	for (int i=0; i<M.rows(); i++){
		for (int j=0; j<M.cols(); j++){
			cout << abs(M(i,j)) << " /_ " << 180*arg(M(i,j))/pi;
			cout << "ºpu	";
		}
		cout << endl;
	}
	cout << endl;
}

Matrix<Matrix3cd, Dynamic, Dynamic> ProdMatrixOverload(Matrix<Matrix3cd, Dynamic, Dynamic> A, Matrix<Matrix3cd, Dynamic, Dynamic> B){
	Matrix<Matrix3cd, Dynamic, Dynamic> P(A.rows(), B.cols());
	if (A.cols() != B.rows()){return P;}
	for (int i=0; i<A.rows(); i++){
		for (int j=0; j<B.cols(); j++){
			P(i,j) = Matrix3cd::Zero(3,3);
			for (int k=0; k<A.cols(); k++){
				P(i,j) += A(i,k)*B(k,j);
			}
		}
	}
	return P;
} // Multiply two overloaded matrix (a matrix of matrix)

VectorXcd Vfaultsym(MatrixXcd Ybus, string genpath, int bus, complex<double> Zf){
	int ref = bus - 1;
	VectorXcd Vbus = Voltage(Ybus);
	MatrixXcd Yth = Ythev(Ybus, genpath);
	VectorXcd Zth = Zthev(Yth, ref);
	complex<double> Ifault = -Vbus(ref)/(Zth(ref)+Zf);
	VectorXcd dV = Ifault*Zth;
	VectorXcd V = Vbus + dV;
	cout << "A corrente de falta na barra " << bus << " é:" << endl << abs(Ifault) << " /_" << 180*arg(Ifault)/3.14159265358979 << "º" << endl << endl;
	return V;
}

VectorXcd VfaultAsym(MatrixXcd Ybus, string genpath, int bus, complex<double> Zf){
	int ref = bus - 1;
	Matrix<Matrix3cd, Dynamic, Dynamic> Ysim(Ybus.rows(), Ybus.cols());
	VectorXcd Vbus = Voltage(Ybus);
	MatrixXcd Yth = Ythev(Ybus, genpath);
	VectorXcd Zth = Zthev(Yth, ref);
	complex<double> Ifault = -Vbus(ref)/(Zth(ref)+Zf);
	VectorXcd dV = Ifault*Zth;
	VectorXcd V = Vbus + dV;
	cout << "A corrente de falta é:" << endl << abs(Ifault) << " /_" << 180*arg(Ifault)/3.14159265358979 << "º" << endl << endl;
	return V;
}

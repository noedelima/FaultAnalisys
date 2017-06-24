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

MatrixXcd Ifault(MatrixXcd Ybus, string genpath, int bus, complex<double> Zf){
	int ref = bus - 1;
	VectorXcd If = VectorXcd::Zero(Ybus.rows());
	VectorXcd Vbus = Voltage(Ybus);
	MatrixXcd Yth = Ythev(Ybus, genpath);
	VectorXcd Zth = Zthev(Yth, ref);
	complex<double> Ifault = -Vbus(ref)/(Zth(ref)+Zf);
	If(ref) = Ifault;
	VectorXcd dV = Ifault*Zth;
	VectorXcd V = Vbus + dV;
	MatrixXcd Iflow = MatrixXcd::Zero(Ybus.rows(),Ybus.cols());
	cout << "A corrente de falta é:" << endl << abs(Ifault) << " /_" << 180*arg(Ifault)/3.14159 << "º" << endl << endl;
	//cout << "O vetor de variação das tensões é:" << endl;
	//printmatrix(dV);
	cout << "O vetor de tensões de falta é:" << endl;
	printmatrix(V);
	cout << endl << endl;
	cout << "Assim, temos as seguintes correntes de falta:" << endl;
	for (int i=0; i<Ybus.rows(); i++){
		for (int j=i+1; j<Ybus.cols(); j++){
			if (abs(Yth(i,j))>0.001){
				Iflow(i,j) = Yth(i,j)*(V(i)-V(j));
				Iflow(j,i) = Iflow(i,j);
				cout << "A variação de Icc da linha " << i+1 << " para a linha " << j+1 << " é "
						<< abs(Iflow(i,j)) << " /_ " << 180*arg(Iflow(i,j))/pi;
				cout << "º pu" << endl;
			}
		}
	}
	cout << endl << endl;
	return Iflow;
}

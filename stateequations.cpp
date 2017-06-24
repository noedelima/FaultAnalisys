/*
 * stateequations.cpp
 *
 *  Created on: 21 de jun de 2017
 *      Author: noedelima
 */

#include "stateequations.hpp"

/*
 * Here is the state equations and functions to load flow
 * Voltage initializes voltage at all buses as 1 pu and 0 degree
 * State uses voltage vector and bus admitance matrix to calculate power flow
 * Current uses voltage vector and bus admitance matrix to calculate current flow
 */
VectorXcd Voltage(MatrixXcd Ybus){
	VectorXcd V = VectorXcd(Ybus.rows());
	for (int i=0; i<Ybus.rows(); i++){V(i) = 1;}
	return V;
} // Initialize voltage vector with 1 pu magnitude and 0 phase angle

MatrixXcd State(MatrixXcd Ybus, VectorXcd V){
	MatrixXcd Power = MatrixXcd::Zero(Ybus.rows(),Ybus.cols());
	for (int i=0; i<Ybus.rows(); i++){
		for (int j=0; j<Ybus.cols(); j++){
			double G = Ybus(i,j).real();
			double B = Ybus(i,j).imag();
			double gsh = 0; // shunt real
			double bsh = 0; // shunt imag
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			double Tij = arg(V(i)) - arg(V(j));
			Power(i,i).real( Power(i,i).real() + Vi*Vj*(G*cos(Tij)+B*sin(Tij)));
			Power(i,i).imag(Power(i,i).imag() + Vi*Vj*(G*sin(Tij)-B*cos(Tij)));
			if(i!=j){
				Power(i,j).real(pow(Vi,2)*(G+gsh) - Vi*Vj*(G*cos(Tij)+B*sin(Tij)));
				Power(i,j).imag(-pow(Vi,2)*(B+bsh) - Vi*Vj*(G*sin(Tij)-B*cos(Tij)));
			}
		}
	}
	return Power;
} // Actualize Power Flow from Ybus and state vector of voltages

MatrixXcd Current(MatrixXcd Ybus, VectorXcd V){
	MatrixXcd Current = MatrixXcd::Zero(Ybus.rows(),Ybus.cols());
	for (int i=0; i<Ybus.rows(); i++){
		for (int j=0; j<Ybus.cols(); j++){
			double G = Ybus(i,j).real();
			double B = Ybus(i,j).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			double Tij = arg(V(i)) - arg(V(j));
			Current(i,j) = pow((pow(G,2)+pow(B,2))*(pow(Vi,2)+pow(Vj,2)-2*Vi*Vj*cos(Tij)), 0.5);
		}
	}
	for (int i=0; i<Ybus.rows(); i++){
		for (int j=0; j<Ybus.cols(); j++){
			if (i!=j){Current(i,i) -= Current(i,j);}
		}
	}
	return Current;
} // Actualize Current Flow from Ybus and state vector of voltages

Vector3cd Symmetrical(Vector3cd X){
	Vector3cd Y = VectorXcd::Zero(X.rows());
	complex<double> a(-0.5,sqrt(3)/2); // vector alpha defined as e^(j120degree)
	Matrix3cd A; // Symmetrical components transformation matrix
	A << 	1, 	1,		 	1,
			1, 	pow(a,2), 	a,
			1, 	a,		 	pow(a,2);
	Y = A*X;
	return Y;
} // Act an symmetrical components transform

Vector3cd backSymmetrical(Vector3cd X){
	Vector3cd Y = VectorXcd::Zero(X.rows());
	complex<double> a(-0.5,sqrt(3)/2); // vector alpha defined as e^(j120degree)
	Matrix3cd A; // Symmetrical components transformation matrix
	A << 	1, 	1,		 	1,
			1, 	pow(a,2), 	a,
			1, 	a,		 	pow(a,2);
	Y = A.inverse()*X;
	return Y;
} // Act an inverse symmetrical components transform

/*
 * Here is the derive functions of state equations to power system state estimator
 * The complete theory can be found in Power System State Estimator, Abur and Expósito
 */
complex<double> dPijdTk(MatrixXcd Ybus, VectorXcd V, int i, int j, int k){
	complex<double> dPijdTk = 0;
	if (i == j){
		if (i == k){
			double Vi = abs(V(i));
			for (int t = 0; t < Ybus.rows(); t++){
				double G = Ybus(i,t).real();
				double B = Ybus(i,t).imag();
				double Tit = arg(V(i)) - arg(V(t));
				double Vt = abs(V(t));
				dPijdTk += Vt*(-G*sin(Tit)+B*cos(Tit));
			}
			dPijdTk = dPijdTk*Vi - pow(Vi,2)*Ybus(i,i).imag();
		} // dPi/dTi
		else{
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			double Tik = arg(V(i)) - arg(V(k));
			double Vi = abs(V(i));
			double Vk = abs(V(k));
			dPijdTk = Vi*Vk*(G*sin(Tik)-B*cos(Tik));
		} // dPi/dTj
	} // Bus Real Power
	else{
		if (i == k){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dPijdTk = Vi*Vj*(G*sin(Tij)-B*cos(Tij));
		} // dPij/dTi
		else if(k == j){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dPijdTk = -Vi*Vj*(G*sin(Tij)-B*cos(Tij));
		} // dPij/dTj
	} // Line Injection of Real Power
	return dPijdTk;
} //derive real power from i bus to j bus relative to theta voltage angle in k bus

complex<double> dQijdTk(MatrixXcd Ybus, VectorXcd V, int i, int j, int k){
	complex<double> dQijdTk = 0;
	if (i == j){
		if (i == k){
			double Vi = abs(V(i));
			for (int t = 0; t < Ybus.rows(); t++){
				double G = Ybus(i,t).real();
				double B = Ybus(i,t).imag();
				double Tit = arg(V(i)) - arg(V(t));
				double Vt = abs(V(t));
				dQijdTk += Vt*(G*cos(Tit)+B*sin(Tit));
			}
			dQijdTk = dQijdTk*Vi - pow(Vi,2)*Ybus(i,i).real();
		} // dQi/dTi
		else{
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			double Tik = arg(V(i)) - arg(V(k));
			double Vi = abs(V(i));
			double Vk = abs(V(k));
			dQijdTk = Vi*Vk*(-G*cos(Tik)-B*sin(Tik));
		} // dQi/dTj
	} // Bus Real Power
	else{
		if (i == k){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dQijdTk = -Vi*Vj*(G*cos(Tij)+B*sin(Tij));
		} // dQij/dTi
		else if(k == j){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dQijdTk = Vi*Vj*(G*cos(Tij)+B*sin(Tij));
		} // dQij/dTj
	} // Line Injection of Real Power
	return dQijdTk;
} //derive reactive power from i bus to j bus relative to theta voltage angle in k bus

complex<double> dPijdVk(MatrixXcd Ybus, VectorXcd V, int i, int j, int k){
	complex<double> dPijdVk = 0;
	if (i == j){
		if (i == k){
			double Vi = abs(V(i));
			for (int t = 0; t < Ybus.rows(); t++){
				double G = Ybus(i,t).real();
				double B = Ybus(i,t).imag();
				double Tit = arg(V(i)) - arg(V(t));
				double Vt = abs(V(t));
				dPijdVk += Vt*(G*cos(Tit)+B*sin(Tit));
			}
			dPijdVk += Vi*Ybus(i,i).real();
		} // dPi/dVi
		else{
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			double Tik = arg(V(i)) - arg(V(k));
			double Vi = abs(V(i));
			double Vk = abs(V(k));
			dPijdVk = Vi*(G*cos(Tik)+B*sin(Tik));
		} // dPi/dVj
	} // Bus Real Power
	else{
		if (i == k){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dPijdVk = -Vj*(G*cos(Tij)+B*sin(Tij)) + 2*Vi*(G+gsh);
		} // dPij/dVi
		else if(k == j){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dPijdVk = -Vi*(G*cos(Tij)+B*sin(Tij));
		} // dPij/dVj
	} // Line Injection of Real Power
	return dPijdVk;
} //derive active power from i bus to j bus relative to voltage magnitude in k bus

complex<double> dQijdVk(MatrixXcd Ybus, VectorXcd V, int i, int j, int k){
	complex<double> dQijdVk = 0;
	if (i == j){
		if (i == k){
			double Vi = abs(V(i));
			for (int t = 0; t < Ybus.rows(); t++){
				double G = Ybus(i,t).real();
				double B = Ybus(i,t).imag();
				double Tit = arg(V(i)) - arg(V(t));
				double Vt = abs(V(t));
				dQijdVk += Vt*(G*sin(Tit)-B*cos(Tit));
			}
			dQijdVk -= Vi*Ybus(i,i).imag();
		} // dQi/dVi
		else{
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			double Tik = arg(V(i)) - arg(V(k));
			double Vi = abs(V(i));
			double Vk = abs(V(k));
			dQijdVk = Vi*(G*sin(Tik)-B*cos(Tik));
		} // dQi/dVj
	} // Bus Real Power
	else{
		if (i == k){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dQijdVk = -Vi*(G*sin(Tij)-B*cos(Tij)) - 2*Vi*(B+bsh);
		} // dQij/dVi
		else if(k == j){
			double G = Ybus(i,k).real();
			double B = Ybus(i,k).imag();
			//double gsh = 0; // shunt real
			//double bsh = 0; // shunt imag
			double Tij = arg(V(i)) - arg(V(j));
			double Vi = abs(V(i));
			double Vj = abs(V(j));
			dQijdVk = -Vi*(G*sin(Tij)-B*cos(Tij));
		} // dQij/dVj
	} // Line Injection of Real Power
	return dQijdVk;
} //derive reactive power from i bus to j bus relative to voltage magnitude in k bus

complex<double> dVidTj(MatrixXcd Y, VectorXcd V, int i, int j){
	complex<double> dVidTj = 0;
	return dVidTj;
} //derive voltage magnitude in i bus relative to voltage phase angle in j bus

complex<double> dVidVj(MatrixXcd Y, VectorXcd V, int i, int j){
	complex<double> dVidVj = 0;
	if (i == j){
		dVidVj = 1;
	}
	return dVidVj;
} //derive voltage magnitude in i bus relative to voltage magnitude in j bus


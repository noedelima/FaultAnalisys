/*
 * impedance.cpp
 *
 *  Created on: 17 de jun de 2017
 *      Author: noedelima
 */

#include "impedance.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "topology.hpp"
using namespace std;
using namespace boost;

MatrixXcd Ythev(MatrixXcd Y, string genpath){
	queue<generator> generating = generators(genpath);
	while(generating.size()>0){
		generator actual = *generating.pop();
		int i = actual.id-1;
		Y(i,i) += pow(actual.Zd, -1);
	}
	return Y;
} // Include to reference connections in Ybus Matrix

MatrixXcd kron(MatrixXcd A){
	int col = A.cols();
	int lin = A.rows();
	MatrixXcd A11, A12, A21;
	A11 = A.block(0, 0, lin-1, col-1);
	A12 = A.block(0, col-1, lin-1, 1);
	A21 = A.block(lin-1, 0, 1, col-1);
	complex<double> A22 = A(lin-1, col-1);
	A11 -= pow(A22, -1)*(A12*A21);
	return A11;
} // A Kron Reduction as part of Gauss method

VectorXcd Zthev(MatrixXcd Ythevenin, int n){
	int dim = Ythevenin.cols();
	VectorXcd I = VectorXcd::Zero(dim);
	VectorXcd Zth = VectorXcd::Zero(dim);
	I(n) = 1;
	Zth = Ythevenin.fullPivLu().solve(I);
	return Zth;
} // Solve a column of Zbus using LU decomposition of Ybus

MatrixXcd Zt1(MatrixXcd Z, generator G){
	int dim = Z.cols();
	//int i = G.id - 1; // Don't matter
	complex<double> zb = G.Zd;
	MatrixXcd M = MatrixXcd::Zero(dim+1, dim+1);
	M.block(0,0,Z.rows(),Z.cols()) << Z;
	M(dim,dim) = zb;
	return M;
} // New connection from a new bus to the reference bus

MatrixXcd Zt2(MatrixXcd Z, branch B){
	int dim = Z.cols();
	int i = B.id_fr - 1;
	//int j = B.id_to - 1; // Don't matter
	complex<double> zb = B.serieimpedance;
	complex<double> zsh = pow(B.halfshuntadmitance,-1);
	MatrixXcd M = MatrixXcd::Zero(dim+1, dim+1);
	M.block(0,0,Z.rows(),Z.cols()) << Z;
	M.block(0,dim,dim,1) << Z.block(0,i,dim,1);
	M.block(dim,0,1,dim) << Z.block(i,0,1,dim);
	M(i,i) += zsh;
	M(dim,dim) = Z(i,i) + zb + zsh;
	return M;
} // New connection from new bus to an existent bus

MatrixXcd Zt3(MatrixXcd Z, generator G){
	int dim = Z.cols();
	int i = G.id - 1;
	complex<double> zb = G.Zd;
	MatrixXcd M = MatrixXcd::Zero(dim+1, dim+1);
	M.block(0,0,Z.rows(),Z.cols()) << Z;
	M.block(0,dim,dim,1) << Z.block(0,i,dim,1);
	M.block(dim,0,1,dim) << Z.block(i,0,1,dim);
	M(dim,dim) = Z(i,i) + zb;
	Z = kron(M);
	return Z;
} // New connection from an existent bus to the reference bus

MatrixXcd Zt4(MatrixXcd Z, branch B){
	int dim = Z.cols();
	int i = B.id_fr - 1;
	int j = B.id_to - 1;
	complex<double> zb = B.serieimpedance;
	complex<double> zsh = pow(B.halfshuntadmitance,-1);
	MatrixXcd M = MatrixXcd::Zero(dim+1, dim+1);
	M.block(0,0,Z.rows(),Z.cols()) << Z;
	M.block(0,dim,dim,1) << (Z.block(0,i,dim,1) - Z.block(0,j,dim,1));
	M.block(dim,0,1,dim) << (Z.block(i,0,1,dim) - Z.block(j,0,1,dim));
	M(dim,dim) = Z(i,i) + Z(j,j) - Z(i,j) - Z(j,i) + zb;
	double tol = 0.1; // Minimal parallel admitance to include
	if (zsh.real() > tol || zsh.imag() > tol){
		Z = kron(M);
		M.block(0,0,Z.rows(),Z.cols()) << Z;
		M.block(0,dim,dim,1) << Z.block(0,i,dim,1);
		M.block(dim,0,1,dim) << Z.block(i,0,1,dim);
		M(dim,dim) = Z(i,i) + zsh;
		Z = kron(M);
		M.block(0,0,Z.rows(),Z.cols()) << Z;
		M.block(0,dim,dim,1) << Z.block(0,j,dim,1);
		M.block(dim,0,1,dim) << Z.block(j,0,1,dim);
		M(dim,dim) = Z(j,j) + zsh;
	}
	Z = kron(M);
	return Z;
} // New connection from an existent bus to another existent bus

MatrixXcd Zbus(string branchpath, string genpath){
	double inf_imped = pow(10,12);
	database<bus, int> buslist = buses(branches(branchpath));
	int dim = buslist.size();
	MatrixXcd Z = inf_imped*(MatrixXcd::Identity(dim, dim));
	queue<generator> generating = generators(genpath);
	while(generating.size()>0){
		generator actual = *generating.pop();
		Z = Zt3(Z, actual);
	}
	queue<branch> connections = branches(branchpath);
	while(connections.size()>0){
		branch actual = *connections.pop();
		Z = Zt4(Z, actual);
	}
	return Z;
}

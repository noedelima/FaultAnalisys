/*
 * topology.cpp
 *
 *  Created on: 28 de mai de 2017
 *      Author: noedelima
 */

#include "topology.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;


queue<branch> branches(string netname){
	string output;
	queue<branch> branches;
	ifstream network(netname);
	if (!network){
		cout << "Falha na abertura do arquivo " << netname << endl;
		return branches;
	}
	cout << "Documento " << netname << " aberto com sucesso!" << endl;
	int i=100, line=0; // Max char to read in each line, and a line counter
	while(!network.eof()){
		line++;
		output = "";
		char buffer[i];
		vector<string> vecsplit;
		network.getline(buffer, i);
		if (line==1 || network.eof()){continue;}
		output += buffer;
		split(vecsplit, output, is_any_of(","));
		string tag = vecsplit[0];
		int i, j;
		double Rserie=0, Xserie=0, Bshunt=0;
		i =  lexical_cast<int>(vecsplit[1]);
		j =  lexical_cast<int>(vecsplit[2]);
		if (i != j){
			Rserie = lexical_cast<double>(vecsplit[3]);
			Xserie = lexical_cast<double>(vecsplit[4]);
			Bshunt = lexical_cast<double>(vecsplit[5]);
			branch* connection = new branch(tag, i, j, Rserie, Xserie, Bshunt);
			branches.push(connection);
		}
	}
	network.close();
	return branches;
}

queue<generator> generators(string genname){
	string output;
	queue<generator> gen;
	ifstream network(genname);
	if (!network){
		cout << "Falha na abertura do arquivo " << genname << endl;
		return gen;
	}
	cout << "Documento " << genname << " aberto com sucesso!" << endl;
	int i=100, line=0; // Max char to read in each line, and a line counter
	while(!network.eof()){
		line++;
		output = "";
		char buffer[i];
		vector<string> vecsplit;
		network.getline(buffer, i);
		if (line==1 || network.eof()){continue;}
		output += buffer;
		split(vecsplit, output, is_any_of(","));
		string tag = vecsplit[0];
		int loc =  lexical_cast<int>(vecsplit[1]);
		double Xserie, Rserie=0; //, Xserie=0, Bshunt=0;
		Rserie = lexical_cast<double>(vecsplit[2]);
		Xserie = lexical_cast<double>(vecsplit[3]);
		//Bshunt = lexical_cast<double>(vecsplit[4]);
		generator* connection = new generator(tag, loc, Rserie, Xserie);
		gen.push(connection);
	}
	network.close();
	return gen;
}

database<bus, int> buses(queue<branch> branches){
	database<bus, int> buslist;
	while(branches.size()){
		branch* temp = branches.pop();
		if(!buslist.find(temp->id_fr)){
			bus* include = new bus("bus", temp->id_fr);
			buslist.push(include, temp->id_fr);
		}
		if(!buslist.find(temp->id_to)){
			bus* include = new bus("bus", temp->id_to);
			buslist.push(include, temp->id_to);
		}
	}
	return buslist;
}

MatrixXcd Ybus(string branchpath){
	database<bus, int> buslist = buses(branches(branchpath));
	int dim = buslist.size();
	MatrixXcd Y = MatrixXcd::Zero(dim, dim);
	queue<branch> connections = branches(branchpath);
	while(connections.size()>0){
		branch actual = *connections.pop();
		int i = actual.id_fr-1;
		int j = actual.id_to-1;
		complex<double> yserie = pow(actual.serieimpedance, -1);
		complex<double> yshunt = actual.halfshuntadmitance;
		Y(i,j) -= yserie;
		Y(j,i) -= yserie;
		Y(i,i) += (yserie + yshunt);
		Y(j,j) += (yserie + yshunt);
	}
	return Y;
}

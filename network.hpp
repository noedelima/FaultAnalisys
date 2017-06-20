//============================================================================
// Name        : network.hpp
// Author      : Noe de Lima Bezerra
// Date        : may, 17, 2017
// Version     : 0.0
// Copyright   : Public License v. 2.0.
// Description : Power System network elements
//============================================================================

#ifndef NETWORK_HPP_
#define NETWORK_HPP_

#include <iostream>
#include <string>
#include <complex>
#include "DataStruct/database.hpp"
using namespace std;

class voltagemeas{
public:
	string tag;
	int id;
	double measure;
	voltagemeas(string get_tag, int locate, double value){
		tag = get_tag;
		id = locate;
		measure = value;
	}
};

class powermeas{
public:
	string tag;
	int id_fr;
	int id_to;
	complex<double> measure;
	string type;
	powermeas(string get_tag, int get_from, int get_to, double P, double Q){
		tag = get_tag;
		id_fr = get_from;
		id_to = get_to;
		measure = complex<double>(P, Q);
		if (id_fr == id_to){
			type = "power injection";
		}
		else{
			type = "power flow";
		}
	}
};

class generator{
public:
	string tag;
	int id;
	complex<double> Zd;
	generator(string get_tag, int locate, double get_Rd, double get_Xd){
		tag = get_tag;
		id = locate;
		Zd.real(get_Rd);
		Zd.imag(get_Xd);
	}
};

class load{
public:
	string tag;
	int id;
	complex<double> Yload;
	complex<double> power;
	load(string get_tag, int locate, complex<double> Yl=0, complex<double> S=0){
		tag = get_tag;
		id = locate;
		Yload = Yl;
		power = S;
	}
};

class bus{
public:
	string tag;
	int id;
	queue<load> staticload;
	queue<generator> Gen;
	string type;
	bus(){
		tag = "";
		id = 0;
		type = "";
	}
	bus(string get_tag, int locate, string bustype = ""){
		tag = get_tag;
		id = locate;
		type = bustype;
	}
};

class branch{
public:
	string tag;
	int id_fr;
	int id_to;
	complex<double> serieimpedance;
	complex<double> halfshuntadmitance;
	branch(){
		tag = "";
		id_fr = 0;
		id_to = 0;
		serieimpedance = complex<double>(0, 0);
		halfshuntadmitance = complex<double>(0, 0);
	}
	branch(string get_tag, int get_from, int get_to, double Rserie, double Xserie, double Bshunt){
		tag = get_tag;
		id_fr = get_from;
		id_to = get_to;
		serieimpedance = complex<double>(Rserie, Xserie);
		halfshuntadmitance = complex<double>(0, Bshunt);
	}
	void Zserie(complex<double> Z){
		serieimpedance = Z;
	}
	void Yshunt(complex<double> Y){
		halfshuntadmitance = Y;
	}
};


#endif /* NETWORK_HPP_ */

//============================================================================
// Name        : main.cpp
// Author      : Noe de Lima Bezerra
// Date        : may, 17, 2017
// Version     : 0.0
// Copyright   : Public License v. 2.0.
// Description : Power System Analisys
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include "topology.hpp"
#include "impedance.hpp"
#include "stateequations.hpp"
#include "faultanalysis.hpp"
using namespace std;

int main() {
	string net = "2"; // Choose system to run
	cout << "!!!Hello World!!!" << endl;
	IOFormat HeavyFmt(FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
	IOFormat matlabFmt(FullPrecision, 0, ", ", "; ", "[", "]", "", "");
	string path, extension, netname, genname;
	path = "/home/noedelima/Documentos/11\ -\ Arquivos\ e\ Utilitários\ Diversos/workspace/PowerSystemAnalisys/src/data/";
	extension = ".csv";
	ofstream fileout(path+"output.txt", ofstream::app);
	fileout << "Início da execução do programa: " << endl << __DATE__ << " >> " << __TIME__ << endl;
	cout << "Qual o nome do arquivo da rede? " << endl;
	netname = "network" + net; //cin >> netname;
	cout << "Qual o nome do arquivo de geração e carga? " << endl;
	genname = "load_gen" + net; //cin >> genname;
	queue<branch> connections = branches(path+netname+extension);
	fileout << "A lista de conexões entre barras tem " << connections.size() << " elementos:" << endl;
	while(connections.size()>0){
		branch actual = *connections.pop();
		fileout << "Conexão da barra " << actual.id_fr << " para a barra " << actual.id_to << " identificada como " << actual.tag << endl;
	}
	queue<generator> generating = generators(path+genname+extension);
	fileout << "A lista de conexões com a referência tem " << generating.size() << " elementos:" << endl;
	while(generating.size()>0){
		generator actual = *generating.pop();
		fileout << "Conexão na barra " << actual.id << " identificado como " << actual.tag << endl;
	}
	MatrixXcd Yb = Ybus(path+netname+extension);
	MatrixXcd Yth = Ythev(Yb, path+genname+extension);
	fileout << "A matriz Ybus da rede é:" << endl << endl << Yb.format(HeavyFmt) << endl << endl;
	fileout << "A matriz  de admitância equivalente Thevenin da rede é:" << endl << endl << Yth.format(HeavyFmt) << endl << endl;
	fileout << "A Zbarra obtida por inversão é:" << endl << endl << Yth.inverse().format(HeavyFmt) << endl << endl;
	fileout << "Aa colunas ímpares da Zbarra são:" << endl;
	for (int i=0; i<Yth.cols(); i+=2){
		fileout << "Coluna " << i+1 << ":";
		fileout << endl << Zthev(Yth, i).format(HeavyFmt) << endl;
	}
	fileout << endl;
	MatrixXcd Zbdir = Zbus(path+netname+extension, path+genname+extension);
	fileout << "Matriz Zbus por montagem direta:" << endl << Zbdir.format(HeavyFmt) << endl;
	fileout << "Fim da execução do programa: " << endl << __DATE__ << " >> " << __TIME__ << endl;
	fileout.close();
	cout << "Continuação do programa" << endl;
	cout << endl;
	Ifault(Yb, path+genname+extension, 2, 0);
	cout << "Saída normal" << endl;
    return 0;
}

/*
 * topology.hpp
 *
 *  Created on: 13 de jun de 2017
 *      Author: noedelima
 */

#ifndef SRC_TOPOLOGY_HPP_
#define SRC_TOPOLOGY_HPP_
#include "DataStruct/database.hpp"
#include "network.hpp"
#include <Eigen/Eigen>
using namespace Eigen;

queue<branch> branches(string netname);
queue<generator> generators(string genname);
database<bus, int> buses(queue<branch> branches);
MatrixXcd Ybus(string branchpath);
MatrixXcd Zbus(string branchpath, string genpath);

#endif /* SRC_TOPOLOGY_HPP_ */

#ifndef _LMC_UTITLITY_INCLUDE_ENCODINGUTITLITY_H_
#define _LMC_UTITLITY_INCLUDE_ENCODINGUTITLITY_H_

#include <set>
#include <string>
#include <unordered_map>
#include <eigen3/Eigen/Dense>
#include "Element.hpp"
using namespace std;
using namespace Eigen;

unordered_map<string, RowVectorXd> GetOneHotEncodeHashmap(set<Element> elementSet);

#endif //_LMC_UTITLITY_UTITLITY_H_


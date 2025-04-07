#ifndef LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_
#define LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "Eigen/Dense"

template <typename T>
void print1DVector(const std::vector<T>& vec) {
    std::cout << "{";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;  // Removed the extra comma before the closing brace
}


void print2DVector(const std::vector<std::vector<size_t>> &vec);

void print2DStringVector(const std::vector<std::vector<std::string>> &vec);

void print3DVector(const std::vector<std::vector<std::vector<size_t>>> &vec);

void printOneHotEncodeHashmap(const std::unordered_map<std::string, Eigen::RowVectorXd> &hashmap);

#endif // LMC_CONSTANT_INCLUDE_PRINTUTILITY_H_

#include "PrintUtility.h"



void print2DVector(const std::vector<std::vector<size_t>> &vec) 
{
    std::cout << "{\n";
    for (const auto &row : vec) {
        std::cout << "  {";
        for (size_t i = 0; i < row.size(); ++i) {
            std::cout << row[i] << ", ";
        }
        std::cout << "},\n";
    }
    std::cout << "}\n";
};

void print2DStringVector(const std::vector<std::vector<std::string>> &vec) 
{
    std::cout << "{\n";
    for (const auto &row : vec) {
        std::cout << "  {";
        for (size_t i = 0; i < row.size(); ++i) {
            std::cout << row[i] << ", ";
        }
        std::cout << "},\n";
    }
    std::cout << "}\n";
};

void print3DVector(const std::vector<std::vector<std::vector<size_t>>>& vec) 
{
    std::cout << "{\n";
    for (const auto& matrix : vec) {  // Iterate over first dimension
        std::cout << "  {\n";
        for (const auto& row : matrix) {  // Iterate over second dimension
            std::cout << "    {";
            for (size_t i = 0; i < row.size(); ++i) {  // Iterate over third dimension
                std::cout << row[i] + 1;
                if (i < row.size() - 1) std::cout << ", ";
            }
            std::cout << "},\n";
        }
        std::cout << "  },\n";
    }
    std::cout << "}\n";
}

void printOneHotEncodeHashmap(const std::unordered_map<std::string, Eigen::RowVectorXd>& hashmap) 
{
    std::cout << "{\n";
    for (const auto& pair : hashmap) {
        std::cout << "  \"" << pair.first << "\": [";
        for (size_t i = 0; i < pair.second.size(); ++i) {
            std::cout << pair.second[i];
            if (i < pair.second.size() - 1) std::cout << ", ";
        }
        std::cout << "],\n";
    }
    std::cout << "}\n";
}

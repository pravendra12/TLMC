#include "PrintUtility.h"


void print2DStringVector(const std::vector<std::vector<std::string>> &vec)
{
    std::cout << "{\n";
    for (const auto &row : vec)
    {
        std::cout << "  {";
        for (size_t i = 0; i < row.size(); ++i)
        {
            std::cout << row[i] << ", ";
        }
        std::cout << "},\n";
    }
    std::cout << "}\n";
};

void print3DVector(const std::vector<std::vector<std::vector<size_t>>> &vec)
{
    std::cout << "{\n";
    for (const auto &matrix : vec)
    { // Iterate over first dimension
        std::cout << "  {\n";
        for (const auto &row : matrix)
        { // Iterate over second dimension
            std::cout << "    {";
            for (size_t i = 0; i < row.size(); ++i)
            { // Iterate over third dimension
                std::cout << row[i] + 1;
                if (i < row.size() - 1)
                    std::cout << ", ";
            }
            std::cout << "},\n";
        }
        std::cout << "  },\n";
    }
    std::cout << "}\n";
}

void printOneHotEncodeHashmap(const std::unordered_map<std::string, Eigen::RowVectorXd> &hashmap)
{
    std::cout << "{\n";
    for (const auto &pair : hashmap)
    {
        std::cout << "  \"" << pair.first << "\": [";
        for (size_t i = 0; i < pair.second.size(); ++i)
        {
            std::cout << pair.second[i];
            if (i < pair.second.size() - 1)
                std::cout << ", ";
        }
        std::cout << "],\n";
    }
    std::cout << "}\n";
}
/*
void writeToDataFile(const pair<size_t, size_t> &latticeIdJumpPair,
                     const unordered_map<Vector3d, size_t, Vector3dHash> &positionLatticeIdMap,
                     vector<vector<size_t>> &equivalentLatticeIdVector)
// const vector<pair<double, vector<size_t>>> &sortable_classes)
{
    ofstream outFile;
    string filename = to_string(latticeIdJumpPair.first) + "_" +
                      to_string(latticeIdJumpPair.second) + "_reducedPositionLCE.dump";

    outFile.open(filename);

    outFile << "ITEM: TIMESTEP\n0\n";
    outFile << "ITEM: NUMBER OF ATOMS\n"
            << positionLatticeIdMap.size() << "\n";
    outFile << "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n";
    outFile << "ITEM: ATOMS id type x y z group\n";

    // Create a lookup table for quick ID-to-class_value mapping
    unordered_map<size_t, int> idToClassValue;
    int i = 0;
    for (const auto &eqSites : equivalentLatticeIdVector)
    {
        int group = i;
        for (size_t site : eqSites)
        {
            idToClassValue[site] = group;
        }
        i += 1;
    }

    for (const auto &entry : positionLatticeIdMap)
    {
        size_t id = entry.second;
        outFile << id << " 1 ";
        for (const auto &coord : entry.first)
        {
            outFile << coord << " ";
        }

        // Add the class value or -1 if not found
        if (idToClassValue.count(id))
        {
            outFile << idToClassValue[id];
        }
        else
        {
            outFile << -1.0; // or some default/error value
        }

        outFile << "\n";
    }

    outFile.close();
}
    */
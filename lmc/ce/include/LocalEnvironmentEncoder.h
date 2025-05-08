#ifndef LMC_CE_INCLUDE_LOCALENVIRONMENTENCODER_H_
#define LMC_CE_INCLUDE_LOCALENVIRONMENTENCODER_H_

#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <unordered_map>
#include "Config.h"
#include "Eigen/Dense"
#include "CorrelationFunction.h"
#include "PrintUtility.h"

using namespace std;
using namespace Eigen;

/**
 * @brief Computes the local environment encoding for a given configuration.
 *
 * This function generates a row vector encoding that represents the local
 * environment of a configuration based on the provided parameters. The encoding
 * is determined using the specified basis type and orbit encoding map.
 *
 * @param config The configuration object containing the lattice structure and
 *               other relevant data.
 * @param elementSet A set of elements to consider for encoding.
 * @param basisType A string specifying the type of basis to use for encoding.
 *                  Example: "polynomial", "fourier", etc.
 * @param orbitEncodingMap A map that defines the encoding of orbits. The key is
 *                         a string representing the orbit type, and the value
 *                         is a vector of vectors, where each inner vector
 *                         contains indices representing a specific orbit.
 *                         Example format:
 *                         {
 *                             "1": {{1}, {14}},
 *                             "2": {{3}},
 *                             "12": {{1, 3}, {3, 14}}
 *                         }
 * @param symmetricallSortedVector A vector of size_t values representing a
 *                                  symmetrically sorted sequence of indices.
 *
 * @return A RowVectorXd object containing the computed local environment encoding.
 */
RowVectorXd GetLocalEnvironmentEncoding(
    const Config &config,
    const set<Element> &elementSet,
    const string &basisType,
    const map<string, vector<vector<size_t>>> orbitEncodingMap,
    const vector<size_t> symmetricallSortedVector);

#endif // LMC_CE_INCLUDE_LOCALENVIRONMENTENCODER_H_
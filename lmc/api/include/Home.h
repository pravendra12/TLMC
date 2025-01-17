/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/6/23 3:05 PM                                                            *
 **************************************************************************************************/

/*! \file  Home.h
 *  \brief Provides the core API for linking simulation parameters with Monte Carlo simulation methods
 */

#ifndef LMC_API_INCLUDE_HOME_H_
#define LMC_API_INCLUDE_HOME_H_

#include "Parameter.h"
#include "CanonicalMcSerial.h"
#include "KineticMcChainOmpi.h"
#include "KineticMcFirstMpi.h"

// #include "KineticMcFirstOmp.h"
// #include "SimulatedAnnealing.h"
// #include "CanonicalMcOmp.h"
// #include "Traverse.h"

namespace api {

  /*! \brief Prints the details of the given simulation parameters. 
   *  \param parameter : The `Parameter` object containing simulation settings to 
   *                     be printed.
   */
  void Print(const Parameter &parameter);
  
  /*! \brief Executes a Monte Carlo simulation using the provided parameters.
   *  \param parameter : The `Parameter` object containing simulation settings to 
   *                     be printed.
   */
  void Run(const Parameter &parameter);
  
  /*! \brief Constructs a serial implementation of a Canonical Monte Carlo 
   *         simulation.
   *  \param parameter : The `Parameter` object containing the configuration for 
   *                     the simulation.
   *  \return          : A `mc::CanonicalMcSerial` object configured with the 
   *                     provided parameters.
   */
  mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter &parameter);
  
  /*! \brief Constructs a first-order MPI implementation of a Kinetic Monte Carlo 
   *         simulation.
   *  \param parameter : The `Parameter` object containing the configuration for 
   *                     the simulation.
   *  \return          : A `mc::KineticMcFirstMpi` object configured with the 
   *                     provided parameters.
   */
  mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter);
  
  /*! \brief Constructs a Second-order OMP implementation of a Kinetic Monte Carlo 
   *         simulation.
   *  \param parameter : The `Parameter` object containing the configuration for 
   *                     the simulation.
   *  \return          : A `mc::KineticMcChainOmpi` object configured with the 
   *                     provided parameters.
   */
  mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter &parameter);
  
// mc::KineticMcFirstOmp BuildKineticMcFirstOmpFromParameter(const Parameter &parameter);
// ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter);
// mc::CanonicalMcOmp BuildCanonicalMcOmpFromParameter(const Parameter &parameter);
// ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter);
} // api

#endif //LMC_API_INCLUDE_HOME_H_

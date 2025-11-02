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
#include "Traverse.h"
#include "SubLatticeOccupancy.h"
#include "VacancyMigrationPredictorTLMC.h"
#include "EnergyPredictorTLMC.h"
#include "KRAPredictorTLMC.h"
#include "LVFEPredictorTLMC.h"
#include "ClusterExpansionParameters.h"
#include "KineticMcFirstMpi.h"
#include "CanonicalMcSerial.h"
#include "TiledSupercell.h"
#include "Cube.h"
#include "CanonicalMcOmp.h"
#include "ConvertAtomVectorsToConfigs.h"


using namespace std;

namespace api
{

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

  void RunCanonicalMcSerialFromParameter(const Parameter &parameter);

  void RunCanonicalMcOmpFromParameter(const Parameter &parameter);

  void RunKineticMcFirstMpiFromParameter(const Parameter &parameter);

  void ProfileEnergyPredictor(const Parameter &parameter);

  void RunAnsysFromParameter(const Parameter &parameter);


} // api

#endif // LMC_API_INCLUDE_HOME_H_

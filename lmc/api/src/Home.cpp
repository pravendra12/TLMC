/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                    *
 * @Last Modified time: 11/28/24 8:30 PM                                                            *
 **************************************************************************************************/

/*! \file  Home.h
 *  \brief File for the Home class implementation.
 */

#include "Home.h"

namespace api
{

  void Print(const Parameter &parameter)
  {
    cout << "Parameters" << endl;
    cout << "simulation_method: " << parameter.method << endl;

    if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcFirstOmp" ||
        parameter.method == "KineticMcChainOmp" || parameter.method == "KineticMcChainMpi" ||
        parameter.method == "KineticMcChainOmpi")
    {

      cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
           << endl;
      cout << "config_filename: " << parameter.config_filename_ << endl;
      cout << "map_filename: " << parameter.map_filename_ << endl;
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      cout << "lattice_param: " << parameter.lattice_param_ << endl;
      cout << "structure_type: " << parameter.structure_type_ << endl;
      cout << "log_dump_steps: " << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: " << parameter.config_dump_steps_ << endl;
      cout << "maximum_steps: " << parameter.maximum_steps_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
           << endl;
      cout << "temperature: " << parameter.temperature_ << endl;
      cout << "restart_steps: " << parameter.restart_steps_ << endl;
      cout << "restart_energy: " << parameter.restart_energy_ << endl;
      cout << "restart_time: " << parameter.restart_time_ << endl;
      cout << "time_temperature_filename: " << parameter.time_temperature_filename_
           << endl;
      cout << "rate_corrector: " << parameter.rate_corrector_ << endl;
      cout << "cube_size: " << parameter.cube_size_ << endl;
      cout << "vacancy_lattice_id: " << parameter.vacancayLatticeId_ << endl;
      cout << "small_config_id: " << parameter.smallConfigId_ << endl;
      cout << "vacancy_trajectory: " << parameter.vacancy_trajectory_
           << endl;
    }

    else if (parameter.method == "SimulatedAnnealing")
    {
      // General configuration information
      cout << "config_filename: "
           << parameter.config_filename_ << endl;
      cout << "json_coefficients_filename: "
           << parameter.json_coefficients_filename_ << endl;

      // Lattice and structure information
      cout << "lattice_param: "
           << parameter.lattice_param_ << endl;
      cout << "structure_type: "
           << parameter.structure_type_ << endl;

      // Simulation parameters
      cout << "maximum_steps: "
           << parameter.maximum_steps_ << endl;
      cout << "log_dump_steps: "
           << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: "
           << parameter.config_dump_steps_ << endl;
      cout << "restart_steps: "
           << parameter.restart_steps_ << endl;
      cout << "restart_energy: "
           << parameter.restart_energy_ << endl;

      // Cluster expansion settings
      cout << "max_cluster_size: "
           << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: "
           << parameter.max_bond_order_ << endl;
      // Supercell size used for CE model fitting
      cout << "supercell_size: "
           << parameter.supercell_size_ << endl;

      // Cutoffs
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;

      // Temperature settings
      cout << "initial_temperature: "
           << parameter.initial_temperature_ << endl;
    }

    else if (parameter.method == "CanonicalMcSerial" || parameter.method == "CanonicalMcOmp")
    {
      cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
           << endl;
      cout << "config_filename: " << parameter.config_filename_ << endl;
      cout << "map_filename: " << parameter.map_filename_ << endl;
      cout << "supercell_size: " << parameter.supercell_size_ << endl;
      cout << "lattice_param: " << parameter.lattice_param_ << endl;
      cout << "structure_type: " << parameter.structure_type_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "log_dump_steps: " << parameter.log_dump_steps_ << endl;
      cout << "config_dump_steps: " << parameter.config_dump_steps_ << endl;
      cout << "maximum_steps: " << parameter.maximum_steps_ << endl;
      cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
           << endl;
      cout << "temperature: " << parameter.temperature_ << endl;
      cout << "restart_steps: " << parameter.restart_steps_ << endl;
      cout << "restart_energy: " << parameter.restart_energy_ << endl;
      // cout << "initial_temperature: " << parameter.initial_temperature_ << endl;
      // cout << "decrement_temperature: " << parameter.decrement_temperature_ << endl;
    }

    else if (parameter.method == "Ansys" || parameter.method == "Reformat")
    {
      cout << "initial_steps: " << parameter.initial_steps_ << endl;
      cout << "increment_steps: " << parameter.increment_steps_ << endl;
      cout << "cutoffs: ";
      transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                ostream_iterator<string>(cout, " "),
                [](auto cutoff)
                { return to_string(cutoff); });
      cout << endl;
      cout << "log_type: " << parameter.log_type_ << endl;
      cout << "config_type: " << parameter.config_type_ << endl;

      cout << "Parameters for CE Encoding" << endl;
      cout << "extract_encoding: " << parameter.extract_encoding_ << endl;
      cout << "max_bond_order: " << parameter.max_bond_order_ << endl;
      cout << "max_cluster_size: " << parameter.max_cluster_size_ << endl;
    }
    else if (parameter.method == "ConvertAtomVectorsToConfigs")
    {
      cout << "path_tlmc_output: " << parameter.path_tlmc_output_ << endl;
    }
  }

  void Run(const Parameter &parameter)
  {

    if (parameter.method == "CanonicalMcSerial")
    {
      api::RunCanonicalMcSerialFromParameter(parameter);
    }
    else if (parameter.method == "KineticMcFirstMpi")
    {
      api::RunKineticMcFirstMpiFromParameter(parameter);
    }
    else if (parameter.method == "ProfileEnergyPredictor")
    {
      api::ProfileEnergyPredictor(parameter);
    }
    else if (parameter.method == "CanonicalMcOmp")
    {
      api::RunCanonicalMcOmpFromParameter(parameter);
    }
    else if (parameter.method == "ConvertAtomVectorsToConfigs")
    {
      ConvertAtomVectorsToConfigs(parameter.path_tlmc_output_);
    }
    else if (parameter.method == "ConvertCFGToXYZ")
    {
      if (parameter.input_filepath_.empty() || parameter.output_filepath_.empty())
      {
        throw std::runtime_error(
            "Error in ConvertCFGToXYZ: Both input and output file paths must be provided.");
      }

      // Optional: you can also check if the files exist
      if (!fs::exists(parameter.input_filepath_))
      {
        throw std::runtime_error(
            "Error in ConvertCFGToXYZ: Input file does not exist: " + parameter.input_filepath_);
      }

      auto cfg = Config::ReadCfg(parameter.input_filepath_);

      Config::WriteXyzExtended(
          parameter.output_filepath_,
          cfg,
          {},
          {});

      std::cout << "[INFO] Converted " << parameter.input_filepath_
                << " to " << parameter.output_filepath_ << std::endl;
    }
  }

  void RunCanonicalMcSerialFromParameter(const Parameter &parameter)
  {
    Config smallConfig = Config::GenerateSupercell(
        parameter.supercell_size_,
        parameter.lattice_param_,
        "X",
        parameter.structure_type_);

    // Declare the cube object
    Cube cubeObj(parameter.cube_size_);

    // Declare the tiled supercell
    TiledSupercell tiledSupercell(
        smallConfig,
        cubeObj);

    vector<uint64_t> atomIndexVector = TiledSupercell::ReadAtomicIndicesFromFile(
        parameter.atomic_indices_filename_);

    // Update the atom vector
    tiledSupercell.UpdateAtomVector(atomIndexVector);

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    smallConfig.UpdateNeighborList({maxClusterCutoff});

    // will be dependent on the size of cutoff list using which smallConfig was updated
    tiledSupercell.UpdateNeighbourLists(1);

    // Generate a small config for declaring symCE
    const size_t primSize = 2;

    auto primConfig = Config::GenerateSupercell(
        primSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictorTLMC symCEEnergyPredictor(
        ceParams,
        tiledSupercell,
        primConfig);

    // Again update the neighbor list
    smallConfig.UpdateNeighborList(parameter.cutoffs_);
    tiledSupercell.UpdateNeighbourLists(parameter.cutoffs_.size());

    // Declare LVFE Predictor
    LVFEPredictorTLMC lvfePredictor(
        ceParams,
        tiledSupercell);

    // Declare Energy Predictor
    EnergyPredictorTLMC energyChangePredictor(
        symCEEnergyPredictor,
        lvfePredictor);

    // Update to first nn
    tiledSupercell.UpdateNeighbourLists(1);

    cout << "Finish config reading. Start CMC." << endl;

    mc::CanonicalMcSerial canonicalMcSerial(
        tiledSupercell,
        parameter.log_dump_steps_,
        parameter.config_dump_steps_,
        parameter.maximum_steps_,
        parameter.thermodynamic_averaging_steps_,
        parameter.restart_steps_,
        parameter.restart_steps_,
        parameter.temperature_,
        energyChangePredictor);

    canonicalMcSerial.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void RunCanonicalMcOmpFromParameter(const Parameter &parameter)
  {
    Config smallConfig = Config::GenerateSupercell(
        parameter.supercell_size_,
        parameter.lattice_param_,
        "X",
        parameter.structure_type_);

    // Declare the cube object
    Cube cubeObj(parameter.cube_size_);

    // Declare the tiled supercell
    TiledSupercell tiledSupercell(
        smallConfig,
        cubeObj);

    if (!parameter.config_filename_.empty() && parameter.atomic_indices_filename_.empty())
    {
      cout << "[INFO] Using configuration file: " << parameter.config_filename_ << endl;

      Config largeConfig = Config::ReadCfg(parameter.config_filename_);
      tiledSupercell.UpdateAtomVector(largeConfig);
    }
    else if (parameter.config_filename_.empty() && !parameter.atomic_indices_filename_.empty())
    {
      cout << "[INFO] Using atomic indices file: " << parameter.atomic_indices_filename_ << endl;

      vector<uint64_t> atomIndexVector = TiledSupercell::ReadAtomicIndicesFromFile(
          parameter.atomic_indices_filename_);
      tiledSupercell.UpdateAtomVector(atomIndexVector);
    }
    else if (!parameter.config_filename_.empty() && !parameter.atomic_indices_filename_.empty())
    {
      cerr << "[ERROR] Both configuration and atomic indices files are provided. "
           << "Please specify only one input source." << endl;
      return;
    }
    else
    {
      cerr << "[ERROR] Neither configuration file nor atomic indices file provided in parameters."
           << endl;
      return;
    }

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    smallConfig.UpdateNeighborList({maxClusterCutoff});

    // will be dependent on the size of cutoff list using which smallConfig was updated
    tiledSupercell.UpdateNeighbourLists(1);

    // Generate a small config for declaring symCE
    const size_t primSize = 2;

    auto primConfig = Config::GenerateSupercell(
        primSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictorTLMC symCEEnergyPredictor(
        ceParams,
        tiledSupercell,
        primConfig);

    // Again update the neighbor list
    smallConfig.UpdateNeighborList(parameter.cutoffs_);
    tiledSupercell.UpdateNeighbourLists(parameter.cutoffs_.size());

    // Declare LVFE Predictor
    LVFEPredictorTLMC lvfePredictor(
        ceParams,
        tiledSupercell);

    // Declare Energy Predictor
    EnergyPredictorTLMC energyChangePredictor(
        symCEEnergyPredictor,
        lvfePredictor);

    // Update to first nn
    tiledSupercell.UpdateNeighbourLists(1);

    cout << "Finish config reading. Start CMC." << endl;

    mc::CanonicalMcOmp canonicalMcOmp(
        tiledSupercell,
        parameter.log_dump_steps_,
        parameter.config_dump_steps_,
        parameter.maximum_steps_,
        parameter.thermodynamic_averaging_steps_,
        parameter.restart_steps_,
        parameter.restart_steps_,
        parameter.temperature_,
        energyChangePredictor);

    canonicalMcOmp.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void RunKineticMcFirstMpiFromParameter(const Parameter
                                             &parameter)
  {
    Config smallConfig = Config::GenerateSupercell(
        parameter.supercell_size_,
        parameter.lattice_param_,
        "X",
        parameter.structure_type_);

    // Declare the cube object
    Cube cubeObj(parameter.cube_size_);

    // Declare the tiled supercell
    TiledSupercell tiledSupercell(
        smallConfig,
        cubeObj);

    vector<uint64_t> atomIndexVector = TiledSupercell::ReadAtomicIndicesFromFile(
        parameter.atomic_indices_filename_);

    // Update the atom vector
    tiledSupercell.UpdateAtomVector(atomIndexVector);

    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    smallConfig.UpdateNeighborList({maxClusterCutoff});

    // will be dependent on the size of cutoff list using which smallConfig was updated
    tiledSupercell.UpdateNeighbourLists(1);

    // Generate a small config for declaring symCE
    const size_t primSize = 2;

    auto primConfig = Config::GenerateSupercell(
        primSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictorTLMC symCEEnergyPredictor(
        ceParams,
        tiledSupercell,
        primConfig);

    // Again update the neighbor list
    smallConfig.UpdateNeighborList(parameter.cutoffs_);
    tiledSupercell.UpdateNeighbourLists(parameter.cutoffs_.size());

    // Declare LVFE Predictor
    LVFEPredictorTLMC lvfePredictor(
        ceParams,
        tiledSupercell);

    // Declare Energy Predictor
    EnergyPredictorTLMC energyChangePredictor(
        symCEEnergyPredictor,
        lvfePredictor);

    // Declare KRA Predictor
    KRAPredictorTLMC eKRAPredictor(
        ceParams,
        tiledSupercell);

    // Declare vacancy migration predictor
    VacancyMigrationPredictorTLMC vacancyMigrationPredictor(
        eKRAPredictor,
        energyChangePredictor);

    // After this config can be again update as it will not be used again

    // Update to first nn
    tiledSupercell.UpdateNeighbourLists(1);

    cout << "Finish config reading. Start KMC." << endl;

    mc::KineticMcFirstMpi kmcFirstMpi(tiledSupercell,
                                      parameter.log_dump_steps_,
                                      parameter.config_dump_steps_,
                                      parameter.maximum_steps_,
                                      parameter.thermodynamic_averaging_steps_,
                                      parameter.restart_steps_,
                                      parameter.restart_energy_,
                                      parameter.restart_time_,
                                      parameter.temperature_,
                                      vacancyMigrationPredictor,
                                      parameter.time_temperature_filename_,
                                      parameter.rate_corrector_,
                                      parameter.vacancy_trajectory_);

    kmcFirstMpi.Simulate();

    cout << "Simulation Completed" << endl;
  }

  void ProfileEnergyPredictor(const Parameter &parameter)
  {
    // Read CE Parameters
    ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

    auto allowedElements = ceParams.GetAllowedElements();

    vector<double> compositionVector;
    compositionVector.reserve(allowedElements.size());

    int n = allowedElements.size();
    int base = 100 / n;
    int remainder = 100 % n;

    // Initialize all parts to 'base'
    compositionVector.assign(n, base);

    // Add 1 to the first 'remainder' parts
    for (int i = 0; i < remainder; ++i)
    {
      compositionVector[i] += 1.0;
    }

    Config smallConfig = Config::GenerateAlloySupercell(
        parameter.supercell_size_,
        parameter.lattice_param_,
        parameter.structure_type_,
        allowedElements,
        compositionVector,
        1);

    // Declare the cube object
    Cube cubeObj(parameter.cube_size_);

    // Declare the tiled supercell
    TiledSupercell tiledSupercell(
        smallConfig,
        cubeObj);

    size_t largeSupercellSize = parameter.supercell_size_ * parameter.cube_size_;

    Config largeConfig = Config::GenerateAlloySupercell(
        largeSupercellSize,
        parameter.lattice_param_,
        parameter.structure_type_,
        allowedElements,
        compositionVector,
        1);

    // Update the atom vector
    tiledSupercell.UpdateAtomVector(largeConfig);

    // Used to update the symmetric CE
    double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
    smallConfig.UpdateNeighborList({maxClusterCutoff});

    // will be dependent on the size of cutoff list using which smallConfig was updated
    tiledSupercell.UpdateNeighbourLists(1);

    // Generate a small config for declaring symCE
    const size_t primSize = 2;

    auto primConfig = Config::GenerateSupercell(
        primSize,
        parameter.lattice_param_,
        "Mo", // Does not matter as the lattice param and structure type is important
        parameter.structure_type_);

    // Declare Symmetric CE
    SymmetricCEPredictorTLMC symCEEnergyPredictor(
        ceParams,
        tiledSupercell,
        primConfig);

    // Again update the neighbor list
    smallConfig.UpdateNeighborList(parameter.cutoffs_);
    tiledSupercell.UpdateNeighbourLists(parameter.cutoffs_.size());

    // Declare LVFE Predictor
    LVFEPredictorTLMC lvfePredictor(
        ceParams,
        tiledSupercell);

    // Declare Energy Predictor
    EnergyPredictorTLMC energyChangePredictor(
        symCEEnergyPredictor,
        lvfePredictor);

    // Select a lattice pair with different elements
    // The central Id is choosen so no worry of periodic boundary
    // All first nn are in the same cell
    size_t latticeId1 = smallConfig.GetCentralAtomLatticeId();
    size_t smallConfigIdx = 1;

    auto latticeSite1 = LatticeSiteMapping(latticeId1, 1);
    auto element1 = tiledSupercell.GetElementAtSite(latticeSite1);

    LatticeSiteMapping latticeSite2;
    for (auto nnId : smallConfig.GetNeighborLatticeIdVectorOfLattice(latticeId1, 1))
    {
      auto nnSite = LatticeSiteMapping(nnId, smallConfigIdx);
      auto element2 = tiledSupercell.GetElementAtSite(nnSite);

      if (element1 != element2)
      {
        latticeSite2 = nnSite;
        break;
      }
    }

    // Profiling
    int maxThreads = omp_get_max_threads();
    cout << "Max Number of Threads: " << maxThreads << endl;

    vector<int> threadList = {1, 2, 4, 8, 12, 16, 32, 48};

    // --- Without vacancy ---
    cout << "=== Profiling without vacancy ===" << endl;
    for (const auto threadCount : threadList)
    {
      if (threadCount <= maxThreads)
      {
        cout << "[Without Vacancy] Threads: " << threadCount << " | ";
        energyChangePredictor.ProfileEnergyChange(
            tiledSupercell,
            make_pair(latticeSite1, latticeSite2),
            threadCount);
      }
    }

    // --- With vacancy ---
    tiledSupercell.SetElementAtSite(latticeSite1, Element("X"));
    cout << "\n=== Profiling with vacancy ===" << endl;
    for (const auto threadCount : threadList)
    {
      if (threadCount <= maxThreads)
      {
        cout << "[With Vacancy] Threads: " << threadCount << " | ";
        energyChangePredictor.ProfileEnergyChange(
            tiledSupercell,
            make_pair(latticeSite1, latticeSite2),
            threadCount);
      }
    }
  }

}

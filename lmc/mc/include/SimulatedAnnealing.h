#ifndef LMC_MC_INCLUDE_SIMULATEDANNEALING_H_
#define LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

#include <iostream>
#include <random>
#include "Config.h"
#include "McAbstract.h"
#include "PotentialEnergyEstimator.h"
#include "EnergyPredictor.h"
#include "ClusterExpansionParameters.h"

using namespace std;

namespace mc
{
  /**
   * @class SimulatedAnnealing
   * @brief This class implements the simulated annealing algorithm for Monte Carlo simulations.
   *
   * The class handles the setup and execution of a simulated annealing process, which includes
   * temperature updates, lattice jumps, and event selection using the Metropolis-Hastings algorithm.
   * It also integrates a cluster expansion-based potential energy estimator to predict energy changes
   * during the simulation. The class provides functions to manage simulation parameters,
   * temperature cooling, and event selection for atomic configurations.
   */
  class SimulatedAnnealing : public McAbstract
  {
  public:
    /**
     * @brief Constructs a SimulatedAnnealing object with the given configuration.
     *
     * @param config The configuration object for the simulation.
     * @param supercell_config The configuration of the supercell.
     * @param log_dump_steps The number of steps between logging simulation data.
     * @param config_dump_steps The number of steps between configuration dumps.
     * @param maximum_steps The maximum number of steps for the simulation.
     * @param restart_steps The number of steps before restarting the simulation.
     * @param restart_energy The energy threshold for simulation restart.
     * @param initial_temperature The initial temperature for the simulated annealing process.
     * @param element_set Set of elements in the system.
     * @param max_cluster_size Maximum size of the clusters in the simulation.
     * @param max_bond_order Maximum bond order for clusters.
     * @param json_coefficients_filename File name for cluster expansion coefficients in JSON forminitial_temperature_at.
     */
    SimulatedAnnealing(Config config,
                       const unsigned long long int logDumpSteps,
                       const unsigned long long int configDumpSteps,
                       const unsigned long long int maximumSteps,
                       const unsigned long long int restartSteps,
                       const double restartEnergy,
                       const double initialTemperature,
                       EnergyPredictor &energyChangePredictor);

    /**
     * @brief Starts the simulated annealing simulation.
     *
     * This function runs the main loop of the simulated annealing process, where it
     * updates the temperature, generates lattice jumps, selects events, and checks for
     * restart conditions.
     */
    void Simulate();

  private:
    /**
     * @brief Dumps the current configuration to an output.
     *
     * This function saves the current atomic configuration to an output file,
     * providing a way to track the system's evolution during the simulation.
     */
    void Dump() const;

    /**
     * @brief Updates the temperature using the exponential cooling function.
     *
     * The temperature is updated according to the exponential cooling scheme,
     * which is commonly used in simulated annealing algorithms to gradually lower
     * the temperature during the simulation.
     */
    void UpdateTemperature();

    /**
     * @brief Generates a random pair of lattice indices representing a possible jump.
     *
     * This function generates a pair of lattice site indices that represent a
     * potential jump in the Monte Carlo simulation.
     *
     * @return A pair of size_t values representing the lattice indices for a jump.
     */
    pair<size_t, size_t> GenerateLatticeIdJumpPair();

    /**
     * @brief Selects an event using the Metropolis-Hastings algorithm.
     *
     * This function decides whether to accept or reject an event based on the
     * Metropolis-Hastings criterion, which considers the energy change (dE)
     * associated with the proposed event and the current temperature.
     *
     * @param lattice_id_jump_pair A pair of lattice indices representing the jump.
     * @param dE The energy change associated with the event.
     */
    void SelectEvent(const pair<size_t, size_t> &lattice_id_jump_pair, double dE);

    /**
     * @brief An object that estimates the potential energy change using cluster expansion.
     *
     * This object is used to calculate the energy change associated with a
     * proposed event during the Monte Carlo simulation, using a cluster expansion model.
     */
    EnergyPredictor &energyChangePredictor_;

    /**
     * @brief A uniform integer distribution used to select atom indices randomly.
     *
     * This distribution is used to randomly select atoms during the simulation.
     */
    mutable uniform_int_distribution<size_t> atom_index_selector_;

    /**
     * @brief The initial temperature for the simulated annealing process.
     *
     * This value determines the starting temperature, which controls the likelihood
     * of accepting higher-energy configurations early in the simulation.
     */
    const double initialTemperature_;
  };
}

#endif // LMC_MC_INCLUDE_SIMULATEDANNEALING_H_

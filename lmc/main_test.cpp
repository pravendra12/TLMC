#include "Home.h"
#include "Config.h"
#include "ClusterExpansion.h"
#include "SymmetrySpglib.h"
#include "GetEquivalentClusters.h"
#include "SymmetrySpglib.h"
#include "TrainingUtility.h"
#include "ClusterExpansionParameters.h"
#include "LVFEPredictor.h"
#include "KRAPredictor.h"

using namespace std;
using namespace Eigen;
namespace fs = filesystem;

int main()
{
  ClusterExpansionParameters ceParams(
      "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/script/coefficientFile_MoTa_V3.1.json",
      true);

#include <vector>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>

  // Seed random generator once
  std::srand(std::time(nullptr));

  // Generate supercell
  auto cfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");
  cfg.UpdateNeighborList({13});

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");

  EnergyPredictor energyPredictor(
      ceParams,
      cfg,
      primCfg);

  exit(1);

  // Randomly assign Mo or Ta to each atom
  size_t numAtoms = cfg.GetNumAtoms();
  std::vector<std::string> elements = {"Mo", "Ta"};

  for (size_t i = 0; i < numAtoms; ++i)
  {
    // Pick randomly between "Mo" and "Ta"
    std::string elem = elements[std::rand() % elements.size()];
    cfg.SetElementOfAtom(i, Element(elem));
  }

  cfg.SetElementOfLattice(0, Element("Ta"));
  cfg.SetElementOfLattice(100, Element("Ta"));
  cfg.SetElementOfLattice(111, Element("Ta"));

  cfg.SetElementOfLattice(cfg.GetCentralAtomLatticeId(), Element("X"));

  LVFEPredictor lvfePredictor(ceParams,
                              cfg);

  size_t vacancyId = cfg.GetVacancyLatticeId();
  auto firstNN = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1); // 1NN

  for (auto neighborId : firstNN)
  {
    cfg.SetElementOfLattice(neighborId, Element("Ta")); // change element to Ta
  }

  size_t neighborId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];
  auto jumpPair = std::make_pair(vacancyId, neighborId);

  // LVFE: Effective vacancy formation energy
  double evfeBefore = lvfePredictor.GetEffectiveVacancyFormationEnergy(cfg, vacancyId);
  cout << "Effective Vacancy Formation Energy (initial): " << evfeBefore << endl;

  // LVFE: Migration energy (forward)
  double deForwardBefore = lvfePredictor.GetDeForVacancyMigration(cfg, jumpPair);
  cout << "Migration Energy (forward, before jump): " << deForwardBefore << endl;

  // Perform the jump
  cfg.LatticeJump(jumpPair);

  // LVFE: Migration energy (backward)
  double deBackwardAfter = lvfePredictor.GetDeForVacancyMigration(cfg, std::make_pair(neighborId, vacancyId));
  cout << "Migration Energy (backward, after jump): " << deBackwardAfter << endl;

  // LVFE: Effective vacancy formation energy at the new site
  double evfeAfter = lvfePredictor.GetEffectiveVacancyFormationEnergy(cfg, neighborId);
  cout << "Effective Vacancy Formation Energy (after jump): " << evfeAfter << endl;

  // KRA predictor
  KRAPredictor kraPredictor(ceParams, cfg);

  // Forward barrier using KRA
  double kRAForward = kraPredictor.GetKRA(cfg, jumpPair);
  cout << "KRA Forward Barrier: " << kRAForward << endl;

  // Backward barrier using KRA
  double kRABackward = kraPredictor.GetKRA(cfg, std::make_pair(neighborId, vacancyId));
  cout << "KRA Backward Barrier: " << kRABackward << endl;

  VacancyMigrationPredictor vacMigrationPredictor(
      ceParams,
      cfg);

  auto barrierDe = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);

  cout << barrierDe.first << endl;
  cout << barrierDe.second << endl;

  cfg.LatticeJump(jumpPair);

  barrierDe = vacMigrationPredictor.GetBarrierAndDeltaE(cfg, jumpPair);

  cout << barrierDe.first << endl;
  cout << barrierDe.second << endl;

  //  cfg.UpdateNeighborList({13.1});
  //
  //  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");
  //
  //  EnergyPredictor energyPredictor(
  //      ceParams,
  //      cfg,
  //      primCfg);
  //
  // IterateDirectoryToExtractData();
}

/*
int main()
{

  auto primCfg = Config::GenerateSupercell(2, 3.2, "Mo", "BCC");
  vector<string> allowedElements = {"Mo", "Ta"};
  vector<double> clusterCutoffs = {13, 8, 6};

  auto supercellCfg = Config::GenerateSupercell(10, 3.2, "Mo", "BCC");

  supercellCfg.UpdateNeighborList({3});
  auto centralId = supercellCfg.GetCentralAtomLatticeId();

  auto nnId = supercellCfg.GetNeighborLatticeIdVectorOfLattice(centralId, 1)[0];

  // Number of Ta atoms you want to add
  int numTa = 1000; // example

  // Get total number of sites
  int numSites = supercellCfg.GetNumAtoms();

  // Random number generator
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> dis(0, numSites - 1);

  // Keep track of which sites are already replaced
  vector<size_t> replacedSites;

  // Randomly replace Mo with Ta
  for (int i = 0; i < numTa; ++i)
  {
    size_t site;
    do
    {
      site = dis(gen);
    } while (find(replacedSites.begin(), replacedSites.end(), site) != replacedSites.end());

    replacedSites.push_back(site);
    supercellCfg.SetElementOfLattice(site, Element("Ta")); // assuming such a function exists
  }

  vector<double> maxCutoffVector = {
      *max_element(clusterCutoffs.begin(), clusterCutoffs.end())};
  supercellCfg.UpdateNeighborList(maxCutoffVector);

  string predictorFile = "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/icetECIs.json";
  EnergyPredictor energyPredictor(
      predictorFile,
      supercellCfg,
      primCfg,
      allowedElements,
      clusterCutoffs);

  auto energyBefore = energyPredictor.GetEnergyOfSiteWithVacancy(supercellCfg,
                                                                 centralId);

  auto elvfeBefore = energyPredictor.GetEffectiveVacancyFormationEnergy(
      supercellCfg,
      centralId);

  cout << "After JUMP" << endl;
  supercellCfg.SetElementOfLattice(centralId, supercellCfg.GetElementOfLattice(nnId));

  auto energyAfter = energyPredictor.GetEnergyOfSiteWithVacancy(supercellCfg,
                                                                nnId);

  auto elvfeAfter = energyPredictor.GetEffectiveVacancyFormationEnergy(supercellCfg, nnId);


  cout << "Energy Before : " << energyBefore << endl;

  cout << "Energy After  : " << energyAfter << endl;

  cout << "ELVFE Before  : " << elvfeBefore << endl;

  cout << "ELVFE After   : " << elvfeAfter << endl;
  cout << "Change in Energy : " << energyAfter - energyBefore << endl;
  cout << "Change in ELVFE : " << elvfeAfter - elvfeBefore << endl;

  return 0;
}
  */
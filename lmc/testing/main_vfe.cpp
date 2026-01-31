#include "Home.h"
#include "ComputeVFE.h"
#include <filesystem>
#include <vector>
#include <string>
#include "SymmetricCE.h"
namespace fs = std::filesystem;

int main()
{
  ClusterExpansionParameters ceParams("//media/sf_Phd/MoTa/dft/neb/coefficients/MoTaX_fitted_coefficients_V3.8.json");

  auto cfgBase = Config::GenerateAlloySupercell(10, 3.2, "BCC", {"Mo", "Ta"}, {50, 50}, 0);

  double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
  cfgBase.UpdateNeighborList({maxClusterCutoff});

  const size_t primSize = 2;

  auto primConfig = Config::GenerateSupercell(
      primSize,
      3.2,
      "Mo",
      "BCC");

  // Declare Symmetric CE
  SymmetricCEPredictor symCEPredictor(
      ceParams,
      cfgBase,
      primConfig);

  std::ofstream fout("//media/sf_Phd/MoTa/vacancy_concentration/data/raw/vfe_output_ceMoTaX_v3.8.txt"); // overwrites; use ios::app to append
  if (!fout)
  {
    std::cerr << "ERROR: could not open vfe_output_ceMoTaX_v3.8.txt for writing\n";
    return 1;
  }

  fout << "seed,vfe_Mo_eV,vfe_Ta_eV\n";
  fout << std::setprecision(12);

  for (int i = 0; i < 1000; i++)
  {
    auto cfg = Config::GenerateAlloySupercell(10, 3.2, "BCC", {"Mo", "Ta"}, {50, 50}, i);

    auto vfeMap = ComputeVFE(cfg, ceParams.GetChemicalPotentialsMap(), symCEPredictor);

    // default to NaN if missing
    double vfeMo = std::numeric_limits<double>::quiet_NaN();
    double vfeTa = std::numeric_limits<double>::quiet_NaN();

    auto itMo = vfeMap.find(Element("Mo"));
    if (itMo != vfeMap.end())
      vfeMo = itMo->second;

    auto itTa = vfeMap.find(Element("Ta"));
    if (itTa != vfeMap.end())
      vfeTa = itTa->second;

    fout << i << "," << vfeMo << "," << vfeTa << "\n";

    // break; // remove if you want all 200
  }

  fout.close();

  return 0;
}

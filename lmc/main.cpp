#include "Home.h"
#include "KRAPredictor.h"

int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);

  // Read CE Parameters
  ClusterExpansionParameters ceParams(parameter.json_coefficients_filename_);

  Config smallConfig = Config::GenerateSupercell(
      parameter.supercell_size_,
      parameter.lattice_param_,
      "X",
      parameter.structure_type_);

  // Again update the neighbor list
  smallConfig.UpdateNeighborList(parameter.cutoffs_);

  // Declare KRA Predictor
  KRAPredictor eKRAPredictor(
      ceParams,
      smallConfig);

  // api::Run(parameter);
}

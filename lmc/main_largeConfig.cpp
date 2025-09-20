#include "TiledSupercell.h"
#include "Cube.h"

int main()
{
  string pathToAtomInfo = "/home/pravendra3/Documents/TLMC/bin/testingTLMC";

  auto smallCfg = Config::ReadConfig("/home/pravendra3/Documents/TLMC/bin/testingTLMC/Mo50Ta50_10x10x10.cfg.gz");
  Cube cubeObj(2);

  for (int i = 0; i <= 100; i++)
  {
    auto atomVector = TiledSupercell::ReadAtomVectorInfoFromFile(pathToAtomInfo + "/" + to_string(i) + ".txt");

    auto largeCfg = TiledSupercell::MakeSupercellFromAtomInfo(
        smallCfg,
        cubeObj,
        atomVector);

    Config::WriteConfig(
        pathToAtomInfo + "/configs/" + to_string(i) + ".cfg.gz",
        largeCfg);
  }
}

//
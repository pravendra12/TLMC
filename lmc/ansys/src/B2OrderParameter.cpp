#include "B2OrderParameter.h"

B2OrderParameter::B2OrderParameter(const Config &config) : config(config)
{
  InitializeAlphaLatticeSites();
  InitializeBetaLatticeSites();
}

double B2OrderParameter::GetB2OrderParameter(const Element &element)
{
  // fractional occupancy of element at alpha and beta sites
  auto alphaOccupancy = GetAlphaSiteOccupancy(element);
  auto betaOccupancy = GetBetaSiteOccupancy(element);

  double b2OrderParameter = (alphaOccupancy - betaOccupancy) / (alphaOccupancy + betaOccupancy);

  return b2OrderParameter;
}

double B2OrderParameter::GetAlphaSiteOccupancy(const Element &element)
{
  // Number of occupied alpha sites by element
  size_t numOccupiedASites = 0;
  size_t numAlphaSites = alphaLatticeSites.size();

  for (const auto &id : alphaLatticeSites)
  {
    auto siteElement = config.GetElementOfLattice(id);

    if (siteElement == element)
    {
      numOccupiedASites += 1;
    }
  }

  double alphaSiteOccupancy = double(numOccupiedASites) / double(numAlphaSites);

  return alphaSiteOccupancy;
}

double B2OrderParameter::GetBetaSiteOccupancy(const Element &element)
{
  // Number of occupied beta sites by element
  size_t numOccupiedBSites = 0;
  size_t numBetaSites = betaLatticeSites.size();

  for (const auto &id : betaLatticeSites)
  {
    auto siteElement = config.GetElementOfLattice(id);

    if (siteElement == element)
    {
      numOccupiedBSites += 1;
    }
  }

  double betaSiteOccupancy = double(numOccupiedBSites) / double(numBetaSites);

  return betaSiteOccupancy;
}

unordered_set<size_t> B2OrderParameter::GetAlphaLatticeSites()
{
  return alphaLatticeSites;
}

unordered_set<size_t> B2OrderParameter::GetBetaLatticeSites()
{
  return betaLatticeSites;
}

bool B2OrderParameter::isB2Ordered(const Config &config, const size_t atomId)
{
  // Element vacancy("X");
  
  // if (config.GetElementOfAtom(atomId) == vacancy)
  // {
  //   return false;
  // }

  Element centerElement = config.GetElementOfAtom(atomId);
  auto firstNN = config.GetNeighborAtomIdVectorOfAtom(atomId, 1);

  Element neighborElement = config.GetElementOfAtom(firstNN[0]);

  for (size_t i = 0; i < firstNN.size(); ++i)
  {
    if (config.GetElementOfAtom(firstNN[i]) != neighborElement)
    {
      return false;
    }
  }

  if (neighborElement == centerElement)
    return false;

  // second NN should all be equal to centerElement for perfect B2

  auto secondNN = config.GetNeighborAtomIdVectorOfAtom(atomId, 2);

  // case of defected B2

  for (size_t i = 0; i < secondNN.size(); ++i)
  {
    if (config.GetElementOfAtom(secondNN[i]) != centerElement)
    {
      return false;
    }
  }

  return true;
}

void B2OrderParameter::InitializeAlphaLatticeSites()
{

  auto secondNNList = config.GetNeighborLists()[1];

  for (size_t id1 = 0; id1 < secondNNList.size(); id1++)
  {
    const auto &secondNN = secondNNList[id1];

    for (size_t id2 : secondNN)
    {
      bool id1Valid = true;
      bool id2Valid = true;

      // Check if id1 has any bond order of 1 with sites already in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (config.GetDistanceOrder(id1, id3) == 1)
        {
          id1Valid = false;
          break;
        }
      }

      // Check if id2 has any bond order of 1 with sites already in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (config.GetDistanceOrder(id2, id3) == 1)
        {
          id2Valid = false;
          break;
        }
      }

      // Add valid sites to alphaLatticeSites
      if (id1Valid)
      {
        alphaLatticeSites.emplace(id1);
      }
      if (id2Valid)
      {
        alphaLatticeSites.emplace(id2);
      }
    }
  }
}

void B2OrderParameter::InitializeBetaLatticeSites()
{
  // betaSites = allSites - alphaSites
  for (size_t id = 0; id < numLattice; id++)
  {
    // O(1) Time Complexity
    if (alphaLatticeSites.find(id) == alphaLatticeSites.end())
    {
      betaLatticeSites.emplace(id);
    }
  }
}

#include "EncodingUtility.h"

unordered_map<string, RowVectorXd> GetOneHotEncodeHashmap(set<Element> elementSet) 
{
  if (elementSet.contains(Element("X")))
  {
    elementSet.erase(Element("X"));
  }

  size_t numElements = elementSet.size();
  unordered_map<string, RowVectorXd> encodeDict;
  
  size_t ct1 = 0;
  for (const auto &element : elementSet) 
  {
    RowVectorXd elementEncode = RowVectorXd::Zero(numElements);
    elementEncode(ct1) = 1.0;
    encodeDict[element.GetElementString()] = elementEncode;
    ++ct1;
  }
  
  size_t numPairs = numElements * numElements;
  size_t ct2 = 0;
  for (const auto &element1 : elementSet) 
  {
    for (const auto &element2 : elementSet) 
    {
      RowVectorXd elementEncode = RowVectorXd::Zero(numPairs);
      elementEncode(ct2) = 1.0;
      encodeDict[element1.GetElementString() + element2.GetElementString()] = elementEncode;
      ++ct2;
    }
  }
  
  size_t numPairsSymmetry = (numElements + 1) * numElements / 2;
  size_t ct3 = 0;
  for (auto it1 = elementSet.cbegin(); it1 != elementSet.cend(); ++it1) 
  {
    for (auto it2 = it1; it2 != elementSet.cend(); ++it2) 
    {
      RowVectorXd elementEncode = RowVectorXd::Zero(numPairsSymmetry);
      elementEncode(ct3) = 1.0;
      encodeDict[it1->GetElementString() + '-' + it2->GetElementString()] = elementEncode;
      ++ct3;
    }
  }
  
  return encodeDict;
}


#include "ClusterExpansionParameters.h"

ClusterExpansionParameters::ClusterExpansionParameters(
    const string &coefficientFilename,
    const bool debug) : predictorFilename_(coefficientFilename)
{

  // Open the JSON file
  ifstream ifs(coefficientFilename);
  if (!ifs.is_open())
  {
    cerr << "Error in `ClusterExpansionParameters`: Unable to open file "
         << coefficientFilename << endl;
    throw runtime_error("Unable to open JSON file: " + coefficientFilename);
  }

  ifs >> allParameters_;

  if (debug)
  {
    DebugAllFunctions();
  }
}

string ClusterExpansionParameters::GetCoefficientFile() const
{
  return predictorFilename_;
}

string ClusterExpansionParameters::GetBasisType() const
{
  string jsonKey = "basisType";

  string basisType = ReadParameterFromJson<string>(allParameters_, jsonKey);

  return basisType;
}

size_t ClusterExpansionParameters::GetMaxBondOrderCE() const
{
  string jsonKey = "maxBondOrderCE";

  size_t maxBondOrderCE = ReadParameterFromJson<size_t>(allParameters_, jsonKey);

  return maxBondOrderCE;
}

size_t ClusterExpansionParameters::GetMaxClusterSizeCE() const
{
  string jsonKey = "maxClusterSizeCE";

  size_t maxClusterSizeCE = ReadParameterFromJson<size_t>(allParameters_, jsonKey);

  return maxClusterSizeCE;
}

set<Element> ClusterExpansionParameters::GetElementSetCE() const
{
  string jsonKey = "ce";
  string jsonSubKey = "elementSet";

  vector<string> elementVector = ReadParametersFromJson<string>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  set<Element> elementSet;
  for (const auto &ele : elementVector)
  {
    elementSet.insert(Element(ele));
  }

  return elementSet;
}

VectorXd ClusterExpansionParameters::GetECIs() const
{
  string jsonKey = "ce";
  string jsonSubKey = "eci";

  vector<double> eciVector = ReadParametersFromJson<double>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  VectorXd eciVectorXd = Map<VectorXd>(eciVector.data(), eciVector.size());

  return eciVectorXd;
}

size_t ClusterExpansionParameters::GetMaxBondOrderOfClusterKRA() const
{
  string jsonKey = "maxBondOrderOfClusterKRA";

  size_t maxBondOrderOfClusterKRA = ReadParameterFromJson<size_t>(allParameters_, jsonKey);

  return maxBondOrderOfClusterKRA;
}

size_t ClusterExpansionParameters::GetMaxBondOrderKRA() const
{
  string jsonKey = "maxBondOrderKRA";

  size_t maxBondOrderKRA = ReadParameterFromJson<size_t>(allParameters_, jsonKey);

  return maxBondOrderKRA;
}

size_t ClusterExpansionParameters::GetMaxClusterSizeKRA() const
{
  string jsonKey = "maxClusterSizeKRA";

  size_t maxClusterSizeKRA = ReadParameterFromJson<size_t>(allParameters_, jsonKey);

  return maxClusterSizeKRA;
}

set<Element> ClusterExpansionParameters::GetElementSetKRA() const
{
  string jsonKey = "kra";
  string jsonSubKey = "elementSet";

  vector<string> elementVector = ReadParametersFromJson<string>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  set<Element> elementSet;
  for (const auto &ele : elementVector)
  {
    elementSet.insert(Element(ele));
  }

  return elementSet;
}

/*
unordered_map<Element, VectorXd, boost::hash<Element>> ClusterExpansionParameters::GetKECIs() const
{
  unordered_map<Element, VectorXd, boost::hash<Element>> keciMap;

  json keciJson = ReadParameterFromJson<json>(allParameters_["kra"], "keci");

  // Iterate over each element in the JSON object
  for (auto &[key, value] : keciJson.items())
  {
    vector<double> vec = value.get<vector<double>>();
    keciMap[Element{key}] = Map<VectorXd>(vec.data(), vec.size());
  }

  return keciMap;
}
  */

VectorXd ClusterExpansionParameters::GetKECIs() const
{
  string jsonKey = "kra";
  string jsonSubKey = "keci";

  vector<double> eciVector = ReadParametersFromJson<double>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  VectorXd eciVectorXd = Map<VectorXd>(eciVector.data(), eciVector.size());

  return eciVectorXd;
}

void ClusterExpansionParameters::DebugAllFunctions()
{
  cout << "------ Testing Cluster Expansion Parameters Reading ------\n\n";

  cout << "BasisType: " << GetBasisType() << "\n\n";

  cout << "Cluster Expansion:\n";
  cout << "  MaxBondOrderCE: " << GetMaxBondOrderCE() << "\n";
  cout << "  MaxClusterSizeCE: " << GetMaxClusterSizeCE() << "\n";

  cout << "  ElementSetCE: ";
  for (const auto &el : GetElementSetCE())
    cout << el << " ";
  cout << "\n";

  cout << "  ECIs: ";
  VectorXd ecis = GetECIs();
  for (int i = 0; i < ecis.size(); ++i)
    cout << setw(10) << ecis[i] << " ";
  cout << "\n\n";

  cout << "KRA:\n";
  cout << "  MaxBondOrderOfClusterKRA: " << GetMaxBondOrderOfClusterKRA() << "\n";
  cout << "  MaxBondOrderKRA: " << GetMaxBondOrderKRA() << "\n";
  cout << "  MaxClusterSizeKRA: " << GetMaxClusterSizeKRA() << "\n";

  cout << "  ElementSetKRA: ";
  for (const auto &el : GetElementSetKRA())
    cout << el << " ";
  cout << "\n";

  cout << "  KECIs:\n";
  VectorXd kecis = GetECIs();
  for (int i = 0; i < kecis.size(); ++i)
    cout << setw(10) << kecis[i] << " ";
  cout << "\n\n";

  cout << "------------------------------------------------------------" << endl;
}

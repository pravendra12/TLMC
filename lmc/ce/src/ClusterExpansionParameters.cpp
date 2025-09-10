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

size_t ClusterExpansionParameters::GetMaxBondOrder(
    const string &jsonKey) const
{
  string jsonSubKey = "maxBondOrder";

  size_t maxBondOrder = ReadSingleParameterFromJson<size_t>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  return maxBondOrder;
}

size_t ClusterExpansionParameters::GetMaxClusterSize(
    const string &jsonKey) const
{
  string jsonSubKey = "maxClusterSize";

  size_t maxClusterSize = ReadSingleParameterFromJson<size_t>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  return maxClusterSize;
}

vector<double> ClusterExpansionParameters::GetECIs(
    const string &ceType) const
{
  string jsonSubKey = "ecis";

  vector<double> eciVector = ReadParametersFromJson<double>(
      allParameters_,
      ceType, // ce or symCE
      jsonSubKey);

  return eciVector;
}

vector<double> ClusterExpansionParameters::GetClusterCutoffs() const
{
  string jsonKey = "symCE";
  string jsonSubKey = "clusterCutoffs";

  vector<double> clusterCutoffs = ReadParametersFromJson<double>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  return clusterCutoffs;
}

vector<string> ClusterExpansionParameters::GetAllowedElements() const
{
  string jsonKey = "symCE";
  string jsonSubKey = "allowedElements";

  vector<string> allowedElements = ReadParametersFromJson<string>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  return allowedElements;
}

unordered_map<string, double> ClusterExpansionParameters::GetChemicalPotentialsMap() const
{
  unordered_map<string, double> chemicalPotentialsMap;

  json chemPotJson = ReadParameterFromJson<json>(allParameters_["symCE"], "chemicalPotentials");

  // Iterate over each element in the JSON object
  for (auto &[key, value] : chemPotJson.items())
  {
    chemicalPotentialsMap[key] = value.get<double>();
  }

  return chemicalPotentialsMap;
}

set<Element> ClusterExpansionParameters::GetElementSet(
    const string &jsonKey) const
{
  // Works for kra and lvfe
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

Vector3d ClusterExpansionParameters::GetReferenceJumpDirection() const
{
  string jsonKey = "kra";
  string jsonSubKey = "referenceJumpDirection";

  vector<double> refJumpDirection = ReadParametersFromJson<double>(
      allParameters_,
      jsonKey,
      jsonSubKey);

  if (refJumpDirection.size() != 3)
  {
    throw runtime_error(
        "Error in `ClusterExpansionParameters::GetReferenceJumpDirection`: referenceJumpDirection must have exactly 3 elements, but JSON contains " + to_string(refJumpDirection.size()));
  }

  return Map<Vector3d>(refJumpDirection.data());
}

unordered_map<Element, VectorXd, boost::hash<Element>> ClusterExpansionParameters::GetKECIsMap() const
{
  unordered_map<Element, VectorXd, boost::hash<Element>> keciMap;

  json keciJson = ReadParameterFromJson<json>(allParameters_["kra"], "kecis");

  // Iterate over each element in the JSON object
  for (auto &[key, value] : keciJson.items())
  {
    vector<double> vec = value.get<vector<double>>();
    keciMap[Element{key}] = Map<VectorXd>(vec.data(), vec.size());
  }

  return keciMap;
}

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
  cout << "=== ClusterExpansionParameters Debug ===" << endl;

  try
  {
    cout << "Coefficient file: " << GetCoefficientFile() << endl;
    cout << "Basis type: " << GetBasisType() << endl;

    // CE
    cout << "\n-- CE --" << endl;
    cout << "Max bond order: " << GetMaxBondOrder("ce") << endl;
    cout << "Max cluster size: " << GetMaxClusterSize("ce") << endl;
    auto ceEcis = GetECIs("ce");
    cout << "ECIs: ";
    for (auto val : ceEcis)
      cout << val << " ";
    cout << endl;

    // LVFE
    cout << "\n-- LVFE --" << endl;
    cout << "Max bond order: " << GetMaxBondOrder("lvfe") << endl;
    cout << "Max cluster size: " << GetMaxClusterSize("lvfe") << endl;
    auto lvfeEcis = GetECIs("lvfe");
    cout << "ECIs: ";
    for (auto val : lvfeEcis)
      cout << val << " ";
    cout << endl;

    // SymCE
    cout << "\n-- SymCE --" << endl;
    auto symCutoffs = GetClusterCutoffs();
    cout << "Cluster cutoffs: ";
    for (auto val : symCutoffs)
      cout << val << " ";
    cout << endl;

    auto allowedElements = GetAllowedElements();
    cout << "Allowed elements: ";
    for (auto &ele : allowedElements)
      cout << ele << " ";
    cout << endl;

    auto chemPotentials = GetChemicalPotentialsMap();
    cout << "Chemical potentials: ";
    for (auto &[ele, val] : chemPotentials)
      cout << ele << "=" << val << " ";
    cout << endl;

    // KRA
    cout << "\n-- KRA --" << endl;
    cout << "Max bond order: " << GetMaxBondOrder("kra") << endl;
    cout << "Max cluster size: " << GetMaxClusterSize("kra") << endl;

    auto kraElementSet = GetElementSet("kra");
    cout << "KRA element set: ";
    for (const auto &ele : kraElementSet)
      cout << ele.GetElementString() << " ";
    cout << endl;

    auto refJumpDir = GetReferenceJumpDirection();
    cout << "Reference jump direction: "
              << refJumpDir.transpose() << endl;

    auto keciMap = GetKECIsMap();
    cout << "KECIs Map:" << endl;
    for (auto &[ele, vec] : keciMap)
    {
      cout << "  " << ele.GetElementString() << ": ";
      for (int i = 0; i < vec.size(); ++i)
        cout << setprecision(6) << vec[i] << " ";
      cout << endl;
    }
  }
  catch (const exception &e)
  {
    cerr << "Exception during DebugAllFunctions: " << e.what() << endl;
  }

  cout << "=== End Debug ===" << endl;
}

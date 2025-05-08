#include "SymmetrySpglib.h"

// step 1 : get all the symmetry operation
std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> GetSymmetryMatrixPure()
{
  // Override lattice, positions, and types for testing
  Eigen::Matrix3d cell;
  cell << 1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;

  double lattice[3][3];
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      lattice[i][j] = cell(j, i); // Transpose to row-major
    }
  }

  size_t num_atom = 2;
  double positionArray[2][3] = {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}};
  std::vector<int> types = {1, 1};

  const int max_operations = 96;
  int rotation[max_operations][3][3];
  double translation[max_operations][3];
  int num_sym = spg_get_symmetry(rotation, translation, max_operations,
                                 lattice, positionArray, types.data(),
                                 static_cast<int>(num_atom), 1e-4);

  std::cout << "Found " << num_sym << " symmetry operations.\n";

  // Store the symmetry operations in Eigen format
  vector<pair<Matrix3d, Vector3d>> sym_ops;
  for (int i = 0; i < num_sym; ++i)
  {
    // Convert rotation matrix to Eigen::Matrix3d
    Matrix3d rot;
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        rot(j, k) = rotation[i][j][k];
      }
    }

    // Convert translation vector to Eigen::Vector3d
    Vector3d trans(translation[i][0], translation[i][1], translation[i][2]);

    sym_ops.push_back({rot, trans});
  }

  return sym_ops;
}

// steps 2 : apply symmetry operation to lattice clusters

LatticeCluster apply_symmetry(const Config &config,
                              const LatticeCluster &latticeCluster)
{

  // Get symmmetry operation

  auto symOpsVector = GetSymmetryMatrixPure();


  LatticeCluster result;

  auto latticeIdVector = latticeCluster.GetLatticeIdVector();

  for (const auto latticeId : latticeIdVector)
  {
    Vector3d relativePos = config.GetRelativePositionOfLattice(latticeId);

    Vector3d new_pos;

    // apply symmetry operations

    for (auto symOps : symOpsVector)
    {
      // rotation

      new_pos = symOps.first * relativePos;

      new_pos += symOps.second;

      // cout << new_pos.transpose() << endl;
    }

    // result.push_back(new_pos);
  }
  return result;
}

vector<LatticeCluster> apply_symmetry_cluster(const Config &config,
                                              const LatticeCluster &latticeCluster)
{
  vector<LatticeCluster> transformedClusters;
  auto symOpsVector = GetSymmetryMatrixPure();
  Matrix3Xd relative_positions = config.GetRelativePositionMatrix();



  for (const auto &symOps : symOpsVector)
  {
    const auto &rotation = symOps.first;
    const auto &translation = symOps.second;

    vector<Vector3d> transformedPositions;
    auto latticeIdVector = latticeCluster.GetLatticeIdVector();

    for (const auto latticeId : latticeIdVector)
    {
      Vector3d relativePos = config.GetRelativePositionOfLattice(latticeId);
      Vector3d new_pos = rotation * relativePos + translation;

      // Wrap into [0,1) unit cell
      for (int i = 0; i < 3; ++i)
      {
        new_pos[i] = fmod(new_pos[i] + 1.0, 1.0);
      }

      transformedPositions.push_back(new_pos);

      cout << new_pos.transpose() << endl;
    }


    

    // convert teh transformed position back to lattice ids

    // Map transformed positions back to lattice IDs
    std::vector<int> new_cluster_ids;

    for (const Vector3d &trans_pos : transformedPositions) {
        int matched_id = -1;

        for (int i = 0; i < config.GetNumLattices(); ++i) {
            Vector3d ref_pos = relative_positions.col(i);
            Vector3d diff = trans_pos - ref_pos;

            // Account for periodicity
            for (int j = 0; j < 3; ++j)
                diff[j] -= round(diff[j]);

            if (diff.norm() < 1e-4) {
                matched_id = i;
                break;
            }
        }

        if (matched_id == -1) {
            std::cerr << "Could not map transformed position back to lattice site!" << std::endl;
            continue;
        }

        new_cluster_ids.push_back(matched_id);
    }


    

    // LatticeCluster ();

    // auto latticeClusterType = IdentifyLatticeClusterType(config, cluster);
    // latticeClusterHashset.emplace(latticeClusterType, cluster);

    // transformedCluster.SetPositions(transformedPositions);
    // transformedClusters.push_back(transformedCluster);
  }

  return transformedClusters;
}

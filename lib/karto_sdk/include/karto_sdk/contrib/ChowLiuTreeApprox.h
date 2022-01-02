#ifndef KARTO_SDK__CHOWLIUTREEAPPROX_H_
#define KARTO_SDK__CHOWLIUTREEAPPROX_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "karto_sdk/Karto.h"
#include "karto_sdk/Mapper.h"
#include "karto_sdk/Types.h"

namespace karto
{
namespace contrib
{

/** An uncertain, gaussian-distributed 2D pose. */
struct UncertainPose2
{
  Pose2 mean;
  Matrix3 covariance;
};

/**
 * Returns the target pose relative to the source pose,
 * accounting for their joint distribution covariance.
 */
UncertainPose2 ComputeRelativePose2(
    Pose2 source_pose, Pose2 target_pose,
    Eigen::Matrix<double, 6, 6> joint_pose_covariance);

/** Marginalizes a variable from a sparse information matrix. */
Eigen::SparseMatrix<double> ComputeMarginalInformationMatrix(
    const Eigen::SparseMatrix<double> & information_matrix,
    const Eigen::Index discarded_variable_index,
    const Eigen::Index variables_dimension);

/**
 * Computes a Chow Liu tree approximation to a given clique
 * (i.e. a graphical representation of joint probability distribution).
 */
std::vector<Edge<LocalizedRangeScan> *> ComputeChowLiuTreeApproximation(
  const std::vector<Vertex<LocalizedRangeScan> *> & clique,
  const Eigen::SparseMatrix<double> & covariance_matrix);

}  // namespace contrib
}  // namespace karto

#endif // KARTO_SDK__CHOWLIUTREEAPPROX_H_

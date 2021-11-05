/*
 * slam_toolbox
 * Copyright (c) 2019, Samsung Research America
 *
 * THE WORK (AS DEFINED BELOW) IS PROVIDED UNDER THE TERMS OF THIS CREATIVE
 * COMMONS PUBLIC LICENSE ("CCPL" OR "LICENSE"). THE WORK IS PROTECTED BY
 * COPYRIGHT AND/OR OTHER APPLICABLE LAW. ANY USE OF THE WORK OTHER THAN AS
 * AUTHORIZED UNDER THIS LICENSE OR COPYRIGHT LAW IS PROHIBITED.
 *
 * BY EXERCISING ANY RIGHTS TO THE WORK PROVIDED HERE, YOU ACCEPT AND AGREE TO
 * BE BOUND BY THE TERMS OF THIS LICENSE. THE LICENSOR GRANTS YOU THE RIGHTS
 * CONTAINED HERE IN CONSIDERATION OF YOUR ACCEPTANCE OF SUCH TERMS AND
 * CONDITIONS.
 *
 */

/* Author: Steven Macenski */

#ifndef SLAM_TOOLBOX__EXPERIMENTAL__SLAM_TOOLBOX_LIFELONG_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__SLAM_TOOLBOX_LIFELONG_HPP_

#include <memory>
#include "slam_toolbox/slam_toolbox_common.hpp"

namespace slam_toolbox
{

class LifelongSlamToolbox : public SlamToolbox
{
public:
  explicit LifelongSlamToolbox(rclcpp::NodeOptions options);
  ~LifelongSlamToolbox() {}

  // computation metrics
  double computeObjectiveScore(
    const double & intersect_over_union, const double & area_overlap,
    const double & reading_overlap, const int & num_constraints,
    const double & initial_score, const int & num_candidates) const;
  static double computeIntersect(LocalizedRangeScan * s1, LocalizedRangeScan * s2);
  static double computeIntersectOverUnion(LocalizedRangeScan * s1, LocalizedRangeScan * s2);
  static double computeAreaOverlapRatio(
    LocalizedRangeScan * ref_scan,
    LocalizedRangeScan * candidate_scan);
  static double computeReadingOverlapRatio(
    LocalizedRangeScan * ref_scan,
    LocalizedRangeScan * candidate_scan);
  static void computeIntersectBounds(
    LocalizedRangeScan * s1, LocalizedRangeScan * s2, double & x_l,
    double & x_u, double & y_l, double & y_u);

protected:
  void laserCallback(
    sensor_msgs::msg::LaserScan::ConstSharedPtr scan) override;
  bool deserializePoseGraphCallback(
    const std::shared_ptr<rmw_request_id_t> request_header,
    const std::shared_ptr<slam_toolbox::srv::DeserializePoseGraph::Request> req,
    std::shared_ptr<slam_toolbox::srv::DeserializePoseGraph::Response> resp) override;

  void evaluateNodeDepreciation(LocalizedRangeScan * range_scan);
  void removeFromSlamGraph(Vertex<LocalizedRangeScan> * vertex);
  double computeScore(
    LocalizedRangeScan * reference_scan, Vertex<LocalizedRangeScan> * candidate,
    const double & initial_score, const int & num_candidates);
  ScoredVertices computeScores(Vertices & near_scans, LocalizedRangeScan * range_scan);
  Vertices FindScansWithinRadius(LocalizedRangeScan * scan, const double & radius);
  void updateScoresSlamGraph(const double & score, Vertex<LocalizedRangeScan> * vertex);
  void checkIsNotNormalized(const double & value);

  /*****************************************************************************/
  /*******************************Implementation********************************/
  std::vector<int> getGridPosition(float x, float y, float resolution);
  std::vector<float> getLaserHit(std::vector<float> const& robot_pose, float range, float angle);
  std::pair<std::vector<int>, std::vector<int>> Bresenham(int x_1, int y_1, int x_2, int y_2);
  std::vector<float> calculateIntersection(std::vector<float> laser_start, std::vector<float> laser_end, std::vector<float> cell_start, std::vector<float> cell_end);
  int getSign(int n_1, int n_2);
  float calculateProbability(float range_1, float range_2);
  float calculateDistance(float x_1, float y_1, float x_2, float y_2);
  void scannerTest();
  std::vector<float> getCellPosition(std::vector<int> grid_cell, float resolution);
  /*****************************************************************************/

  bool use_tree_;
  double iou_thresh_;
  double removal_score_;
  double overlap_scale_;
  double constraint_scale_;
  double candidates_scale_;
  double iou_match_;
  double nearby_penalty_;
};

}  // namespace slam_toolbox

#endif  // SLAM_TOOLBOX__EXPERIMENTAL__SLAM_TOOLBOX_LIFELONG_HPP_

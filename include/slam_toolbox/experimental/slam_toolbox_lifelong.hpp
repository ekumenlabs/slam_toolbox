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
#include <cmath>
#include <math.h>
#include <tuple>
#include "slam_toolbox/slam_toolbox_common.hpp"
#include "slam_toolbox/experimental/theoretic_information.hpp"

namespace slam_toolbox
{

class LifelongSlamToolbox : public SlamToolbox
{
  typedef std::tuple<int, int, int> map_tuple;
  typedef std::pair<map_tuple, float> map_pair;

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

public:
  // Cell occupancy struct for unordered map

  struct HashTuple
  {
      std::size_t operator() (map_tuple const& key) const
      {
          std::size_t hash = 5381u;
          hash = (hash << 5) + hash + std::get<0>(key);
          hash = (hash << 5) + hash + std::get<1>(key);
          hash = (hash << 5) + hash + std::get<2>(key);
          return hash;
      }
  };
  struct Occupancy
  {
    int fr, oc ,un;
    bool operator==(Occupancy const& st) const
    {
      return (st.fr == fr) && (st.oc == oc) && (st.un == un);
    }
    
    struct CombinationsHash
    {
      std::size_t operator()(Occupancy const& key) const 
      {
        std::size_t hash = 5381u;
        hash = (hash << 5) + hash + key.fr;
        hash = (hash << 5) + hash + key.oc;
        hash = (hash << 5) + hash + key.un;
        return hash;
      }
    };
  };

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

  /*******************************Implementation********************************/
  void scannerTest();

  // Grid operations
  void initializeGrids();
  void calculateLimits(
    std::vector<float> & initial_x, std::vector<float> & initial_y, std::vector<float> & final_x, std::vector<float> & final_y, 
    float & limit_x, float & limit_y, float & min_x, float & max_x, float & min_y, float & max_y, std::vector<int> & robot_grid_pos, 
    std::vector<int> & final_grid_pos);

  // Grid and position information
  std::pair<std::vector<int>, std::vector<int>> Bresenham(int x_1, int y_1, int x_2, int y_2);
  std::vector<int> getGridPosition(float x, float y);
  std::vector<float> laserHitDistance(std::vector<float> const& robot_pose, float range, float angle);
  std::vector<float> calculateIntersection(std::vector<float> laser_start, std::vector<float> laser_end, std::vector<float> cell_start, std::vector<float> cell_end);
  int getSign(int n_1, int n_2);

  // Measurements calculations <P(free), P(Occ), P(Unk)>
  float probabilityFromObservation(float range_1, float range_2);
  float euclideanDistance(float x_1, float y_1, float x_2, float y_2);

  // Mutual information 
  float measurementOutcomeEntropy(map_tuple const& meas_outcome);

  void recoverProbability();
  float calculateLogs(float probability);
  float probabilityFromLogs(float log);
  float calculateEntropy(float probability);

  float calculateMapMutualInformation();
  void updateCellMutualInformation(float mut_inf_val);

  // Measurement outcomes probabilities
  void appendCellProbabilities(std::vector<float>& meas_outcomes);
  void computeProbabilities(std::vector<std::vector<float>>& meas_outcm);
  std::vector<std::vector<float>> retreiveMeasurementOutcomes();
  std::vector<int> unhashIndex(int hash);

  // Data structures 
  std::unordered_map<map_tuple, float, HashTuple> m_map_out;
  std::unordered_map<Occupancy, float, Occupancy::CombinationsHash> m_un_cmb;
  std::map<std::vector<int>, std::vector<std::vector<float>>> m_cell_probabilities;
  std::vector<std::vector<float>> m_mutual_grid;
  std::vector<std::vector<int>> m_grid;
  float m_map_dist;
  float m_resolution;
  int m_cell_x;
  int m_cell_y;
  int m_num_cells;

  // Robot information - Defining values just for testing
  std::vector<std::vector<float>> robot_poses {{5.6f, 6.0f, M_PI/2}, {3.5f, 9.0f, 0.0f}};
  std::vector<std::vector<float>> laser_ranges {{1.65f, 5.0f, 5.0f, 5.0f, 5.0f}, {5.0f, 5.0f, 4.0f, 5.0f, 5.0f}};
  std::vector<float> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

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

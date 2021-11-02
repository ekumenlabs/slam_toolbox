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

#include <algorithm>
#include <memory>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include "slam_toolbox/experimental/slam_toolbox_lifelong.hpp"

#define PI 3.14159265

namespace slam_toolbox
{

/*****************************************************************************/
void LifelongSlamToolbox::checkIsNotNormalized(const double & value)
/*****************************************************************************/
{
  if (value < 0.0 || value > 1.0) {
    RCLCPP_FATAL(get_logger(),
      "All stores and scales must be in range [0, 1].");
    exit(-1);
  }
}

/*****************************************************************************/
LifelongSlamToolbox::LifelongSlamToolbox(rclcpp::NodeOptions options)
: SlamToolbox(options)
/*****************************************************************************/
{
  use_tree_ = false;
  use_tree_ = this->declare_parameter("lifelong_search_use_tree", use_tree_);
  iou_thresh_ = 0.10;
  iou_thresh_ = this->declare_parameter("lifelong_minimum_score", iou_thresh_);
  iou_match_ = 0.85;
  iou_match_ = this->declare_parameter("lifelong_iou_match", iou_match_);
  removal_score_ = 0.10;
  removal_score_ = this->declare_parameter("lifelong_node_removal_score",
      removal_score_);
  overlap_scale_ = 0.5;
  overlap_scale_ = this->declare_parameter("lifelong_overlap_score_scale",
      overlap_scale_);
  constraint_scale_ = 0.05;
  constraint_scale_ = this->declare_parameter("lifelong_constraint_multiplier",
      constraint_scale_);
  nearby_penalty_ = 0.001;
  nearby_penalty_ = this->declare_parameter("lifelong_nearby_penalty",
      nearby_penalty_);
  candidates_scale_ = 0.03;
  candidates_scale_ = this->declare_parameter("lifelong_candidates_scale",
      candidates_scale_);

  checkIsNotNormalized(iou_thresh_);
  checkIsNotNormalized(constraint_scale_);
  checkIsNotNormalized(removal_score_);
  checkIsNotNormalized(overlap_scale_);
  checkIsNotNormalized(iou_match_);
  checkIsNotNormalized(nearby_penalty_);
  checkIsNotNormalized(candidates_scale_);

  RCLCPP_WARN(get_logger(), "Lifelong mapping mode in SLAM Toolbox is considered "
    "experimental and should be understood before proceeding. Please visit: "
    "https://github.com/SteveMacenski/slam_toolbox/wiki/"
    "Experimental-Lifelong-Mapping-Node for more information.");

  // in lifelong mode, we cannot have interactive mode enabled
  enable_interactive_mode_ = false;

  scannerTest();
}

/*****************************************************************************/
/**
  * The are implementation notes from Camilo OCT 21
  * 
  * First step
    * Get the occupancy probability for each cell in the map. This solution is based
    * on the Log-Odds rather than the Bayesian filter, only for computation performance
    * Initialization for this elements should be made before we start the following loop:
      * Gather the sensor data. 
      * Inverse measurement model (Tranforming scan to measure occupancy belief map)
      * Log-Odds calculation for the current reading and update the grid or map.
      * Probability map calculation and update for current reading (This is done for the cells)\
  * 
  * NOTE
  * Once we get the probabilities we can calculate the entropy of each grid cell. I am assuming
  * this should be done inside the loop as we want to calculate this entropy at each time.  

  * Log odds assume values from −∞ to ∞. The Bayes filter for updating beliefs in
  * log odds representation is computationally elegant. It avoids truncation problems that
  * arise for probabilities close to 0 or 1.
  * 
  * 
  * We are using std::vector !!
  * 
  * 
  * 
  * According tothe the last meeting we handle, I will start with a prrof of concept
  * Create a grid
  * Calculate the scan position (WRO a given cell). Consider where is my initial cell 
  * Calculate the probability of seing an grid given a set of measurements
  * Create the histogram with Algrotihm 1
/*****************************************************************************/

/*****************************************************************************/
/*******************************Implementation********************************/
  
void LifelongSlamToolbox::scannerTest()
{
  RCLCPP_WARN(get_logger(), "<-------- Scanner Test -------->");
 
  float resolution = 0.5f; // Cell resolution - 1 meter
  float map_dist = 20.0f; // Total map distance
  
  int number_cells = map_dist / resolution;
  std::vector<std::vector<int>> grid(number_cells); // 20 meters in Y dimension
  for (int i=0; i<number_cells; ++i)
  {
    grid[i].resize(number_cells); // 20 meters in X dimension
    for (int j=0; j<number_cells; ++j)
    {
      grid[i][j] = 0;
    }
  }

  // I will have 5 lasers for each reading (Located at 0, +-25, +-50)
  std::vector<float> robot_pose{5.5f, 6.0f, -PI/2};

  // This is the initial point
  std::vector<int> robot_grid_pos = getGridPosition(robot_pose[0], robot_pose[1], resolution);
  std::cout << "Robot position: " << robot_grid_pos[0] << ", " << robot_grid_pos[1] << std::endl;

  // Creating the laser scan with 5 beams ------- Angles will be -55, -25, 0, 25, 55 //----// +-55 = +-0.87266 : +-25 = +-0.43633
  // std::vector<float> ranges{5.0f, 5.0f, 5.0f, 5.0f, 5.0f}; // Maximum sensor range is 5 meters
  // std::vector<float> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

  std::vector<float> ranges{5.0f, 5.0f, 5.0f, 5.0f, 5.0f}; // Maximum sensor range is 5 meters
  std::vector<float> angles{0.0f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

  // Current yaw + beam angle: -PI/2 (-1.570795) -0.87266 = 2.44345 (-55 degrees)
  // for (int i = 0; i < ranges.size(); ++i)
  for (int i = 0; i < 1; ++i) // One reading only
  {
    std::cout << "........ New laser ........" << std::endl;
    std::cout << "Distance: " << ranges[i] << ", Angle: " << angles[i] << std::endl;

    std::vector<float> laser_grid = getLaserHit(robot_pose, ranges[i], angles[i]);
    std::vector<int> final_grid_pos = getGridPosition(laser_grid[0], laser_grid[1], resolution);
    std::cout << "Laser end: " << final_grid_pos[0] << ", " << final_grid_pos[1] << std::endl;
    // robot_grid_pos[0] // X1 - robot_grid_pos[1] // Y1
    // final_grid_pos[0] // X2 - final_grid_pos[1] // Y2
    std::vector<int> cells_x, cells_y;
    std::pair<std::vector<int>, std::vector<int>> res_pair = Bresenham(robot_grid_pos[0], robot_grid_pos[1], final_grid_pos[0], final_grid_pos[1]);
    /**
     * Need to append here the last point, the algorithm is not doing it
    */

    // Cells visited by this laser beam
    cells_x = res_pair.first;
    cells_y = res_pair.second;

    // Laser beam discretization
    int samples = 100;
    float step = ranges[i] / samples;
    float current_dis = 0.0f;

    // Vectors for discretization of the laser beam
    std::vector<float> x_coord, y_coord;
    x_coord.reserve(samples + 1);
    y_coord.reserve(samples + 1);

    while(current_dis < (ranges[i] + step))
    {
      // Discretization loop
      float x_p = current_dis * sin(angles[i]);
      float y_p = current_dis * cos(angles[i]);
      x_coord.emplace_back(x_p);
      y_coord.emplace_back(y_p);
      current_dis += step;
    }

    int count_idx = 0;
    float dist_carry = 0.0f;

    std::vector<float> initial_point(2), final_point(2);
      
    initial_point[0] = 0.0f; 
    initial_point[1] = 0.0f;

    // Move alongside all the cells that the laser beam hits
    for (int i = 0; i < cells_x.size(); ++i)
    {
      // std::cout << cells_x[i] << ", " << cells_y[i] << std::endl;
      // Converting the cell into distance (Global)
      float limit_x = cells_x[i] * resolution;
      float limit_y = cells_y[i] * resolution;

      std::cout << "******** New cell ********" << std::endl;
      
      std::cout << "Initial point: " << initial_point[0] << ", " << initial_point[1] << std::endl;      
      std::cout << "starting at: " << count_idx << std::endl;
      
      for (int j = count_idx; j < x_coord.size() + 1; ++j)
      {
        // Adding the robot position to the laser beam reading - Shift X and Y
        float read_x = x_coord[j] + robot_pose[0];
        float read_y = y_coord[j] + robot_pose[1];

        // std::cout << "Limit X: " << limit_x << ", " << limit_x + resolution << std::endl;
        // std::cout << "Limit Y: " << limit_y << ", " << limit_y + resolution << std::endl;
        std::cout << "Count: " << j << std::endl;

        // Evaluate to what cell corresponds the current reading
        if (((abs(read_x) >= abs(limit_x)) && (abs(read_x) <= abs(limit_x + resolution))) &&
            ((abs(read_y) >= abs(limit_y)) && (abs(read_y) <= abs(limit_y + resolution))))
        {
          // std::cout << "Reading: " << read_x << ", " << read_y << std::endl;
          // std::cout << "In cell: " << cells_x[i] << ", " << cells_y[i] << std::endl;
          final_point[0] = x_coord[j]; 
          final_point[1] = y_coord[j];
          // std::cout << "Count: " << j << std::endl;
          ++count_idx;

        }
        else 
        {
          // This one is working 
          break;
        }        
      }
      std::cout << "Final point: " << final_point[0] << ", " << final_point[1] << std::endl;
      // From the start of the cell to the end - This is an approximation - From the robot pose
      float distance = calculateDistance(initial_point[0], initial_point[1], final_point[0], final_point[1]);
      dist_carry += distance;
      std::cout << "Distance: " << distance << ", " << dist_carry << std::endl;

      // Final point becomes the initial for the next movement  
      initial_point[0] = final_point[0]; 
      initial_point[1] = final_point[1];

    }
    std::cout << "*************************" << std::endl;

  }
}


float LifelongSlamToolbox::calculateDistance(float x_1, float y_1, float x_2, float y_2)
{
  /* 
    Calculates the euclidean distance between two points
  */
  float diff_x = x_2 - x_1;
  float diff_y = y_2 - y_1;

  return sqrt(diff_x*diff_x + diff_y*diff_y);
}

float LifelongSlamToolbox::calculateProbability(float range)
{
  /* 
    Calculates the probability of a cell being observed by a given measurement
  */
  float max_range = 5.0f;
  float lambda = 0.285f;
  float nu = 1.0f / lambda;
  if (range <= max_range)
  {
    return nu*lambda*exp(-lambda*range);
  }
  return 0.0f;
}

std::vector<int> LifelongSlamToolbox::getGridPosition(float x, float y, float resolution)
{
  /* 
    Maps the distance into grid coordinates 
  */
  int x_cell = ceil((1 / resolution) * x); 
  int y_cell = ceil((1 / resolution) * y); 

  return {x_cell, y_cell};
}

std::vector<float> LifelongSlamToolbox::getLaserHit(std::vector<float> const& robot_pose, float distance, float angle)
{
  /* 
    Returns the distance where the laser beam hits something 
  */ 
  float angle_r = atan2(sin(robot_pose[2] + angle), cos(robot_pose[2] + angle));
  float x_occ = distance * cos(angle_r) + robot_pose[0]; // This is X
  float y_occ = -distance * sin(angle_r) + robot_pose[1]; // This is Y

  std::vector<int> grid_pos = getGridPosition(x_occ, y_occ, 0.5f);

  return {x_occ, y_occ};
}


std::pair<std::vector<int>, std::vector<int>> LifelongSlamToolbox::Bresenham(int x_1, int y_1, int x_2, int y_2)
{
  /* 
    Returns the set of cells hit by a laser beam 
  */
	std::vector<int> x_bres;
	std::vector<int> y_bres;

	int x = x_1;
	int y = y_1;
	
	int delta_x = abs(x_2 - x_1);
	int delta_y = abs(y_2 - y_1);

  int s_x = getSign(x_1, x_2);
  int s_y = getSign(y_1, y_2);
  bool interchange = false;

	if (delta_y > delta_x)
  {
    int temp = delta_x;
    delta_x = delta_y;
    delta_y = temp;    
		interchange = true;
  }
  else
  {
		interchange = false;
  }

  int a_res = 2 * delta_y;
	int b_res = 2 * (delta_y - delta_x);
	int e_res = (2 * delta_y) - delta_x;

  x_bres.push_back(x);
  y_bres.push_back(y);

  for (int i = 1; i < delta_x; ++i) 
  {
    if (e_res < 0)
    {
      if (interchange)
      {
        y += s_y;
      }
      else
      {
        x += s_x;
      }
      e_res += a_res;
    }
    else 
    {
      y += s_y;
      x += s_x;
      e_res += b_res;
    }
    x_bres.push_back(x);
    y_bres.push_back(y);
  }
  return std::pair<std::vector<int>, std::vector<int>>{x_bres, y_bres};
}

int LifelongSlamToolbox::getSign(int n_1, int n_2)
{
  /* 
    Returns the sign of an operation, used for Bresenham algorithm
  */
  int difference = n_2 - n_1;
  
  if (difference == 0) { return 0; }
  else if (difference < 0) { return -1; }
  else { return 1; }
}
/*****************************************************************************/


/*****************************************************************************/
void LifelongSlamToolbox::laserCallback(
  sensor_msgs::msg::LaserScan::ConstSharedPtr scan)
/*****************************************************************************/
{
  // no odom info
  Pose2 pose;
  if (!pose_helper_->getOdomPose(pose, scan->header.stamp)) {
    RCLCPP_WARN(get_logger(), "Failed to compute odom pose");
    return;
  }

  // ensure the laser can be used
  LaserRangeFinder * laser = getLaser(scan);

  if (!laser) {
    RCLCPP_WARN(get_logger(), "Failed to create laser device for"
      " %s; discarding scan", scan->header.frame_id.c_str());
    return;
  }

  // LTS additional bounded node increase parameter (rate, or total for run or at all?)
  // LTS pseudo-localization mode. If want to add a scan, but
  // not deleting a scan, add to local buffer?
  // LTS if (eval() && dont_add_more_scans) {addScan()} else {localization_add_scan()}
  // LTS if (eval() && ctr / total < add_rate_scans) {addScan()} else {localization_add_scan()}
  LocalizedRangeScan * range_scan = addScan(laser, scan, pose);
  evaluateNodeDepreciation(range_scan);
}

/*****************************************************************************/
void LifelongSlamToolbox::evaluateNodeDepreciation(
  LocalizedRangeScan * range_scan)
/*****************************************************************************/
{
  if (range_scan) {
    boost::mutex::scoped_lock lock(smapper_mutex_);

    const BoundingBox2 & bb = range_scan->GetBoundingBox();
    const Size2<double> bb_size = bb.GetSize();
    double radius = sqrt(bb_size.GetWidth() * bb_size.GetWidth() +
        bb_size.GetHeight() * bb_size.GetHeight()) / 2.0;
    Vertices near_scan_vertices = FindScansWithinRadius(range_scan, radius);

    ScoredVertices scored_verices =
      computeScores(near_scan_vertices, range_scan);

    ScoredVertices::iterator it;
    for (it = scored_verices.begin(); it != scored_verices.end(); ++it) {
      if (it->GetScore() < removal_score_) {
        RCLCPP_DEBUG(get_logger(),
          "Removing node %i from graph with score: %f and old score: %f.",
          it->GetVertex()->GetObject()->GetUniqueId(),
          it->GetScore(), it->GetVertex()->GetScore());
        removeFromSlamGraph(it->GetVertex());
      } else {
        updateScoresSlamGraph(it->GetScore(), it->GetVertex());
      }
    }
  }
}

/*****************************************************************************/
Vertices LifelongSlamToolbox::FindScansWithinRadius(
  LocalizedRangeScan * scan, const double & radius)
/*****************************************************************************/
{
  // Using the tree will create a Kd-tree and find all neighbors in graph
  // around the reference scan. Otherwise it will use the graph and
  // access scans within radius that are connected with constraints to this
  // node.

  if (use_tree_) {
    return
      smapper_->getMapper()->GetGraph()->FindNearByVertices(
      scan->GetSensorName(), scan->GetBarycenterPose(), radius);
  } else {
    return
      smapper_->getMapper()->GetGraph()->FindNearLinkedVertices(scan, radius);
  }
}

/*****************************************************************************/
double LifelongSlamToolbox::computeObjectiveScore(
  const double & intersect_over_union,
  const double & area_overlap,
  const double & reading_overlap,
  const int & num_constraints,
  const double & initial_score,
  const int & num_candidates) const
/*****************************************************************************/
{
  // We have some useful metrics. lets make a new score
  // intersect_over_union: The higher this score, the better aligned in scope these scans are
  // area_overlap: The higher, the more area they share normalized by candidate size
  // reading_overlap: The higher, the more readings of the new scan the candidate contains
  // num_constraints: The lower, the less other nodes may rely on this candidate
  // initial_score: Last score of this vertex before update

  // this is a really good fit and not from a loop closure, lets just decay
  if (intersect_over_union > iou_match_ && num_constraints < 3) {
    return -1.0;
  }

  // to be conservative, lets say the overlap is the lesser of the
  // area and proportion of laser returns in the intersecting region.
  double overlap = overlap_scale_ * std::min(area_overlap, reading_overlap);

  // if the num_constraints are high we want to stave off the decay, but not override it
  double contraint_scale_factor = std::min(1.0,
      std::max(0., constraint_scale_ * (num_constraints - 2)));
  contraint_scale_factor = std::min(contraint_scale_factor, overlap);

  //
  double candidates = num_candidates - 1;
  double candidate_scale_factor = candidates_scale_ * candidates;

  // Give the initial score a boost proportional to the number of constraints
  // Subtract the overlap amount, apply a penalty for just being nearby
  // and scale the entire additional score by the number of candidates
  double score =
    initial_score * (1.0 + contraint_scale_factor) -
    overlap -
    nearby_penalty_;

  // score += (initial_score - score) * candidate_scale_factor;

  if (score > 1.0) {
    RCLCPP_ERROR(get_logger(),
      "Objective function calculated for vertex score (%0.4f)"
      " greater than one! Thresholding to 1.0", score);
    return 1.0;
  }

  return score;
}

/*****************************************************************************/
double LifelongSlamToolbox::computeScore(
  LocalizedRangeScan * reference_scan,
  Vertex<LocalizedRangeScan> * candidate,
  const double & initial_score, const int & num_candidates)
/*****************************************************************************/
{
  double new_score = initial_score;
  LocalizedRangeScan * candidate_scan = candidate->GetObject();

  // compute metrics for information loss normalized
  double iou = computeIntersectOverUnion(reference_scan, candidate_scan);
  double area_overlap = computeAreaOverlapRatio(reference_scan, candidate_scan);
  int num_constraints = candidate->GetEdges().size();
  double reading_overlap = computeReadingOverlapRatio(reference_scan, candidate_scan);

  bool critical_lynchpoint = candidate_scan->GetUniqueId() == 0 ||
    candidate_scan->GetUniqueId() == 1;
  int id_diff = reference_scan->GetUniqueId() - candidate_scan->GetUniqueId();
  if (id_diff < smapper_->getMapper()->getParamScanBufferSize() ||
    critical_lynchpoint)
  {
    return initial_score;
  }

  double score = computeObjectiveScore(iou,
      area_overlap,
      reading_overlap,
      num_constraints,
      initial_score,
      num_candidates);

  RCLCPP_INFO(get_logger(), "Metric Scores: Initial: %f, IOU: %f,"
    " Area: %f, Num Con: %i, Reading: %f, outcome score: %f.",
    initial_score, iou, area_overlap, num_constraints, reading_overlap, score);
  return score;
}

/*****************************************************************************/
ScoredVertices LifelongSlamToolbox::computeScores(
  Vertices & near_scans,
  LocalizedRangeScan * range_scan)
/*****************************************************************************/
{
  ScoredVertices scored_vertices;
  scored_vertices.reserve(near_scans.size());

  // must have some minimum metric to utilize
  // IOU will drop sharply with fitment, I'd advise not setting this value
  // any higher than 0.15. Also check this is a linked constraint
  // We want to do this early to get a better estimate of local candidates
  ScanVector::iterator candidate_scan_it;
  double iou = 0.0;
  for (candidate_scan_it = near_scans.begin();
    candidate_scan_it != near_scans.end(); )
  {
    iou = computeIntersectOverUnion(range_scan,
        (*candidate_scan_it)->GetObject());
    if (iou < iou_thresh_ || (*candidate_scan_it)->GetEdges().size() < 2) {
      candidate_scan_it = near_scans.erase(candidate_scan_it);
    } else {
      ++candidate_scan_it;
    }
  }

  for (candidate_scan_it = near_scans.begin();
    candidate_scan_it != near_scans.end(); ++candidate_scan_it)
  {
    ScoredVertex scored_vertex((*candidate_scan_it),
      computeScore(range_scan, (*candidate_scan_it),
      (*candidate_scan_it)->GetScore(), near_scans.size()));
    scored_vertices.push_back(scored_vertex);
  }
  return scored_vertices;
}

/*****************************************************************************/
void LifelongSlamToolbox::removeFromSlamGraph(
  Vertex<LocalizedRangeScan> * vertex)
/*****************************************************************************/
{
  smapper_->getMapper()->RemoveNodeFromGraph(vertex);
  smapper_->getMapper()->GetMapperSensorManager()->RemoveScan(
    vertex->GetObject());
  dataset_->RemoveData(vertex->GetObject());
  vertex->RemoveObject();
  delete vertex;
  vertex = nullptr;
  // LTS what do we do about the contraints that node had about it?Nothing?Transfer?
}

/*****************************************************************************/
void LifelongSlamToolbox::updateScoresSlamGraph(
  const double & score,
  Vertex<LocalizedRangeScan> * vertex)
/*****************************************************************************/
{
  // Saved in graph so it persists between sessions and runs
  vertex->SetScore(score);
}

/*****************************************************************************/
bool LifelongSlamToolbox::deserializePoseGraphCallback(
  const std::shared_ptr<rmw_request_id_t> request_header,
  const std::shared_ptr<slam_toolbox::srv::DeserializePoseGraph::Request> req,
  std::shared_ptr<slam_toolbox::srv::DeserializePoseGraph::Response> resp)
/*****************************************************************************/
{
  if (req->match_type == procType::LOCALIZE_AT_POSE) {
    RCLCPP_ERROR(get_logger(), "Requested a localization deserialization "
      "in non-localization mode.");
    return false;
  }

  return SlamToolbox::deserializePoseGraphCallback(request_header, req, resp);
}

/*****************************************************************************/
void LifelongSlamToolbox::computeIntersectBounds(
  LocalizedRangeScan * s1, LocalizedRangeScan * s2,
  double & x_l, double & x_u, double & y_l, double & y_u)
/*****************************************************************************/
{
  Size2<double> bb1 = s1->GetBoundingBox().GetSize();
  Size2<double> bb2 = s2->GetBoundingBox().GetSize();
  Pose2 pose1 = s1->GetBarycenterPose();
  Pose2 pose2 = s2->GetBarycenterPose();

  const double s1_upper_x = pose1.GetX() + (bb1.GetWidth() / 2.0);
  const double s1_upper_y = pose1.GetY() + (bb1.GetHeight() / 2.0);
  const double s1_lower_x = pose1.GetX() - (bb1.GetWidth() / 2.0);
  const double s1_lower_y = pose1.GetY() - (bb1.GetHeight() / 2.0);

  const double s2_upper_x = pose2.GetX() + (bb2.GetWidth() / 2.0);
  const double s2_upper_y = pose2.GetY() + (bb2.GetHeight() / 2.0);
  const double s2_lower_x = pose2.GetX() - (bb2.GetWidth() / 2.0);
  const double s2_lower_y = pose2.GetY() - (bb2.GetHeight() / 2.0);

  x_u = std::min(s1_upper_x, s2_upper_x);
  y_u = std::min(s1_upper_y, s2_upper_y);
  x_l = std::max(s1_lower_x, s2_lower_x);
  y_l = std::max(s1_lower_y, s2_lower_y);
}

/*****************************************************************************/
double LifelongSlamToolbox::computeIntersect(
  LocalizedRangeScan * s1,
  LocalizedRangeScan * s2)
/*****************************************************************************/
{
  double x_l, x_u, y_l, y_u;
  computeIntersectBounds(s1, s2, x_l, x_u, y_l, y_u);
  const double intersect = (y_u - y_l) * (x_u - x_l);

  if (intersect < 0.0) {
    return 0.0;
  }

  return intersect;
}

/*****************************************************************************/
double LifelongSlamToolbox::computeIntersectOverUnion(
  LocalizedRangeScan * s1,
  LocalizedRangeScan * s2)
/*****************************************************************************/
{
  // this is a common metric in machine learning used to determine
  // the fitment of a set of bounding boxes. Its response sharply
  // drops by box matches.

  const double intersect = computeIntersect(s1, s2);

  Size2<double> bb1 = s1->GetBoundingBox().GetSize();
  Size2<double> bb2 = s2->GetBoundingBox().GetSize();
  const double uni = (bb1.GetWidth() * bb1.GetHeight()) +
    (bb2.GetWidth() * bb2.GetHeight()) - intersect;

  return intersect / uni;
}

/*****************************************************************************/
double LifelongSlamToolbox::computeAreaOverlapRatio(
  LocalizedRangeScan * ref_scan,
  LocalizedRangeScan * candidate_scan)
/*****************************************************************************/
{
  // ref scan is new scan, candidate scan is potential for decay
  // so we want to find the ratio of space of the candidate scan
  // the reference scan takes up

  double overlap_area = computeIntersect(ref_scan, candidate_scan);
  Size2<double> bb_candidate = candidate_scan->GetBoundingBox().GetSize();
  const double candidate_area =
    bb_candidate.GetHeight() * bb_candidate.GetWidth();

  return overlap_area / candidate_area;
}

/*****************************************************************************/
double LifelongSlamToolbox::computeReadingOverlapRatio(
  LocalizedRangeScan * ref_scan,
  LocalizedRangeScan * candidate_scan)
/*****************************************************************************/
{
  const PointVectorDouble & pts = candidate_scan->GetPointReadings(true);
  const int num_pts = pts.size();

  // get the bounds of the intersect area
  double x_l, x_u, y_l, y_u;
  computeIntersectBounds(ref_scan, candidate_scan, x_l, x_u, y_l, y_u);

  PointVectorDouble::const_iterator pt_it;
  int inner_pts = 0;
  for (pt_it = pts.begin(); pt_it != pts.end(); ++pt_it) {
    if (pt_it->GetX() < x_u && pt_it->GetX() > x_l &&
      pt_it->GetY() < y_u && pt_it->GetY() > y_l)
    {
      inner_pts++;
    }
  }

  return static_cast<double>(inner_pts) / static_cast<double>(num_pts);
}

}  // namespace slam_toolbox

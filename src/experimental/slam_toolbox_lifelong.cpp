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

  float unknown_prob = 0.5f;
  float initial_log = calculateLogs(unknown_prob); 
  std::vector<std::vector<int>> grid(number_cells); // 20 meters in Y dimension
  std::vector<std::vector<float>> grid_prob(number_cells); // 20 meters in Y dimension
  std::vector<std::vector<float>> grid_logs(number_cells); // 20 meters in Y dimension

  for (int i=0; i<number_cells; ++i)
  {
    grid[i].resize(number_cells); // 20 meters in X dimension
    grid_prob[i].resize(number_cells);
    grid_logs[i].resize(number_cells);
    for (int j=0; j<number_cells; ++j)
    {
      grid[i][j] = 0;
      grid_prob[i][j] = unknown_prob; 
      grid_logs[i][j] = initial_log;
    }
  }

  std::cout << "Logs---------------: " << grid_logs[5][25] << std::endl;

  // I will have 5 lasers for each reading (Located at 0, +-25, +-50)
  // std::vector<float> robot_pose{5.6f, 6.0f, PI/2};
  // std::vector<float> robot_pose{5.6f, 6.0f, PI/4};
  std::vector<float> robot_pose{5.6f, 6.0f, PI/2};
  // std::vector<float> robot_pose{5.5f, 6.0f, -PI/4};
  // std::vector<float> robot_pose{5.5f, 6.0f, -3*PI/4};
  // -PI is an unhandled case
  // std::vector<float> robot_pose{5.5f, 6.0f, PI};
  // std::vector<float> robot_pose{5.5f, 6.0f, 0.0f};

  // This is the initial point
  std::vector<int> robot_grid_pos = getGridPosition(robot_pose[0], robot_pose[1], resolution);
  std::cout << "Robot position: " << robot_grid_pos[0] << ", " << robot_grid_pos[1] << std::endl;

  // Creating the laser scan with 5 beams ------- Angles will be -50, -25, 0, 25, 50 //----// +-50 = +-0.87266 : +-25 = +-0.43633
  std::vector<float> ranges{1.65f, 5.0f, 5.0f, 5.0f, 5.0f}; // Maximum sensor range is 5 meters
  // std::vector<float> ranges{2.8f, 5.0f, 5.0f, 5.0f, 5.0f}; // Maximum sensor range is 5 meters

  std::vector<float> angles{0.0f, -0.43633f, 0.0f, 0.43633f, 0.87266f};
  // std::vector<float> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};
  // std::vector<float> angles{0.785398f, -0.43633f, 0.0f, 0.43633f, 0.87266f};
  // std::vector<float> angles{0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

  // Current yaw + beam angle: -PI/2 (-1.570795) -0.87266 = 2.44345 (-55 degrees)
  // for (int i = 0; i < ranges.size(); ++i)
  for (int i = 0; i < 1; ++i) // One reading only
  {
    std::cout << "........ New laser ........" << std::endl;
    std::cout << "Distance: " << ranges[i] << ", Angle: " << angles[i] << std::endl;

    // Laser continuous distance
    std::vector<float> laser_grid = getLaserHit(robot_pose, ranges[i], angles[i]);
    // Laser final cell
    std::vector<int> final_grid_pos = getGridPosition(laser_grid[0], laser_grid[1], resolution);
    std::cout << final_grid_pos[0] << ", " << final_grid_pos[1] << std::endl;

    std::vector<float> cell_pos = getCellPosition({11, 12}, resolution);
    std::cout << cell_pos[0] << ", " << cell_pos[1] << std::endl;

    // robot_grid_pos[0] // X1 - robot_grid_pos[1] // Y1
    // final_grid_pos[0] // X2 - final_grid_pos[1] // Y2

    std::vector<int> cells_x, cells_y;
    std::pair<std::vector<int>, std::vector<int>> res_pair = Bresenham(robot_grid_pos[0], robot_grid_pos[1], final_grid_pos[0], final_grid_pos[1]);

    // Cells visited by this laser beam
    cells_x = res_pair.first;
    cells_y = res_pair.second;

    // Adding last hit cell to the set
    cells_x.push_back(final_grid_pos[0]); 
    cells_y.push_back(final_grid_pos[1]);

    std::cout << "Cells" << std::endl;
    for (int c = 0; c < cells_x.size(); ++c) // One reading only
    {
      std::cout << cells_x[c] << ", " << cells_y[c] << std::endl;
    }
    std::cout << "End of cells" << std::endl;

    std::cout << " ...---...---...---...---...---...--- " << std::endl;

    inverseMeasurement(grid_prob, grid_logs, cells_x, cells_y, robot_grid_pos, ranges[i], angles[i], resolution);

    std::cout << " ...---...---...---...---...---...--- " << std::endl;

    std::cout << "Probability: " << grid_prob[9][8] << std::endl;
    std::cout << "Logs: " << grid_logs[9][8] << std::endl;

    // std::cout << "Probability: " << grid_prob[11][13] << std::endl;
    // std::cout << "Probability: " << grid_prob[11][15] << std::endl;

    std::vector<float> initial_point(2), final_point(2);

    initial_point[0] = 0.0f;
    initial_point[1] = 0.0f;

    std::cout << " +++++++++++++ " << std::endl;
    // Visiting the cells this beam travels trough
    for (int j = 0; j < cells_x.size(); ++j)
    {
      if ((robot_grid_pos[0] == cells_x[j]) &&  (robot_grid_pos[1] == cells_y[j]))
      {
        // We do not need to calculate in the first cell 
        continue;
      }

      // Cells visualization
      std::cout << cells_x[j] << ", " << cells_y[j] << std::endl;

      // Converting the cell into distance
      float limit_x = cells_x[j] * resolution;
      float limit_y = cells_y[j] * resolution;

      std::vector<float> initial_x {limit_x, limit_x, limit_x + resolution, limit_x + resolution};
      std::vector<float> initial_y {limit_y, limit_y, limit_y + resolution, limit_y + resolution};

      std::vector<float> final_x {limit_x + resolution, limit_x, limit_x + resolution, limit_x};
      std::vector<float> final_y {limit_y, limit_y + resolution, limit_y, limit_y + resolution};

      float min_x = limit_x;
      float max_x = limit_x + resolution;
      float min_y = limit_y;
      float max_y = limit_y + resolution;

      if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] >= robot_grid_pos[1])
      {
        // X minor and Y greater. WRO final points
        initial_x[2] = limit_x - resolution;
        initial_x[3] = limit_x - resolution;

        final_x[0] = limit_x - resolution;
        final_x[2] = limit_x - resolution;

        min_x = limit_x - resolution;
        max_x = limit_x;
      }

      if (final_grid_pos[0] >= robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
      {
        // X greater and Y minor. WRO final points
        initial_y[2] = limit_y - resolution;
        initial_y[3] = limit_y - resolution;

        final_y[1] = limit_y - resolution;
        final_y[3] = limit_y - resolution;

        min_y = limit_y - resolution;
        max_y = limit_y;
      }

      if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
      {
        // X minor and Y minor. WRO final points
        initial_x[2] = limit_x - resolution;
        initial_x[3] = limit_x - resolution;
        initial_y[2] = limit_y - resolution;
        initial_y[3] = limit_y - resolution;

        final_x[0] = limit_x - resolution;
        final_x[2] = limit_x - resolution;
        final_y[1] = limit_y - resolution;
        final_y[3] = limit_y - resolution;

        min_x = limit_x - resolution;
        max_x = limit_x;
        min_y = limit_y - resolution;
        max_y = limit_y;
      }

      std::vector<float> inter_x, inter_y; 
      for (int k = 0; k < 4; ++k)
      {
        //std::cout << "Points: " << initial_x[k] << ", " << initial_y[k] << ", " << final_x[k] << ", " << final_y[k] << std::endl;
        std::vector<float> intersection = calculateIntersection(robot_pose, laser_grid, {initial_x[k], initial_y[k]}, {final_x[k], final_y[k]});
        if(intersection.size() != 0)
        {
          // If the laser and a cell intersects, we need to make sure it happens in the right bounds
          if ((abs(intersection[0]) >= abs(min_x - 0.001f)) &&
            (abs(intersection[0]) <= abs(max_x + 0.001f)) &&
            (abs(intersection[1]) >= abs(min_y - 0.001f)) &&
            (abs(intersection[1]) <= abs(max_y + 0.001f)))
          {
            /*
              Two points where the beam cuts the cell
              - A laser beam can cut the cell at least 1 time (Enter)
              - A laser beam can cut the cell at most 2 times (Enter an exit)
            */
            std::cout << "Interception at: " << intersection[0] << ", " << intersection[1] << std::endl;
            inter_x.push_back(intersection[0]);
            inter_y.push_back(intersection[1]);
          }
        }
      }
      
      // Enter (d1) and Exit (d2) distances
      std::vector<float> distances;
      for (int k = 0; k < inter_x.size(); ++k)
      {
        /*
          dist_point[0]: d1
          dist_point[1]: d2
        */  
        float dist_point = calculateDistance(robot_pose[0], robot_pose[1], inter_x[k], inter_y[k]);
        distances.push_back(dist_point);

        std::cout << "Distance: " << dist_point << std::endl;
      }

      /*
        For last cell I need to modify something
      */

      // Integral 1: 0 to d1 which is distance from robot pose to first point where the cell is cut
      // Integral 2: d1 which is distance from robot pose to first point where the cell is cut to 
      // d2 which is distance from robot pose to second point where the cell is cut
      // Integral 3: d2 which is distance from robot pose to second point where the cell is cut to z_max

      // I think this one is working well
      float prob_not = calculateProbability(0.0f, distances[0]); // 3.8 - Does not oberserve
      float prob_occ = calculateProbability(distances[0], distances[1]); // 3.9 - Occupied
      float prob_free = calculateProbability(distances[1], 5.0f); // 3.10 - Free
      std::cout << "Probabilities: " << prob_not << ", " << prob_occ << ", " << prob_free << std::endl;
      
    
      // The entropy of a grid cell C given a set of measurement outcomes z is given by 3.5
      // I need to calculate the entropy of a cell given a set of measurements

      // Inverse model
      // Probabiloityu
      // Entropy
      std::cout << " ++++++++++++++++++++++++ " << std::endl;
    }
  }
}


void LifelongSlamToolbox::inverseMeasurement(
  std::vector<std::vector<float>>& grid_prob, 
  std::vector<std::vector<float>>& grid_logs, 
  std::vector<int>& cells_x, 
  std::vector<int>& cells_y, 
  std::vector<int>& robot_grid_pos, 
  float range, float angle, float resolution)
{
  /*
    On this function we update a probability map based 
    We capture the cells that are hitted by the ray and the conditions
    And we update the probability as being occupied or free
  */

  float alpha = 1.0f;
  float max_r = 5.0f;

  for (int i = 0; i < cells_x.size(); ++i) // One reading only
  {
    std::cout << "Cells: " << cells_x[i] << ", " << cells_y[i] << std::endl;
    // This ditances will be calculated to the center of mass
    // mxi - pos_x //--// myi - pos_y
    float dx = (cells_x[i] - robot_grid_pos[0]) * resolution;
    float dy = (cells_y[i] - robot_grid_pos[1]) * resolution;

    float r = sqrt(pow(dx , 2) + pow(dy, 2));
    float phi = atan2(dy, dx) - robot_grid_pos[2];
    
    
    // I can assing only one value and calculate it once (Half the code)
  
    // Cell is occupied -- This condition might be different
    if ((range < max_r) && (abs(r - range) < (alpha / 2.0f)))
    {
      std::cout << "Occupied" << std::endl;
      // Update the probability
      updateCellProbability(grid_prob, 0.7f, cells_x[i], cells_y[i]);
      // Update the log-odds
      updateCellLogs(grid_prob, grid_logs, cells_x[i], cells_y[i], 0.0f);
      // Log-odds to probability
      float prob = logToProbability(grid_logs[cells_x[i]][cells_y[i]]);
    }
    // Cell is free
    else if (r <= range)
    {
      std::cout << "Free" << std::endl;
      // Update the probability
      updateCellProbability(grid_prob, 0.3f, cells_x[i], cells_y[i]);
      // Update the log-odds
      updateCellLogs(grid_prob, grid_logs, cells_x[i], cells_y[i], 0.0f);
      // Log-odds to probability
      float prob = logToProbability(grid_logs[cells_x[i]][cells_y[i]]);
    }

    std::cout << "Relative range: " << r << ", Angle: " << phi << std::endl;
  }
}


/*
  Next steps
  Done - log-odds will be acumulating the probability based on different observations - Bayes filter
  Transform the log-odds into probability 
  Calculate the entropy of a given cell
*/

void LifelongSlamToolbox::calculateCellEntropy()
{
  // Need for this one a matrix with all the entropies
  // Need to define where this calculation will take place
  ;
}


float LifelongSlamToolbox::logToProbability(float log)
{
  /*  
    To transform the Log-odds into probability
  */
  return (exp(log) / (1 + exp(log)));
}

void LifelongSlamToolbox::updateCellLogs(std::vector<std::vector<float>>& grid_prob, std::vector<std::vector<float>>& grid_logs, int cell_x, int cell_y, float initial_log)
{
  /*
    To update the log-odds matrix
  */
  grid_logs[cell_x][cell_y] = grid_logs[cell_x][cell_y] + calculateLogs(grid_prob[cell_x][cell_y]) - initial_log;
  // std::cout << "Cell Logs: " << grid_logs[cell_x][cell_y] << std::endl;
}

void LifelongSlamToolbox::updateCellProbability(std::vector<std::vector<float>>& grid_prob, float probability, int cell_x, int cell_y)
{
  /*
    To perform the probability update at the given cell
  */
  grid_prob[cell_x][cell_y] = probability;
  // std::cout << "Cell probability: " << grid_prob[cell_x][cell_y] << std::endl;
}

float LifelongSlamToolbox::calculateLogs(float probability)
{
  /*
    To calculate the log-odds
  */
  // This is natural algorithm
  return log(probability / (1 - probability));
}

std::vector<float> LifelongSlamToolbox::getCellPosition(std::vector<int> grid_cell, float resolution)
{
  /*
    To get the current cell
  */
  float x = (grid_cell[0] * resolution) + (resolution / 2);
  float y = (grid_cell[1] * resolution) + (resolution / 2);

  return {x, y};
}

std::vector<float> LifelongSlamToolbox::calculateIntersection(
  std::vector<float> laser_start, std::vector<float> laser_end,
  std::vector<float> cell_start, std::vector<float> cell_end)
{
  /*
    Initial point laser beam: laser_start
    Final point laser beam: laser_end
    Initial point cell: cell_start
    Final point cell: cell_end
  */
  float x1 = laser_start[0];
  float x2 = laser_end[0];
  float x3 = cell_start[0];
  float x4 = cell_end[0];

  float y1 = laser_start[1];
  float y2 = laser_end[1];
  float y3 = cell_start[1];
  float y4 = cell_end[1];

  float den = ((x2-x1)*(y4-y3) - (x4-x3)*(y2-y1));
  if (den == 0.0f)
  {
      // Parallel lines or not intersection at all
      return {};
  }
  else
  {
      float x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2-x1)) / den;
      float y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2-y1)) / den;
      // std::cout << x<< ", " << y << std::endl;
      return {x, y};
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

float LifelongSlamToolbox::calculateProbability(float range_1, float range_2)
{
  /*
    Calculates the probability of a cell being observed by a given measurement
    range_1: lower limit, range_2: upper limit
  */
  float max_range = 5.0f;
  float lambda = 0.285f;
  float nu = 1.0f / lambda;

  range_2 = (range_2 > max_range) ? max_range : range_2;

  // https://www.wolframalpha.com/input/?i2d=true&i=Integrate%5Bn*c*Power%5Be%2C-c*x%5D%2C%7Bx%2Ca%2Cb%7D%5D
  return nu * (exp(-lambda*range_1) - exp(-lambda*range_2));
}

std::vector<int> LifelongSlamToolbox::getGridPosition(float x, float y, float resolution)
{
  /*
    Maps the distance into grid coordinates
  */
  int x_cell = floor((1 / resolution) * x);
  int y_cell = floor((1 / resolution) * y);

  return {x_cell, y_cell};
}

std::vector<float> LifelongSlamToolbox::getLaserHit(std::vector<float> const& robot_pose, float range, float angle)
{
  /*
    Returns the distance where the laser beam hits something
    Rigid Body Trasnformation from the global to the sensor frame
  */
  float x_tf = (range * cos(robot_pose[2] + angle)) + robot_pose[0];
  float y_tf = (range * sin(robot_pose[2] + angle)) + robot_pose[1];

  return {x_tf, y_tf};
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

#ifndef M3RSM__H
#define M3RSM__H

#include <stdint.h>
#include <map>
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <chrono>
#include <utility>
#include <string>

#include "Eigen/Core"
#include "rclcpp/rclcpp.hpp"
#include "karto_sdk/Karto.h"
#include "karto_sdk/Mapper.h"

namespace {
  typedef std::vector<std::vector<uint8_t>> DataMatrix;
}

namespace karto
{

class M3RSM_scan;

// TODO: Implement this class.
class LookupTable {

  public:
    // Builds highest resolution lookup table from base scans
    LookupTable(const std::vector<LocalizedRangeScan>& base_scans, int decimation_factor, int kernel_width);
    // Copy constructor
    LookupTable(const LookupTable& table);
    // Builds lower resolution lookup table from higher resolution lookup table
    static LookupTable GetLowerResTable(const LookupTable& higher_res_table);
    double computeCost(const LocalizedRangeScan& scan);
    bool isMaxResolution();
    int getTableSize();

  private:
    DataMatrix data_;
    bool is_max_res_;
    int decimation_factor_;
    int kernel_width_;

};  // class LookupTable

class MockLookupTable {

  public:
    MockLookupTable(const DataMatrix& data, int decimation_factor, int kernel_width);
    // Builds lower resolution lookup table from higher resolution lookup table
    MockLookupTable(const MockLookupTable& table);
    static MockLookupTable GetLowerResTable(const MockLookupTable& higher_res_table);
    double computeCost(const std::vector<std::pair<int, int>>& occupied_cells);
    bool isMaxResolution();
    int getTableSize();
    DataMatrix data_;
    bool is_max_res_;
    int decimation_factor_;
    int kernel_width_;
};  // class MockLookupTable

class LookupTableManager {

  public:
    // TODO: Implement this method when integrating with slam_toolbox codebase.
    // LookupTableManager(const std::vector<LocalizedRangeScan>& base_scans);
    LookupTableManager(const DataMatrix& full_res_matrix);
    MockLookupTable getTable(int resolution);
    int getAmountOfResolutions();

  private:
    // Lookup tables ordered by resolution.
    // Each vector contains a vector of same-resolution lookup tables.
    // TODO: Use this attribute when integrating with slam_toolbox.
    //std::vector<std::vector<LookupTable>> tables_;
    std::vector<MockLookupTable> tables_;
    int n_resolutions_;

};  // class LookupTableManager

class EvaluationTask {
  public:
    EvaluationTask(const LookupTable& lookup_table);
    double computeCost();
    std::vector<EvaluationTask> generateChildrenTasks();
    bool isMaxResolution();
    Pose2 getPose();

  private:
    std::shared_ptr<LookupTable> lookup_table_;
};  // class EvaluationTask

class SearchHeap {

  public:
    SearchHeap(const LocalizedRangeScan& received_scan, const LookupTableManager& lookup_table_manager);
    Pose2 getEstimatedPose(); 

  private:
    std::shared_ptr<LocalizedRangeScan> received_scan_;
    std::shared_ptr<LookupTableManager> lookup_table_manager_;
    std::vector<EvaluationTask> tasks;

};  // class SearchHeap


class M3RScanMatcher {

public:
  M3RScanMatcher();
  Pose2 MatchScan(LocalizedRangeScan received_scan, std::vector<LocalizedRangeScan> potential_scans);

  private:
};  // M3RScanMatcher

} // namespace karto 

#endif  // M3RSM__H

#ifndef SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_

#include <algorithm>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <unordered_map>
#include "lib/karto_sdk/include/karto_sdk/Karto.h"

// #include "utils.hpp"

#include "slam_toolbox/slam_toolbox_common.hpp"
#include "slam_toolbox/experimental/slam_toolbox_lifelong.hpp"

class InformationEstimates
{
    typedef std::tuple<int, int, int> map_tuple;
    typedef std::pair<map_tuple, kt_double> map_pair;

public:
    InformationEstimates(); // Empty contructor for now
    virtual ~InformationEstimates() {}

public:
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

public:
    float calculateMutualInformation(karto::PointVectorDouble const& laser_readings, karto::Pose2 const& karto_pose);

private:
    // Grid operations
    void initializeGrids();
    void updateCellLimits(
        std::vector<kt_double>& initial_x, std::vector<kt_double>& initial_y, std::vector<kt_double>& final_x, std::vector<kt_double>& final_y, 
        kt_double limit_x, kt_double limit_y, std::vector<kt_double>& cell_limits, karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos);

    // Grid and position information
    std::pair<std::vector<int>, std::vector<int>> rayCasting(karto::Vector2<int> const& initial_pt, karto::Vector2<int> const& final_pt);
    karto::Vector2<int> getGridPosition(karto::Vector2<kt_double> const& pose);
    std::vector<kt_double> calculateCellIntersectionPoints(karto::Vector2<kt_double> const & laser_start, karto::Vector2<kt_double> const & laser_end, std::vector<kt_double> cell_start, std::vector<kt_double> cell_end);
    std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
        karto::Vector2<kt_double> const & laser_start, karto::Vector2<kt_double> const & laser_end, 
        karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos,
        kt_double limit_x, kt_double limit_y);
    int signum(int num);
    void clearVisitedCells();

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);

    // Mutual information 
    kt_double calculateInformationContent(kt_double prob);
    kt_double calculateLogOddsFromProbability(kt_double probability);
    kt_double calculateMapMutualInformation();
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);
    void recoverProbability();
    void updateCellMutualInformation(kt_double mut_inf, std::vector<int> cell);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<kt_double>& measurements, std::vector<int> cell);
    std::unordered_map<map_tuple, kt_double, HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm);
    std::vector<std::vector<kt_double>> retreiveCellProbabilities(std::vector<int> cell);


private:
    // Data structures 
    std::unordered_map<map_tuple, kt_double, HashTuple> m_map_out;
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;
    std::vector<std::vector<kt_double>> m_mutual_grid;
    std::vector<std::vector<bool>> m_visited_grid;

    kt_double m_map_dist;
    kt_double m_cell_resol;
    int m_num_cells;

    const kt_double l_free = log(0.3 / (1.0 - 0.3));
    const kt_double l_occ = log(0.7 / (1.0 - 0.7));
    const kt_double l_o = log(0.5 / (1.0 - 0.5));


    kt_double m_max_sensor_range = 5.0;
    kt_double m_obs_lambda = 0.35; 
    kt_double m_obs_nu = 0.28;

public:
    // Setters
    void setMaxSensorRange(kt_double const sensor_range);
    void setObservationLambda(kt_double const lambda);
    void setObservationNu(kt_double const nu);
    void setCellResolution(kt_double const resolution);
    void setMapDistance(kt_double const distance);
};

#endif

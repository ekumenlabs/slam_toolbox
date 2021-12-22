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

class InformationEstimates
{
    typedef std::tuple<int, int, int> map_tuple;
    typedef std::pair<map_tuple, double> map_pair;

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

private:
    // Testing
    void scannerTest();

    // Grid operations
    void initializeGrids();
    void updateCellLimits(
        std::vector<double>& initial_x, std::vector<double>& initial_y, std::vector<double>& final_x, std::vector<double>& final_y, 
        double limit_x, double limit_y, std::vector<double>& cell_limits, std::vector<int> const& robot_grid_pos, std::vector<int> const& final_grid_pos);

    // Grid and position information
    std::pair<std::vector<int>, std::vector<int>> rayCasting(int x_1, int y_1, int x_2, int y_2);
    std::vector<int> getGridPosition(double x, double y);
    std::vector<double> laserHitDistance(std::vector<double> const& robot_pose, double range, double angle);
    std::vector<double> calculateCellIntersectionPoints(std::vector<double> const& laser_start, std::vector<double> const& laser_end, std::vector<double> cell_start, std::vector<double> cell_end);
    int signum(int num);

    std::pair<std::vector<double>, std::vector<double>> computeLineBoxIntersection(
        std::vector<double> const& laser_start, std::vector<double> const& laser_end, 
        std::vector<int> const& robot_grid_pos, std::vector<int> const& final_grid_pos, 
        double limit_x, double limit_y);

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    double calculateScanMassProbabilityBetween(double range_1, double range_2);
    double euclideanDistance(double x_1, double y_1, double x_2, double y_2);

    // Mutual information 
    double calculateInformationContent(double prob);
    double calculateLogOddsFromProbability(double probability);
    double calculateMapMutualInformation();
    double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    double calculateProbabilityFromLogOdds(double log);
    void recoverProbability();
    void updateCellMutualInformation(double mut_inf, std::vector<int> cell);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<double>& measurements, std::vector<int> cell);
    std::unordered_map<map_tuple, double, HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<double>>& meas_outcm);
    std::vector<std::vector<double>> retreiveCellProbabilities(std::vector<int> cell);

    void clearVisitedCells();

private:
    // Data structures 
    std::unordered_map<map_tuple, double, HashTuple> m_map_out;
    std::map<std::vector<int>, std::vector<std::vector<double>>> m_cell_probabilities;
    std::vector<std::vector<double>> m_mutual_grid;
    std::vector<std::vector<int>> m_grid;

    double m_map_dist;
    double m_cell_resol;
    int m_num_cells;

    const double l_free = log(0.3 / (1.0 - 0.3));
    const double l_occ = log(0.7 / (1.0 - 0.7));
    const double l_o = log(0.5 / (1.0 - 0.5));

    // Robot information - Defining values just for testing
    // std::vector<std::vector<double>> robot_poses {{5.6f, 6.0f, M_PI/2}, {7.75f, 9.22f, -M_PI/2}};
    // std::vector<std::vector<double>> robot_poses {{7.75f, 9.22f, -M_PI/2}, {3.5f, 9.0f, 0.0f}};
    std::vector<std::vector<double>> robot_poses {{5.6f, 6.0f, M_PI/2}, {3.5f, 9.0f, 0.0f}};
    std::vector<std::vector<double>> laser_ranges {{1.65f, 5.0f, 5.0f, 5.0f, 5.0f}, {5.0f, 5.0f, 4.0f, 5.0f, 5.0f}};
    // std::vector<std::vector<double>> laser_ranges {{5.0f, 5.0f, 5.0f, 5.0f, 5.0f}, {5.0f, 5.0f, 4.0f, 5.0f, 5.0f}};
    std::vector<double> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

    std::vector<std::vector<bool>> visited_cells;

    double m_max_sensor_range = 5.0;
    double m_obs_lambda = 0.35; 
    double m_obs_nu = 0.28;

    void setMaxSensorRange(double const sensor_range);
    void setObservationLambda(double const lambda);
    void setObservationNu(double const nu);
    void setCellResolution(double const resolution);
    void setMapDistance(double const distance);

    karto::Vector2<kt_double> point;
    // This pose will save actual robot poses
    karto::Pose2 pose;

};

#endif

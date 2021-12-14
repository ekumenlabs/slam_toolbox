#ifndef SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_

#include <algorithm>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <unordered_map>

class TeorethicInformation
{
    typedef std::tuple<int, int, int> map_tuple;
    typedef std::pair<map_tuple, float> map_pair;

public:
    TeorethicInformation(); // Empty contructor for now
    ~TeorethicInformation() {}

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
        std::vector<float> & initial_x, std::vector<float> & initial_y, std::vector<float> & final_x, std::vector<float> & final_y, 
        float & limit_x, float & limit_y, std::vector<float> & cell_limits, std::vector<int> & robot_grid_pos, std::vector<int> & final_grid_pos);

    // Grid and position information
    std::pair<std::vector<int>, std::vector<int>> Bresenham(int x_1, int y_1, int x_2, int y_2);
    std::vector<int> getGridPosition(float x, float y);
    std::vector<float> laserHitDistance(std::vector<float> const& robot_pose, float range, float angle);
    std::vector<float> calculateCellIntersectionPoints(std::vector<float> & laser_start, std::vector<float> & laser_end, std::vector<float> & cell_start, std::vector<float> & cell_end);
    int getSign(int n_1, int n_2);

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    float probabilityFromObservation(float range_1, float range_2);
    float euclideanDistance(float x_1, float y_1, float x_2, float y_2);

    // Mutual information 
    float calculateEntropy(float probability);
    float calculateLogs(float probability);
    float calculateMapMutualInformation();
    float measurementOutcomeEntropy(map_tuple const& meas_outcome);
    float probabilityFromLogs(float log);
    void recoverProbability();
    void updateCellMutualInformation(float mut_inf);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<float>& measurements);
    void computeProbabilities(std::vector<std::vector<float>>& meas_outcm);
    std::vector<std::vector<float>> retreiveMeasurementOutcomes();
    std::vector<int> unhashIndex(int hash);

private:
    // Data structures 
    std::unordered_map<map_tuple, float, HashTuple> m_map_out;
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
};

#endif
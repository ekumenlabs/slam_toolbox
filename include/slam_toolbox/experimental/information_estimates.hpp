#ifndef SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_

#include "utils.hpp"

class InformationEstimates
{
    typedef std::tuple<int, int, int> map_tuple;

public:
    InformationEstimates(); // Empty contructor for now
    virtual ~InformationEstimates() {}

private:
    // Mutual information 
    kt_double calculateInformationContent(kt_double prob);
    kt_double calculateMapMutualInformation();
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);
    void updateCellMutualInformation(kt_double mut_inf, std::vector<int> cell);
    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<kt_double>& measurements, std::vector<int> cell);
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm);
    std::vector<std::vector<kt_double>> retreiveCellProbabilities(std::vector<int> cell);
    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);
    void updateLaserMutualInformation();

private:
    // Data structures 
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

    kt_double m_map_mutual_info;
    kt_double m_laser_mutual_info;

public:
    // Main function
    float calculateMutualInformation(karto::PointVectorDouble const& laser_readings, karto::Pose2 const& karto_pose);

    // Setters
    void setMaxSensorRange(kt_double const sensor_range);
    void setObservationLambda(kt_double const lambda);
    void setObservationNu(kt_double const nu);
    void setCellResolution(kt_double const resolution);
    void setMapDistance(kt_double const distance);
    

    // Getters
    kt_double getMapMutualInformation();
    kt_double getLaserMutualInformation();

};

#endif

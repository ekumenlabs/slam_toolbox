#ifndef SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__THEORETIC_INFORMATION_HPP_

#include "utils.hpp"

class InformationEstimates
{
    typedef std::tuple<int, int, int> map_tuple;

public:
    InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu);
    InformationEstimates();
    virtual ~InformationEstimates() {}

public:
    // Main function
    std::vector<kt_double> findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans);

private:
    // Mutual information
    kt_double calculateInformationContent(kt_double prob);
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);
    kt_double mutualInformationFromScans(std::vector<karto::LocalizedRangeScan*> const& range_scans, bool ignore_scan=false, int scan_idx=0);
    void updateCellMutualInformation(kt_double mut_inf, karto::Vector2<int> const & cell);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<kt_double>& measurements, karto::Vector2<int> const & cell);
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm);
    void insertMeasurementOutcome(map_tuple tuple, kt_double probability, std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple>& map);
    std::vector<std::vector<kt_double>> retrieveCellProbabilities(karto::Vector2<int> const & cell);

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);

private:
    // Data structures
    std::map<karto::Vector2<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;

    kt_double m_max_sensor_range;
    kt_double m_cell_resol;
    kt_double m_obs_lambda;
    kt_double m_obs_nu;

    kt_double m_low_x;
    kt_double m_low_y;
    kt_double m_high_x;
    kt_double m_high_y;

    // Map grids
    Eigen::MatrixXd m_mutual_grid;
    Eigen::MatrixXi m_visited_grid;
};

#endif

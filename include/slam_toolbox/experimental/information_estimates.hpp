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
    std::tuple<int, kt_double> findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans);

private:
    // Mutual information
    kt_double calculateInformationContent(kt_double prob);
    kt_double calculateMapMutualInformation();
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);
    void updateCellMutualInformation(kt_double mut_inf, karto::Vector2<int> const & cell);

    // Measurement outcomes probabilities
    void appendCellProbabilities(std::vector<kt_double>& measurements, karto::Vector2<int> const & cell);
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm);
    std::vector<std::vector<kt_double>> retrieveCellProbabilities(karto::Vector2<int> const & cell);

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);
    kt_double calculateLaserMutualInformation();

private:
    // Data structures
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;

    const kt_double l_free = log(0.3 / (1.0 - 0.3));
    const kt_double l_occ = log(0.7 / (1.0 - 0.7));
    const kt_double l_o = log(0.5 / (1.0 - 0.5));

    kt_double m_max_sensor_range;
    kt_double m_cell_resol;
    kt_double m_obs_lambda;
    kt_double m_obs_nu;
    kt_double m_map_dist;
    int m_num_cells;

    // Map grids
    Eigen::MatrixXd m_info_grid;
    Eigen::MatrixXd m_mutual_grid;
    Eigen::MatrixXi m_visited_grid;
};

#endif

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
    std::tuple<int, kt_double> calculateMutualInformation(std::vector<karto::LocalizedRangeScan> const& range_scans);

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
    std::vector<std::vector<kt_double>> retrieveCellProbabilities(std::vector<int> cell);
    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);
    kt_double calculateLaserMutualInformation(kt_double const & map_info, kt_double const & curr_info);

private:
    // Data structures 
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;
    std::vector<std::vector<kt_double>> m_mutual_grid;
    std::vector<std::vector<bool>> m_visited_grid;

    const kt_double l_free = log(0.3 / (1.0 - 0.3));
    const kt_double l_occ = log(0.7 / (1.0 - 0.7));
    const kt_double l_o = log(0.5 / (1.0 - 0.5));

    kt_double m_max_sensor_range;
    kt_double m_cell_resol;
    kt_double m_obs_lambda; 
    kt_double m_obs_nu;

    kt_double m_map_dist;
    int m_num_cells;

    kt_double m_curr_mut_info;
};

#endif

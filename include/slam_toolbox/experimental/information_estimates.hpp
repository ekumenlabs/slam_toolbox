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
    // std::vector<kt_double> findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans);
    std::vector<kt_double> findMutualInfo(std::vector<karto::LocalizedRangeScan*> const& range_scans);

private:

    void calculateAndAppendCellProbabilities(
        std::vector<karto::Vector2<int>> & visited_cells,
        std::vector<kt_double> const & distances,
        karto::Vector2<int> const & cell
    );

    std::optional<std::vector<kt_double>> calculateBeamAndCellIntersections(
        utils::Segment2<kt_double> const & beam_segment,
        karto::Vector2<int> const & cell
    );

    int findClosestLaserIndexToCell(
        kt_bool & skip_cell_eval,
        kt_double const & angle_to_cell,
        kt_double const & scan_pose_heading,
        karto::LaserRangeFinder *laser_range_finder
    );

    std::optional<kt_double> adjustBeamReadingDistance(
        kt_double const & beam_distance,
        kt_double const & distance_to_cell,
        karto::LaserRangeFinder *laser_range_finder
    );

    karto::Vector2<int> getLaserBeamCell(
        kt_double const & angle_to_cell,
        kt_double const & reading_distance
    );

    std::vector<karto::Vector2<int>> getScanGroupVisitedCells(
        std::vector<karto::LocalizedRangeScan *> const & range_scans,
        karto::LaserRangeFinder *laser_range_finder,
        int const & scan_to_skip
    );

    void calculateScanGroupMutualInformation(
        std::vector<karto::Vector2<int>> const & visited_cells,
        std::vector<kt_double> & scans_mutual_information
    );

    std::optional<int> InformationEstimates::findClosestLaserIndexToCell(
        kt_double const & angle_to_cell,
        kt_double const & scan_pose_heading,
        karto::LaserRangeFinder *laser_range_finder
    );

    void calculateCellProbabilities(
        std::vector<karto::LocalizedRangeScan *> const & range_scans,
        std::vector<karto::Vector2<int>> & visited_cells,
        karto::Vector2<int> const & cell,
        karto::LaserRangeFinder *laser_range_finder,
        int const & scan_to_skip
    );

    void resizeGridFromScans(
        std::vector<karto::LocalizedRangeScan *> const & range_scans
    );

    std::vector<kt_double> getScanGroupMutualInformation(
        std::vector<karto::LocalizedRangeScan *> const & range_scans
    );

    // Mutual information
    kt_double calculateInformationContent(kt_double prob);
    kt_double measurementOutcomeEntropy(map_tuple const& meas_outcome);
    kt_double calculateProbabilityFromLogOdds(kt_double log);

    // Measurement outcomes probabilities
    void appendCellProbabilities(
        std::vector<kt_double>& measurements,
        karto::Vector2<int> const & cell
    );

    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> computeMeasurementOutcomesHistogram(
        std::vector<std::vector<kt_double>>& meas_outcm
    );

    void insertMeasurementOutcome(
        map_tuple tuple,
        kt_double probability,
        std::unordered_map<map_tuple,kt_double, utils::tuple_hash::HashTuple>& map
    );

    // Measurements calculations <P(free), P(Occ), P(Unk)>
    kt_double calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2);

private:
    // Data structures
    std::map<karto::Vector2<int>, std::vector<std::vector<kt_double>>> m_cell_probabilities;

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
    Eigen::MatrixXd m_mutual_grid;
    Eigen::MatrixXi m_visited_grid;

    kt_double m_upper_limit_x;
    kt_double m_upper_limit_y;
    kt_double m_lower_limit_x;
    kt_double m_lower_limit_y;

    kt_double m_map_dim_x;
    kt_double m_map_dim_y;
};

#endif

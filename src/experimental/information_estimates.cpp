#include <math.h>
#include <cmath>
#include "slam_toolbox/experimental/information_estimates.hpp"

#include <iostream>

InformationEstimates::InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu)
    : m_max_sensor_range{sensor_range},
      m_cell_resol{resolution},
      m_obs_lambda{lambda},
      m_obs_nu{nu}
{
}

InformationEstimates::InformationEstimates()
    : m_max_sensor_range{25.0},
      m_cell_resol{0.1},
      m_obs_lambda{0.35},
      m_obs_nu{0.28}
{
}

/*************************************************************************/
//
void InformationEstimates::resizeGridFromScans(std::vector<karto::LocalizedRangeScan *> const & range_scans)
{
    m_lower_limit_x = range_scans[0]->GetCorrectedPose().GetX();
    m_lower_limit_y = range_scans[0]->GetCorrectedPose().GetY();
    m_upper_limit_x = range_scans[0]->GetCorrectedPose().GetX();
    m_upper_limit_y = range_scans[0]->GetCorrectedPose().GetY();

    for (const auto &scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();
        kt_double minimumAngle = scan->GetLaserRangeFinder()->GetMinimumAngle();
        kt_double angularResolution = scan->GetLaserRangeFinder()->GetAngularResolution();

        // Finding the lower limit pose
        m_lower_limit_x = std::min(pose.GetX(), m_lower_limit_x);
        m_lower_limit_y = std::min(pose.GetY(), m_lower_limit_y);

        // Finding the higher limit pose
        m_upper_limit_x = std::max(pose.GetX(), m_upper_limit_x);
        m_upper_limit_y = std::max(pose.GetY(), m_upper_limit_y);
    }

    // Map dimensions
    m_map_dim_x = std::fabs(m_upper_limit_x - m_lower_limit_x) + (2 * m_max_sensor_range);
    m_map_dim_y = std::fabs(m_upper_limit_y - m_lower_limit_y) + (2 * m_max_sensor_range);

    // Get the number of cells
    int n_cells_x = static_cast<int>(m_map_dim_x / m_cell_resol);
    int n_cells_y = static_cast<int>(m_map_dim_y / m_cell_resol);

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);
}
//
/*************************************************************************/


/*************************************************************************/
//
std::optional<int> InformationEstimates::findClosestLaserIndexToCell(
    kt_double const & angle_to_cell,
    kt_double const & scan_pose_heading,
    karto::LaserRangeFinder *laser_range_finder)
{
    const kt_double angle_between = angle_to_cell - scan_pose_heading;
    const kt_double angle_diff = atan2(sin(angle_between), cos(angle_between));

    if (angle_diff < laser_range_finder->GetMinimumAngle() || angle_diff >  laser_range_finder->GetMaximumAngle())
    {
        // We are outside the laser FOV, but we cannot cancel it
        // because we will not evaluate other cells
        return {};
    }

    const kt_double laser_middle_point = ( laser_range_finder->GetMaximumAngle() - laser_range_finder->GetMinimumAngle()) / 2.0;
    return std::optional { (laser_middle_point + angle_between) / laser_range_finder->GetAngularResolution() };
}
//
/*************************************************************************/


/*************************************************************************/
//
std::optional<std::vector<kt_double>> InformationEstimates::calculateBeamAndCellIntersections(
    utils::Segment2<kt_double> const & beam_segment,
    karto::Vector2<int> const & cell)
{
    karto::Vector2<kt_double> current_cell_position{ cell.GetX() * m_cell_resol, cell.GetY() * m_cell_resol };

    std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections =
        utils::grid_operations::computeLineBoxIntersection(
            beam_segment,
            current_cell_position,
            m_cell_resol
        );

    if (intersections.first.size() == 0)
    {
        return {};
    }

    // Enter (d1) and Exit (d2) distances
    std::vector<kt_double> distances;
    for (int k = 0; k < intersections.first.size(); ++k)
    {
        // From robot position to intersection points
        karto::Vector2<kt_double> intersection{
            intersections.first[k],
            intersections.second[k]
        };
        kt_double distance = beam_segment.start.Distance(intersection);
        distances.push_back(distance);
    }
    std::sort(distances.begin(), distances.end());

    return std::optional{ distances };
}
//
/*************************************************************************/

/*************************************************************************/
//
std::optional<kt_double> InformationEstimates::adjustBeamReadingDistance(
    kt_double const & beam_distance,
    kt_double const & distance_to_cell,
    karto::LaserRangeFinder *laser_range_finder)
{
    if ((beam_distance > laser_range_finder->GetRangeThreshold()) || (beam_distance < distance_to_cell))
    {
        return {};
    }
    return std::optional { distance_to_cell };
}
//
/*************************************************************************/


/*************************************************************************/
//
void InformationEstimates::calculateAndAppendCellProbabilities(
    std::vector<karto::Vector2<int>> & visited_cells,
    std::vector<kt_double> const & distances,
    karto::Vector2<int> const & cell)
{
    // Measurement outcomes vector {Pfree, Pocc, Pun}
    std::vector<kt_double> probabilities{
        calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
        calculateScanMassProbabilityBetween(distances[0], distances[1]),
        calculateScanMassProbabilityBetween(0.0f, distances[0])
    };

    visited_cells.push_back( {cell.GetX(), cell.GetY()});

    // Appending new measurement outcomes for the current cell
    appendCellProbabilities(probabilities, cell);
}
//
/*************************************************************************/


/*************************************************************************/
//
void InformationEstimates::calculateCellProbabilities(
    std::vector<karto::LocalizedRangeScan *> const & range_scans,
    std::vector<karto::Vector2<int>> & visited_cells,
    karto::Vector2<int> const & cell,
    karto::LaserRangeFinder *laser_range_finder,
    int const & scan_to_skip)
{
    for (int s = 0; s < range_scans.size(); ++s)
    {
        // Skip the given scan, because we are evaluating the mutual information
        // of the group, excluding it.
        if (scan_to_skip == s)
        {
            continue;
        }

        // Scan pose
        karto::Pose2 scan_pose = range_scans[s]->GetCorrectedPose();
        karto::Pose2 grid_scan_pose {
            scan_pose.GetX() - (m_lower_limit_x - m_max_sensor_range),
            scan_pose.GetY() - (m_lower_limit_y - m_max_sensor_range),
            scan_pose.GetHeading()
        };
        karto::Vector2<kt_double> cell_center{cell.GetX() * m_cell_resol + (m_cell_resol / 2), cell.GetY() * m_cell_resol + (m_cell_resol / 2)};
        kt_double angle_to_cell = atan2(cell_center.GetY() - grid_scan_pose.GetY(), cell_center.GetX() - grid_scan_pose.GetX());
        kt_double distance_to_cell = sqrt(pow(grid_scan_pose.GetX() - cell_center.GetX(), 2) + pow(grid_scan_pose.GetY() - cell_center.GetY(), 2));
        if(distance_to_cell > m_max_sensor_range)
        {
            // Cell is not being seen
            continue;
        }

        std::optional<int> laser_idx = findClosestLaserIndexToCell(
            angle_to_cell,
            scan_pose.GetHeading(),
            laser_range_finder
        );

        if (!laser_idx.has_value())
        {
            continue;
        }

        kt_double beam_read_dist = range_scans[s]->GetRangeReadings()[*laser_idx];
        std::optional<kt_double> beam_distance = adjustBeamReadingDistance(beam_read_dist, distance_to_cell, laser_range_finder);

        if (!beam_distance.has_value())
        {
            continue;
        }

        karto::Vector2<kt_double> beam_point {
            grid_scan_pose.GetX() + (*beam_distance * cos(angle_to_cell)),
            grid_scan_pose.GetY() + (*beam_distance * sin(angle_to_cell))
        };

        utils::Segment2<kt_double> beam_segment {
            grid_scan_pose.GetPosition(),
            beam_point
        };

        std::optional<std::vector<kt_double>> distances = calculateBeamAndCellIntersections(
            beam_segment,
            cell
        );

        if (!distances.has_value())
        {
            continue;
        }

        calculateAndAppendCellProbabilities(
            visited_cells,
            *distances,
            cell
        );
    }
}
//
/*************************************************************************/


/*************************************************************************/
//
std::vector<kt_double> InformationEstimates::getScanGroupMutualInformation(
    std::vector<karto::LocalizedRangeScan *> const & range_scans
)
{
    karto::LaserRangeFinder *laser_range_finder = range_scans[0]->GetLaserRangeFinder();
    kt_double range_threshold = laser_range_finder->GetRangeThreshold();
    kt_double angle_increment = laser_range_finder->GetAngularResolution();
    kt_double min_angle = laser_range_finder->GetMinimumAngle();
    kt_double max_angle = laser_range_finder->GetMaximumAngle();
    kt_double total_range = max_angle - min_angle;

    std::vector<kt_double> scans_mutual_information;

    for (int n = 0; n < range_scans.size(); ++n)
    {
        m_mutual_grid.setZero();
        m_visited_grid.setZero();

        std::vector<karto::Vector2<int>> visited_cells = getScanGroupVisitedCells(range_scans, laser_range_finder, n);
        calculateScanGroupMutualInformation(visited_cells, scans_mutual_information);
    }
    return scans_mutual_information;
}
//
/*************************************************************************/


/*************************************************************************/
//
void InformationEstimates::calculateScanGroupMutualInformation(
    std::vector<karto::Vector2<int>> const & visited_cells,
    std::vector<kt_double> & scans_mutual_information
)
{
    for (auto &cell : visited_cells)
    {
        std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob =
            computeMeasurementOutcomesHistogram(m_cell_probabilities.at(cell));
        kt_double cell_mutual_inf = 0.0f;
        for (auto &pair : meas_out_prob)
        {
            cell_mutual_inf += pair.second * measurementOutcomeEntropy(pair.first);
        }

        // Mutual information of cell x, y given a set of measurements
        m_mutual_grid(cell.GetX(), cell.GetY()) = 1.0 - cell_mutual_inf;
    }
    scans_mutual_information.push_back(m_mutual_grid.sum());
    m_cell_probabilities.clear();
}
//
/*************************************************************************/


/*************************************************************************/
//
std::vector<karto::Vector2<int>> InformationEstimates::getScanGroupVisitedCells(
    std::vector<karto::LocalizedRangeScan *> const & range_scans,
    karto::LaserRangeFinder *laser_range_finder,
    int const & scan_to_skip
)
{
    std::vector<karto::Vector2<int>> visited_cells;

    int x_upper_limit = floor(m_map_dim_x / m_cell_resol);
    int y_upper_limit = floor(m_map_dim_y / m_cell_resol);

    for (int i = 0; i < x_upper_limit; ++i)
    {
        for (int j = 0; j < y_upper_limit; ++j)
        {
            calculateCellProbabilities(range_scans, visited_cells, {i, j}, laser_range_finder, scan_to_skip);
        }
    }
    return visited_cells;
}
//
/*************************************************************************/


/*************************************************************************/
std::vector<kt_double> InformationEstimates::findMutualInfo(std::vector<karto::LocalizedRangeScan *> const &range_scans)
{
    std::vector<kt_double> result_vector;
    resizeGridFromScans(range_scans);

    std::vector<kt_double> mutual_information_vct = getScanGroupMutualInformation(range_scans);
    return mutual_information_vct;
}

void InformationEstimates::appendCellProbabilities(std::vector<kt_double> &measurements, karto::Vector2<int> const &cell)
{
    /**
     * Append measured porbabilities for the given cell
     * Arguments:
     * measurements [std::vector<kt_double>]: Vector of LocalizedRangeScan pointers
     * cell [karto::Vector2<int>]: Cell for appending the data
     * Return:
     * Void
     */
    std::map<karto::Vector2<int>, std::vector<std::vector<kt_double>>>::iterator it_cell;
    it_cell = m_cell_probabilities.find(cell);

    if (it_cell == m_cell_probabilities.end())
    {
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<karto::Vector2<int>, std::vector<std::vector<kt_double>>>(
            cell, {measurements})
        );
        m_visited_grid(cell.GetX(), cell.GetY()) = 1;
    }
    else
    {
        if (m_visited_grid(cell.GetX(), cell.GetY()) == 1) // It must be the longitude
        {
            int idx = it_cell->second.size() - 1;
            if (measurements[2] < it_cell->second[idx][2])
            {
                // Replacing
                it_cell->second[idx][0] = measurements[0];
                it_cell->second[idx][1] = measurements[1];
                it_cell->second[idx][2] = measurements[2];
            }
        }
        else
        {
            it_cell->second.push_back(measurements);
            m_visited_grid(cell.GetX(), cell.GetY()) = 1;
        }
    }
}

kt_double InformationEstimates::calculateInformationContent(kt_double prob)
{
    /**
     * Calculate the information content or self-information based on the probability of cell being occupied
     * Arguments:
     * prob [kt_double]: Probability of being occupied
     * Return:
     * kt_double: Information content
     */
    return -(prob * log2(prob)) - ((1 - prob) * log2(1 - prob));
}

kt_double InformationEstimates::calculateProbabilityFromLogOdds(kt_double log)
{
    /**
     * Map Log-Odds into probability
     * Arguments:
     * log [kt_double]: Log-Odds
     * Return:
     * kt_double: Probability
     */
    return (exp(log) / (1 + exp(log)));
}

void InformationEstimates::insertMeasurementOutcome(map_tuple tuple, kt_double probability, std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> &map)
{
    int key_free, key_occ, key_unk;
    std::tie(key_free, key_occ, key_unk) = tuple;

    if (map[std::make_tuple(key_free, key_occ, key_unk)])
    {
        map[std::make_tuple(key_free, key_occ, key_unk)] += probability;
    }
    else
    {
        map[std::make_tuple(key_free, key_occ, key_unk)] = probability;
    }
}

std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple> InformationEstimates::computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>> &meas_outcm)
{
    /**
     * Implementation of algorithm 1
     * Compute all the possible combinations of a grid cell, given a set of measurement outcomes
     * Arguments:
     * meas_outcm [std::vector<std::vector<kt_double>>]: Vector of measurement outcomes in the form {p_free, p_occ, p_unk}
     * Return:
     * std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple>: Map of combination, it contains the combination and its probability
     */

    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> probabilities_map;
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> temp_map;

    // Clear the temporal map
    temp_map.clear();
    probabilities_map.clear();

    // Root
    if (meas_outcm.size() == 0)
    {
        probabilities_map[std::make_tuple(0, 0, 0)] = 1.0f;
    }

    // First measurement
    probabilities_map[std::make_tuple(1, 0, 0)] = meas_outcm[0][0]; // Probability free
    probabilities_map[std::make_tuple(0, 1, 0)] = meas_outcm[0][1]; // Probability occupied
    probabilities_map[std::make_tuple(0, 0, 1)] = meas_outcm[0][2]; // Probability unknown

    for (int r = 1; r < meas_outcm.size(); ++r)
    {
        // Measurement outcome probability
        kt_double free_prop = meas_outcm[r][0];
        kt_double occ_prop = meas_outcm[r][1];
        kt_double unk_prop = meas_outcm[r][2];

        // Temporal map to only take the last level combination
        temp_map = probabilities_map;
        probabilities_map.clear();

        for (auto &combination : temp_map)
        {
            int key_free, key_occ, key_unk;
            std::tie(key_free, key_occ, key_unk) = combination.first;

            // Adding next level measurement outcomes
            insertMeasurementOutcome(std::make_tuple(key_free + 1, key_occ, key_unk), combination.second * free_prop, probabilities_map);
            insertMeasurementOutcome(std::make_tuple(key_free, key_occ + 1, key_unk), combination.second * occ_prop, probabilities_map);
            insertMeasurementOutcome(std::make_tuple(key_free, key_occ, key_unk + 1), combination.second * unk_prop, probabilities_map);
        }
        temp_map.clear();
    }
    return probabilities_map;
}

kt_double InformationEstimates::measurementOutcomeEntropy(map_tuple const &meas_outcome)
{
    /**
     * Calculate the measurement outcome entropy
     * Calculate Log-Odds from initial probability guess
     * Calculate the probability from those logs
     * Calculate the entropy with the retrieved probability
     * Arguments:
     * meas_outcome [map_tuple]: Measurement outcome in the form {p_free, p_occ, p_unk}
     * Return:
     * kt_double: Measurement outcome entropy
     */
    int num_free, num_occ, num_unk;
    std::tie(num_free, num_occ, num_unk) = meas_outcome;
    kt_double log_occ = (num_free * l_free) + (num_occ * l_occ) - ((num_free + num_occ - 1) * l_o);
    kt_double prob_occ = calculateProbabilityFromLogOdds(log_occ);
    return calculateInformationContent(prob_occ);
}

kt_double InformationEstimates::calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2)
{
    /**
     * Calculate the mass probability of a cell being observed by a given measurement
     * Arguments:
        * range_1 [kt_double]: Lower bound
        * range_2 [kt_double]: Upper bound
     * Return:
        * kt_double: Mass probability
    */
    range_1 = (range_1 > m_max_sensor_range) ? m_max_sensor_range : range_1;
    range_2 = (range_2 > m_max_sensor_range) ? m_max_sensor_range : range_2;

    // 12.2 because that is roughly the max range when lambda is equals to 1
    kt_double lambda = 12.2 / m_max_sensor_range;

    return -(exp(-lambda * range_2) - exp(-lambda * range_1));
}

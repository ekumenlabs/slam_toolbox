#include <math.h>
#include "slam_toolbox/experimental/information_estimates.hpp"

InformationEstimates::InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu)
{
    m_max_sensor_range = sensor_range;
    m_cell_resol = resolution;
    m_obs_lambda = lambda;
    m_obs_nu = nu;
}

InformationEstimates::InformationEstimates()
{
    m_max_sensor_range = 25.0;
    m_cell_resol = 0.1;
    m_obs_lambda = 0.35;
    m_obs_nu = 0.28;
}


kt_double InformationEstimates::findMutualInfo(std::vector<karto::LocalizedRangeScan*> const& range_scans)
{
    /*
        For your reference @hidmic

        The loop I am working on.

        1. Find the grid size based on the ranges of a group of scans.
            2. Select one cell.
                3. Find which scans can see that cell.
                    4. Select one scan from the group.
                        5. Find the beam that sees the selected cell.
                            6. Calculate and append the probability of that cell being seen.

        Note for point 3:
            I am not sure if we should have another loop right before this point.
                - We should extract one scan from the group of scans.
                - Then execute steps 3, 4, 5, 6.
            Following it, we will be able to find what is the laser that contains less mutual information.

            Another alternative for facing this issue could be: Have a external loop for extracting one scan
            at a time, and then use the ramining scans to execute all the steps explained above.
    */

    // I took the first scan for reference. From this one I extract the min an max points (Poses)
    m_low_x = range_scans[0]->GetCorrectedPose().GetX();
    m_low_y = range_scans[0]->GetCorrectedPose().GetY();

    m_high_x = range_scans[0]->GetCorrectedPose().GetX();
    m_high_y = range_scans[0]->GetCorrectedPose().GetY();

    // Iterating through the group of scans to rescale the grid
    for (const auto & scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();
        // std::cout << "Processing pose: " << pose.GetX() << ", " << pose.GetY() << std::endl;

        // Finding the closest pose
        m_low_x = pose.GetX() < m_low_x ? pose.GetX() : m_low_x;
        m_low_y = pose.GetY() < m_low_y ? pose.GetY() : m_low_y;

        // Finding the farthest pose
        m_high_x = pose.GetX() > m_high_x ? pose.GetX() : m_high_x;
        m_high_y = pose.GetY() > m_high_y ? pose.GetY() : m_high_y;
    }

    // Margins for the closest points
    m_low_x -= m_max_sensor_range;
    m_low_y -= m_max_sensor_range;

    // Margins for the farthest points
    m_high_x += m_max_sensor_range;
    m_high_y += m_max_sensor_range;

    // ===== Grid resizing =====
    // Map dimensions
    kt_double dist_x = std::fabs(m_high_x) + std::fabs(m_low_x);
    kt_double dist_y = std::fabs(m_high_y) + std::fabs(m_low_y);

    // Get the number of cells
    int n_cells_x = static_cast<int>(dist_x / m_cell_resol);
    int n_cells_y = static_cast<int>(dist_y / m_cell_resol);

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);

    int low_x = range_scans[0]->GetCorrectedPose().GetX();
    int low_y = range_scans[0]->GetCorrectedPose().GetY();

    int high_x = range_scans[0]->GetCorrectedPose().GetX();
    int high_y = range_scans[0]->GetCorrectedPose().GetY();

    // ===== Main Loop =====
    // I would need to move it N times
    // Cause here I do remove one scan at a time and call the following loop

    // Iterating through the cells
    for (int i = 0; i < n_cells_x; ++i)
    {
        for (int j = 0; j < n_cells_y; ++j)
        {
            // Iterating through the scans to find which one hits this cell
            for (const auto & scan : range_scans)
            {
                karto::Pose2 pose = scan->GetCorrectedPose();
                karto::Vector2<int> min_cell = utils::grid_operations::getGridPosition(min_pt, m_cell_resol);
                karto::Vector2<int> max_cell = utils::grid_operations::getGridPosition(max_pt, m_cell_resol);

                // kt_double cell_y {(min_cell[1] + max_cell[1]) / 2};
                // kt_double cell_x {(min_cell[0] + max_cell[0])/ 2};
                // if (sqrt(pow((pose[0] - cell_x), 2) + pow((pose[1] - cell_y), 2)) > m_max_sensor_range)
                // {
                //     continue;
                // }

                // Lowest X point
                kt_double low_pt_x = pose.GetX() - m_max_sensor_range;
                low_pt_x = low_pt_x < 0.0 ? 0.0 : low_pt_x;
                // Maximum X point
                kt_double max_pt_x = pose.GetX() + m_max_sensor_range;
                max_pt_x = max_pt_x > dist_x ? dist_x : max_pt_x;
                // Lowest Y point
                kt_double low_pt_y = pose.GetY() - m_max_sensor_range;
                low_pt_y = low_pt_y < 0.0 ? 0.0 : low_pt_y;
                // Maximum Y point
                kt_double max_pt_y = pose.GetY() + m_max_sensor_range;
                max_pt_y = max_pt_y > dist_y ? dist_y : max_pt_y;

                karto::Vector2<kt_double> min_pt{low_pt_x, low_pt_y};
                karto::Vector2<kt_double> max_pt{max_pt_x, max_pt_y};

                // This assertion is just to avoid evaluation
                if ((i <= min_pt.GetX() || i >= max_pt.GetX()) || (j <= min_pt.GetY() || j >= max_pt.GetY()))
                {
                    continue;
                }

                // Current scan laser readings
                karto::PointVectorDouble laser_readings = scan->GetPointReadings(true);

                karto::Pose2 robot_pose_raw = scan->GetCorrectedPose();
                karto::Pose2 robot_pose(robot_pose_raw.GetX() - m_low_x, robot_pose_raw.GetY() - m_low_y, robot_pose_raw.GetHeading());
                karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose.GetPosition(), m_cell_resol);

                // ----- Point 3 -----
                // NIT: @kmilo7204 can use the range based for here
                for (int i = 0; i < laser_readings.size(); ++i)
                {
                    karto::Vector2<kt_double> transformed_laser{laser_readings[i].GetX() - m_low_x, laser_readings[i].GetY() - m_low_y};

                    karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(transformed_laser, m_cell_resol);
                    kt_double laser_to_grid_dist =
                        sqrt(pow((transformed_laser[0] - kt_double(beam_grid[0])), 2) + pow((transformed_laser[1] - kt_double(beam_grid[1])), 2));

                    // Compare the point here, the transformed laser and the cell position
                    // If not fulfill the condition I wrote above, continue
                    if (laser_to_grid_dist > m_cell_resol)
                    {
                        continue;
                    }

                    // Inidividual cell limits
                    kt_double limit_x = i * m_cell_resol;
                    kt_double limit_y = j * m_cell_resol;

                    std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                    robot_pose.GetPosition(), transformed_laser, robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

                    // In case we did not find a measurement
                    if (intersections.first.size() == 0)
                    {
                        continue;
                    }

                    // Enter (d1) and Exit (d2) distances
                    std::vector<kt_double> distances;
                    for (int k = 0; k < intersections.first.size(); ++k)
                    {
                        // From robot position to intersection points
                        karto::Vector2<kt_double> intersection{intersections.first[k], intersections.second[k]};
                        kt_double distance = robot_pose.GetPosition().Distance(intersection);
                        distances.push_back(distance);
                    }

                    // Measurement outcomes vector {Pfree, Pocc, Pun}
                    std::vector<kt_double> probabilities {
                        calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
                        calculateScanMassProbabilityBetween(distances[0], distances[1]),
                        calculateScanMassProbabilityBetween(0.0f, distances[0])
                    };

                    // Appending new measurement outcomes for the current cell
                    appendCellProbabilities(probabilities, cell);
                }
            }
        }
    }
}

std::vector<kt_double> InformationEstimates::findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans)
{
    /**
     * Find and return the Laser Scan index which provide the less mutual information in a set of Laser Scans
     * Arguments:
        * range_scans [std::vector<karto::LocalizedRangeScan*>]: Vector of LocalizedRangeScan pointers
     * Return:
        * std::vector<kt_double>: Vector containing the mutual information fo each of the processed laser scans
    */

    // Get the robot poses for creating our new local coordinate system
    for (const auto & scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();

        // Finding the closest pose
        m_low_x = pose.GetX() < m_low_x ? pose.GetX() : m_low_x;
        m_low_y = pose.GetY() < m_low_y ? pose.GetY() : m_low_y;

        // Finding the farthest pose
        m_high_x = pose.GetX() > m_high_x ? pose.GetX() : m_high_x;
        m_high_y = pose.GetY() > m_high_y ? pose.GetY() : m_high_y;
    }

    // Margins for the closest points
    m_low_x -= m_max_sensor_range;
    m_low_y -= m_max_sensor_range;

    // Margins for the farthest points
    m_high_x += m_max_sensor_range;
    m_high_y += m_max_sensor_range;

    // Map dimensions
    kt_double dist_x = std::fabs(m_high_x) + std::fabs(m_low_x);
    kt_double dist_y = std::fabs(m_high_y) + std::fabs(m_low_y);

    // Number of X and Y cells
    int n_cells_x = static_cast<int>(dist_x / m_cell_resol);
    int n_cells_y = static_cast<int>(dist_y / m_cell_resol);

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);

    std::vector<kt_double> scans_mut_info;
    scans_mut_info.reserve(range_scans.size());

    // Total mutual information
    kt_double map_mut_info = mutualInformationFromScans(range_scans);
    // Clearing the cells for the next time it is called
    m_mutual_grid.setZero();


    for (int i = 0; i < range_scans.size(); ++i)
    {
        // Find mutual information when extracting one laser scan
        kt_double temp_mut_info = mutualInformationFromScans(range_scans, true, i);
        scans_mut_info.emplace_back(map_mut_info - temp_mut_info);
        // Clearing the cells for the next time it is called
        m_mutual_grid.setZero();
    }

    return scans_mut_info;
}


kt_double InformationEstimates::mutualInformationFromScans(std::vector<karto::LocalizedRangeScan*> const& range_scans, bool ignore_scan, int scan_idx)
{
    /**
     * Append measured porbabilities for the given cell
     * Arguments:
        * range_scans [std::vector<karto::LocalizedRangeScan*>]: Vector of LocalizedRangeScan pointers
        * ignore_scan [bool]: Indicates the type of calculation
        * scan_idx [int]: Index within the set to be ignored
     * Return:
        * kt_double
    */

    int curr_idx = 0;
    for (auto & scan : range_scans)
    {
        if (curr_idx == scan_idx && ignore_scan)
        {
            ++curr_idx;
            continue;
        }

        karto::Pose2 robot_pose_raw = scan->GetCorrectedPose();

        // Moving the pose to our local system
        karto::Pose2 robot_pose_tf(robot_pose_raw.GetX() - m_low_x, robot_pose_raw.GetY() - m_low_y, robot_pose_raw.GetHeading());

        karto::PointVectorDouble laser_readings = scan->GetPointReadings(true);
        karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose_tf.GetPosition(), m_cell_resol);

        // Set as false the current boolean map
        m_visited_grid.setZero();

        for (int i = 0; i < laser_readings.size(); ++i)
        {
            karto::Vector2<kt_double> tf_laser{laser_readings[i].GetX() - m_low_x, laser_readings[i].GetY() - m_low_y};

            // Laser final cell
            karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(tf_laser, m_cell_resol);

            // Visited cells by this laser beam
            std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(robot_grid, beam_grid);

            for (auto & cell : cells)
            {
                // Inidividual cell limits
                kt_double limit_x = cell.GetX() * m_cell_resol;
                kt_double limit_y = cell.GetY() * m_cell_resol;

                std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                    robot_pose_tf.GetPosition(), tf_laser, robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

                if (intersections.first.size() == 0)
                    continue;

                // Enter (d1) and Exit (d2) distances
                std::vector<kt_double> distances;
                for (int k = 0; k < intersections.first.size(); ++k)
                {
                    // From robot position to intersection points
                    karto::Vector2<kt_double> intersection{intersections.first[k], intersections.second[k]};
                    kt_double distance = robot_pose_tf.GetPosition().Distance(intersection);
                    distances.push_back(distance);
                }

                // Measurement outcomes vector {Pfree, Pocc, Pun}
                std::vector<kt_double> probabilities {
                    calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
                    calculateScanMassProbabilityBetween(distances[0], distances[1]),
                    calculateScanMassProbabilityBetween(0.0f, distances[0])
                };

                // Appending new measurement outcomes for the current cell
                appendCellProbabilities(probabilities, cell);

                // Compute all the possible combinations for the current cell - algorithm 1
                std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(m_cell_probabilities.at(cell));

                kt_double cell_mutual_inf = 0.0f;
                for (auto& pair : meas_out_prob)
                {
                    cell_mutual_inf +=  pair.second * utils::probability_operations::calculateMeasurementOutcomeEntropy(pair.first);
                }

                // Mutual information of cell x, y given a set of measurements
                m_mutual_grid(cell.GetX(), cell.GetY()) = 1.0 - cell_mutual_inf;

            }
        }
        ++curr_idx;
    }
    m_cell_probabilities.clear();

    return m_mutual_grid.sum();
}

void InformationEstimates::appendCellProbabilities(std::vector<kt_double>& measurements, karto::Vector2<int> const & cell)
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
            cell, {measurements}));
        m_visited_grid(cell.GetX(), cell.GetY()) = 1;
    }
    else
    {
        if(m_visited_grid(cell.GetX(), cell.GetY()) == 1)
        {
            // Compare the unknown probability, the smallest it is the most information we will have
            // from the occupied or free state
            int idx = it_cell->second.size();
            if(measurements[2] < it_cell->second[idx][2])
            {
                // Replacing
                it_cell->second[idx][0] = measurements[0];
                it_cell->second[idx][1] = measurements[1];
                it_cell->second[idx][2] = measurements[2];
            }
        }
        else
        {
            // Cell is already in the map, only add the next measurement outcome
            it_cell->second.push_back({measurements[0], measurements[1], measurements[2]});
            m_visited_grid(cell.GetX(), cell.GetY()) = 1;
        }
    }
}

void InformationEstimates::insertMeasurementOutcome(map_tuple tuple, kt_double probability, std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple>& map)
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

std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple> InformationEstimates::computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm)
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

        for (auto & combination : temp_map)
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

    return m_obs_nu * (exp(-m_obs_lambda*range_1) - exp(-m_obs_lambda*range_2));
}

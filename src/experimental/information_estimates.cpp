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

std::vector<kt_double> InformationEstimates::findMutualInfo(std::vector<karto::LocalizedRangeScan *> const &range_scans)
{
    std::cout << "Vector size: " << range_scans.size() << std::endl;
    /*
        Note: In fact, I should calculate the mutual information with all the elements
        and then calculate the mutual information of each element.

        -- Keep this in mind:
        Entre menos informacion mutua apora, mayor va a ser el resultado total del la operacion al extraerlo.
    */
    std::vector<kt_double> result_vector;
    m_low_x = range_scans[0]->GetCorrectedPose().GetX();
    m_low_y = range_scans[0]->GetCorrectedPose().GetY();

    m_high_x = range_scans[0]->GetCorrectedPose().GetX();
    m_high_y = range_scans[0]->GetCorrectedPose().GetY();


    for (const auto &scan : range_scans)
    {
        if (scan == nullptr)
        {
            continue;
        }

        karto::Pose2 pose = scan->GetCorrectedPose();
        karto::PointVectorDouble laser_readings = scan->GetPointReadings(false);

        kt_double minimumAngle = scan->GetLaserRangeFinder()->GetMinimumAngle();
        kt_double angularResolution = scan->GetLaserRangeFinder()->GetAngularResolution();

        kt_double rangeReading = scan->GetRangeReadings()[0];
        int size = scan->GetRangeReadingsVector().size();

        // Finding the lower limit pose
        m_low_x = std::min(pose.GetX(), m_low_x);
        m_low_y = std::min(pose.GetY(), m_low_y);

        // Finding the higher limit pose
        m_high_x = std::max(pose.GetX(), m_high_x);
        m_high_y = std::max(pose.GetY(), m_high_y);
    }

    // Map dimensions
    kt_double dim_x = std::fabs(m_high_x - m_low_x) + (2 * m_max_sensor_range);
    kt_double dim_y = std::fabs(m_high_y - m_low_y) + (2 * m_max_sensor_range);

    // std::cout << "X Dimension: " << dim_x << std::endl;
    // std::cout << "Y Dimension: " << dim_y << std::endl;

    // Get the number of cells
    int n_cells_x = static_cast<int>(dim_x / m_cell_resol);
    int n_cells_y = static_cast<int>(dim_y / m_cell_resol);

    // Resize the grid with new number of cells
    m_mutual_grid.resize(n_cells_x, n_cells_y);
    m_visited_grid.resize(n_cells_x, n_cells_y);

    // for (int n = 0; n < 1; ++n)
    int skip_laser = 0;
    for (int n = 0; n < range_scans.size(); ++n)
    {
        m_mutual_grid.setZero();
        m_visited_grid.setZero();

        karto::Pose2 robot_pose_raw = range_scans[n]->GetCorrectedPose();
        // Get the current robot pose and extract the extra range for locating it at the lowest X and Y
        karto::Pose2 local_grid_robot_pose{robot_pose_raw.GetX() - (m_low_x - m_max_sensor_range), robot_pose_raw.GetY() - (m_low_y - m_max_sensor_range), robot_pose_raw.GetHeading()};
        // std::cout << "Local robot pose: " << local_grid_robot_pose.GetX() << ", " << local_grid_robot_pose.GetY() << std::endl;

        // @NIT
        // This is an ideal point to add an assertion, since the robot pose must be positive in all of the scenarios

        // Minimum X point
        kt_double lower_limit_x = std::max(0.0, local_grid_robot_pose.GetX() - m_max_sensor_range);
        kt_double upper_limit_x = std::min(dim_x, local_grid_robot_pose.GetX() + m_max_sensor_range);
        kt_double lower_limit_y = std::max(0.0, local_grid_robot_pose.GetY() - m_max_sensor_range);
        kt_double upper_limit_y = std::min(dim_y, local_grid_robot_pose.GetY() + m_max_sensor_range);

        karto::Vector2<kt_double> lower_limit{lower_limit_x, lower_limit_y};
        karto::Vector2<kt_double> upper_limit{upper_limit_x, upper_limit_y};

        // @NIT
        // I can use an rvalue reference here
        karto::Vector2<int> lower_limit_cell = utils::grid_operations::getGridPosition(lower_limit, m_cell_resol);
        karto::Vector2<int> upper_limit_cell = utils::grid_operations::getGridPosition(upper_limit, m_cell_resol);

        karto::Vector2<int> local_robot_cell = utils::grid_operations::getGridPosition(local_grid_robot_pose.GetPosition(), m_cell_resol);

        // std::cout << " =================== " << std::endl;
        // std::cout << "Lower limit cell: " << lower_limit_cell.GetX() << ", " << lower_limit_cell.GetY() << std::endl;
        // std::cout << "Upper limit cell: " << upper_limit_cell.GetX() << ", " << upper_limit_cell.GetY() << std::endl;

        karto::LaserRangeFinder *laser_range_finder = range_scans[n]->GetLaserRangeFinder();
        kt_double range_threshold = laser_range_finder->GetRangeThreshold();
        kt_double angle_increment = laser_range_finder->GetAngularResolution();
        kt_double min_angle = laser_range_finder->GetMinimumAngle();
        kt_double max_angle = laser_range_finder->GetMaximumAngle();
        kt_double total_range = max_angle - min_angle;

        std::vector<karto::Vector2<int>> vis_cells;

        bool printed = false;
        for (int i = lower_limit_cell.GetX(); i < upper_limit_cell.GetX(); ++i)
        {
            for (int j = lower_limit_cell.GetY(); j < upper_limit_cell.GetY(); ++j)
            {
                // Here I would need to find which of the scans can see this cell
                // this is with the distance we have found
                // Euclidean distance between cell center and robot position
                for (int l = 0; l < range_scans.size(); ++l)
                {
                    // Evaluate if the current laser that will be evaluated should be excluded
                    if (n == l)
                    {
                        // Current laser will noe be evaluated
                        if (!printed)
                        {
                            // std::cout << "Skipping: " << l << std::endl;
                            printed = true;
                        }
                        continue;
                    }

                    // Evaluate if the laser can see the current cell
                    karto::Vector2<kt_double> cell_center{i * m_cell_resol + (m_cell_resol / 2), j * m_cell_resol + (m_cell_resol / 2)};
                    kt_double distance_to_cell = sqrt(pow(local_grid_robot_pose.GetX() - cell_center.GetX(), 2) + pow(local_grid_robot_pose.GetY() - cell_center.GetY(), 2));
                    if(distance_to_cell > m_max_sensor_range)
                    {
                        // It means the robot cannot see the current cell
                        continue;
                    }

                    // Continue the staff here
                    kt_double angle_to_cell = atan2(cell_center.GetY() - local_grid_robot_pose.GetY(), cell_center.GetX() - local_grid_robot_pose.GetX());
                    kt_double robot_heading = robot_pose_raw.GetHeading();
                    kt_double angle_between = angle_to_cell - robot_heading;

                    kt_double angle_diff = atan2(sin(angle_between), cos(angle_between));
                    if (angle_diff < min_angle || angle_diff > max_angle)
                    {
                        // We are outside the laser FOV, but we cannot cancel it
                        // because we will not evaluate other cells
                        continue;
                    }

                    int laser_idx = ((total_range / 2.0) + angle_between) / angle_increment;
                    kt_double laser_read_dist = range_scans[n]->GetRangeReadings()[laser_idx];
                    kt_double laser_read_raw = range_scans[n]->GetRangeReadings()[laser_idx];
                    if (laser_read_dist > range_threshold)
                    {
                        continue;
                    }

                    if (laser_read_dist < distance_to_cell)
                    {
                        continue;
                    }
                    else
                    {
                        laser_read_dist = distance_to_cell;
                    }

                    kt_double angle_to_cell_proj = angle_to_cell > 0 ? angle_to_cell : (2.0 * M_PI + angle_to_cell);
                    kt_double point_x_raw = local_grid_robot_pose.GetX() + (laser_read_raw * cos(angle_to_cell_proj));
                    kt_double point_y_raw = local_grid_robot_pose.GetY() + (laser_read_raw * sin(angle_to_cell_proj));
                    karto::Vector2<int> laser_beam_cell_raw = utils::grid_operations::getGridPosition({point_x_raw, point_y_raw}, m_cell_resol);

                    // Inidividual cell limits
                    kt_double limit_x = i * m_cell_resol;
                    kt_double limit_y = j * m_cell_resol;

                    karto::PointVectorDouble laser_readings = range_scans[n]->GetPointReadings(true);
                    karto::Vector2<kt_double> transformed_laser{laser_readings[laser_idx].GetX() - m_low_x, laser_readings[laser_idx].GetY() - m_low_y};

                    std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                    local_grid_robot_pose.GetPosition(), {point_x_raw, point_y_raw}, local_robot_cell, laser_beam_cell_raw, limit_x, limit_y, m_cell_resol);

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
                        kt_double distance = local_grid_robot_pose.GetPosition().Distance(intersection);
                        distances.push_back(distance);
                        std::sort(distances.begin(), distances.end());
                    }

                    // Measurement outcomes vector {Pfree, Pocc, Pun}
                    std::vector<kt_double> probabilities{
                        calculateScanMassProbabilityBetween(distances[1], m_max_sensor_range),
                        calculateScanMassProbabilityBetween(distances[0], distances[1]),
                        calculateScanMassProbabilityBetween(0.0f, distances[0])};

                    karto::Vector2<int> cell{i, j};
                    vis_cells.push_back(cell);

                    // Appending new measurement outcomes for the current cell
                    appendCellProbabilities(probabilities, cell);
                }
            }
        }

        // std::cout << "<------------------------->" << std::endl;
        for (auto &cell : vis_cells)
        {
            std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(m_cell_probabilities.at(cell));
            kt_double cell_mutual_inf = 0.0f;
            for (auto &pair : meas_out_prob)
            {
                cell_mutual_inf += pair.second * measurementOutcomeEntropy(pair.first);
            }

            // Mutual information of cell x, y given a set of measurements
            m_mutual_grid(cell.GetX(), cell.GetY()) = 1.0 - cell_mutual_inf;
        }
        // std::cout << "Mutual information: " << m_mutual_grid.sum() << std::endl;
        result_vector.push_back(m_mutual_grid.sum());

        // Puedo tener un vector que me almacene el resultado
        // Entonces el push va a ser secuencial
        // Cada resultado va a demostrar la informacion mutua calculada cuando se toma todo el grupo - el scan de interes
        // Por ejemplo vct[0] = Index 0 es el calculo de la informacion mutua excluyendo el primer scan
        // Por ejemplo vct[1] = Index 1 es el calculo de la informacion mutua excluyendo el segundo scan
        // Por ejemplo vct[2] = Index 2 es el calculo de la informacion mutua excluyendo el tercer scan

        m_cell_probabilities.clear();
    }
    return result_vector;
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

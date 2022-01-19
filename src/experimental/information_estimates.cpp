#include <math.h>
#include "slam_toolbox/experimental/information_estimates.hpp"

#include <iostream>

InformationEstimates::InformationEstimates(kt_double sensor_range, kt_double resolution, kt_double lambda, kt_double nu)
{
    m_max_sensor_range = sensor_range;
    m_cell_resol = resolution;
    m_obs_lambda = lambda;
    m_obs_nu = nu;

    m_map_dist = 200.0f;
    m_num_cells = static_cast<int>(m_map_dist / m_cell_resol);

    m_info_grid.resize(m_num_cells, m_num_cells);
    m_mutual_grid.resize(m_num_cells, m_num_cells);
    m_visited_grid.resize(m_num_cells, m_num_cells);
}

InformationEstimates::InformationEstimates()
{
    m_max_sensor_range = 5.0;
    m_cell_resol = 0.05;
    m_obs_lambda = 0.35;
    m_obs_nu = 0.28;

    m_map_dist = 200.0f;
    m_num_cells = static_cast<int>(m_map_dist / m_cell_resol);

    m_info_grid.resize(m_num_cells, m_num_cells);
    m_mutual_grid.resize(m_num_cells, m_num_cells);
    m_visited_grid.resize(m_num_cells, m_num_cells);
}

std::tuple<int, kt_double> InformationEstimates::findLeastInformativeLaser(std::vector<karto::LocalizedRangeScan*> const& range_scans)
{
    /**
     * Find and return the Laser Scan index which provide the less mutual information in a set of Laser Scans
     * Arguments:
        * range_scans [std::vector<karto::LocalizedRangeScan*>]: Vector of LocalizedRangeScan pointers
     * Return:
        * std::tuple<int, kt_double>: Tuple containing the index of the LozalizedRangeScan and its corresponding mutual information
    */

    std::vector<kt_double> scans_mut_info;
    scans_mut_info.reserve(range_scans.size());

    for (auto & scan : range_scans)
    {
        karto::Pose2 robot_pose = scan->GetCorrectedPose();
        karto::PointVectorDouble laser_readings = scan->GetPointReadings(true);
        karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose.GetPosition(), m_cell_resol);

        // Set as false the current boolean map
        utils::grid_operations::clearVisitedCells(m_visited_grid);

        for (int i = 0; i < laser_readings.size(); ++i)
        {
            // Laser final cell
            karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(laser_readings[i], m_cell_resol);

            // Visited cells by this laser beam
            std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(robot_grid, beam_grid);

            for (auto & cell : cells)
            {
                // Inidividual cell limits
                kt_double limit_x = cell.GetX() * m_cell_resol;
                kt_double limit_y = cell.GetY() * m_cell_resol;

                std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                    robot_pose.GetPosition(), laser_readings[i], robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

                if (intersections.first.size() == 0)
                    continue;

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

                // Get all the measurement outcomes for the current cell
                std::vector<std::vector<kt_double>> cell_prob = retrieveCellProbabilities(cell);

                // Compute all the possible combinations for the current cell - algorithm 1
                std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(cell_prob);

                kt_double cell_mutual_inf = 0.0f;
                for (auto& pair : meas_out_prob)
                {
                    cell_mutual_inf +=  pair.second * measurementOutcomeEntropy(pair.first);
                }

                // Mutual information of cell x, y given a set of measurements
                updateCellMutualInformation(1.0 - cell_mutual_inf, cell);
            }
        }
        // // Mutual information provided by this laser scan
        // kt_double laser_mut_info = calculateLaserMutualInformation();
        // scans_mut_info.emplace_back(laser_mut_info);

        // utils::grid_operations::clearVisitedCells(m_info_grid);
    }

    // Total mutual information
    kt_double map_mut_info = calculateMapMutualInformation();

    // Clearing the cells for the next time it is called
    utils::grid_operations::clearVisitedCells(m_mutual_grid);

    // Outer loop for know the number of times we should calculate
    for (int i = 0; i < range_scans.size(); ++i)
    {
        // Inner loop for calculating the mutual information without a given measurement
        int scan_idx = 0;
        for (auto & scan : range_scans)
        {
            if (scan_idx == i)
            {
                std::cout << "Breaking: " << i << std::endl;
                ++scan_idx;
                continue;
            }

            karto::Pose2 robot_pose = scan->GetCorrectedPose();
            karto::PointVectorDouble laser_readings = scan->GetPointReadings(true);
            karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose.GetPosition(), m_cell_resol);

            // Set as false the current boolean map
            utils::grid_operations::clearVisitedCells(m_visited_grid);

            for (int i = 0; i < laser_readings.size(); ++i)
            {
                // Laser final cell
                karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(laser_readings[i], m_cell_resol);

                // Visited cells by this laser beam
                std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(robot_grid, beam_grid);

                for (auto & cell : cells)
                {
                    // Inidividual cell limits
                    kt_double limit_x = cell.GetX() * m_cell_resol;
                    kt_double limit_y = cell.GetY() * m_cell_resol;

                    std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(
                        robot_pose.GetPosition(), laser_readings[i], robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

                    if (intersections.first.size() == 0)
                        continue;

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

                    // Get all the measurement outcomes for the current cell
                    std::vector<std::vector<kt_double>> cell_prob = retrieveCellProbabilities(cell);

                    // Compute all the possible combinations for the current cell - algorithm 1
                    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(cell_prob);

                    kt_double cell_mutual_inf = 0.0f;
                    for (auto& pair : meas_out_prob)
                    {
                        cell_mutual_inf +=  pair.second * measurementOutcomeEntropy(pair.first);
                    }

                    // Mutual information of cell x, y given a set of measurements
                    updateCellMutualInformation(1.0 - cell_mutual_inf, cell);
                }
            }
            std::cout << "Scan index: " << scan_idx << std::endl;
            ++scan_idx;
        }

        kt_double temp_mut_info = calculateMapMutualInformation();

        // I(M, Z) - I(M, Z \ {Z[n]})
        std::cout << "Appending: " << map_mut_info << ", " << temp_mut_info << std::endl;
        scans_mut_info.emplace_back(map_mut_info - temp_mut_info);
        utils::grid_operations::clearVisitedCells(m_mutual_grid);

        std::cout << "------------------" << std::endl;
    }

    // Finding the less informative laser scan
    std::vector<kt_double>::iterator it_min;
    it_min = std::min_element(scans_mut_info.begin(), scans_mut_info.end());
    int idx = it_min - scans_mut_info.begin();
    std::cout << "Result: " << idx << ", " << *it_min << std::endl;

    return std::make_tuple(idx, *it_min);

    /*
        - Calculate the mutual information of the whole map with the whole set of measurement outcomes
            - This part is already done
        - Calculate all the mutual information when extracting a measurement
            - Create a loop of size measurements - 1
            - Calculate the mutual information when extracting each element of the array
    */
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
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>>::iterator it_cell;

    it_cell = m_cell_probabilities.find({cell.GetX(), cell.GetY()});
    if (it_cell == m_cell_probabilities.end())
    {
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<std::vector<int>, std::vector<std::vector<kt_double>>>(
            {cell.GetX(), cell.GetY()}, {{measurements[0], measurements[1], measurements[2]}}));
        m_visited_grid(cell.GetX(), cell.GetY()) = 1;
    }
    else
    {
        if(m_visited_grid(cell.GetX(), cell.GetY()) == 1)
        {
            // Compare the unknown probability, the smallest it is the most information we will have
            // from the occupied or free state
            int idx = it_cell->second.size() - 1;
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

kt_double InformationEstimates::calculateInformationContent(kt_double prob)
{
    /**
     * Calculate the information content or self-information based on the probability of cell being occupied
     * Arguments:
        * prob [kt_double]: Probability of being occupied
     * Return:
        * kt_double: Information content
    */
    return - (prob * log2(prob)) -  ((1 - prob) * log2(1 - prob));
}

kt_double InformationEstimates::calculateMapMutualInformation()
{
    /**
     * Calculate the mutual information of the current map, this is the summation of all cells mutual information
     * Arguments:
        * Void
     * Return:
        * kt_double: Map mutual information
    */
    return m_mutual_grid.sum();
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
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> temp_map;

    // The number of measurements
    int k = meas_outcm.size();
    int r = 1;

    kt_double p_free = meas_outcm[0][0];
    kt_double p_occ = meas_outcm[0][1];
    kt_double p_un = meas_outcm[0][2];

    temp_map.clear();

    // Root
    temp_map[std::make_tuple(0, 0, 0)] = 1.0f;

    // First measurement
    temp_map[std::make_tuple(1, 0, 0)] = p_free;
    temp_map[std::make_tuple(0, 1, 0)] = p_occ;
    temp_map[std::make_tuple(0, 0, 1)] = p_un;

    for (int r = 1; r < k; ++r)
    {
        std::vector<map_tuple> tup_vct;
        std::vector<kt_double> acc_prob;

        // Measurement outcome probability
        kt_double free_prop = meas_outcm[r][0];
        kt_double occ_prop = meas_outcm[r][1];
        kt_double un_prop = meas_outcm[r][2];

        for (auto& pair : temp_map)
        {
            // Index
            int idx_free, idx_occ, idx_unk;
            std::tie(idx_free, idx_occ, idx_unk) = pair.first;

            if (idx_free + idx_occ + idx_unk == r)
            {
                // Searching for the current combination in this level
                std::vector<map_tuple>::iterator it_comb;
                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(idx_free + 1, idx_occ, idx_unk));

                // Free
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += pair.second * free_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(idx_free + 1, idx_occ, idx_unk));
                    acc_prob.push_back(pair.second * free_prop);
                }

                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(idx_free, idx_occ + 1, idx_unk));

                // Occupied
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += pair.second * occ_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(idx_free, idx_occ + 1, idx_unk));
                    acc_prob.push_back(pair.second * occ_prop);
                }

                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(idx_free, idx_occ, idx_unk + 1));

                // Unobserved
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += pair.second * un_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(idx_free, idx_occ, idx_unk + 1));
                    acc_prob.push_back(pair.second * un_prop);
                }
            }
        }
        // Inserting the elements into the map
        for (int k = 0; k < tup_vct.size(); ++k)
        {
            temp_map[tup_vct[k]] = acc_prob[k];
        }
    }

    // Leaving in the map only the final outcomes
    std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> out_map;
    for (auto& pair : temp_map)
    {
        int idx_free, idx_occ, idx_unk;
        std::tie(idx_free, idx_occ, idx_unk) = pair.first;
        if (idx_free + idx_occ + idx_unk == k)
        {
            out_map[pair.first] = pair.second;
        }
    }
    return out_map;
}

kt_double InformationEstimates::measurementOutcomeEntropy(map_tuple const& meas_outcome)
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

    return m_obs_nu * (exp(-m_obs_lambda*range_1) - exp(-m_obs_lambda*range_2));
}

std::vector<std::vector<kt_double>> InformationEstimates::retrieveCellProbabilities(karto::Vector2<int> const & cell)
{
    /**
     * Retrieve the cell probabilities (Measurement outcomes)
     * Arguments:
        * cell [karto::Vector2<int>]: Cell coordinates
     * Return:
        * std::vector<std::vector<kt_double>>: Vector of cell probabilities in the form {p_free, p_occ, p_unk} (Measurement outcomes)
    */
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>>::iterator it_cells;
    it_cells = m_cell_probabilities.find({cell.GetX(), cell.GetY()});

    return it_cells->second;
}

kt_double InformationEstimates::calculateLaserMutualInformation()
{
    /**
     * Calculate the laser mutual information considering the map mutual information and the current mutual information calculation
     * Arguments:
        * map_info [kt_double]: Map mutual information
        * curr_info [kt_double]: Current mutual information
     * Return:
        * kt_double: Laser mutual information
    */
    return m_info_grid.sum();
}

void InformationEstimates::updateCellMutualInformation(kt_double mut_inf, karto::Vector2<int> const & cell)
{
    /**
     * Note: Result of the summation 3.12
     * Update the mutual information for each individual cell
     * Arguments:
        * mut_inf [kt_double]: Cell mutual information
        * cell [karto::Vector2<int>]: Cell coordinates
     * Return:
        * Void
    */
    m_mutual_grid(cell.GetX(), cell.GetY()) = mut_inf;
    m_info_grid(cell.GetX(), cell.GetY()) = mut_inf;
}

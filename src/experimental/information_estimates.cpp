#include <math.h>
#include "slam_toolbox/experimental/information_estimates.hpp"

InformationEstimates::InformationEstimates()
{
    m_map_mutual_info = 0.0;
    m_laser_mutual_info = 0.0;
    m_cell_resol = 0.05f; // Map resolution
    m_map_dist = 200.0f; // Total map distance
    m_num_cells = static_cast<int>(m_map_dist / m_cell_resol);

    utils::grid_operations::initializeGrid<kt_double>(m_mutual_grid, m_num_cells, m_num_cells);
    utils::grid_operations::initializeGrid<bool>(m_visited_grid, m_num_cells, m_num_cells);
}

float InformationEstimates::calculateMutualInformation(karto::PointVectorDouble const& laser_readings, karto::Pose2 const& robot_pose)
{
    m_laser_mutual_info = m_laser_mutual_info;
    karto::Vector2<int> robot_grid = utils::grid_operations::getGridPosition(robot_pose.GetPosition(), m_cell_resol);

    // Set as false the current boolean map
    utils::grid_operations::clearVisitedCells(m_visited_grid);

    for (int i = 0; i < laser_readings.size(); ++i)
    {
        // Laser final cell
        karto::Vector2<int> beam_grid = utils::grid_operations::getGridPosition(laser_readings[i], m_cell_resol);

        // Ray tracing for getting the visited cells
        std::vector<int> cells_x, cells_y;
        std::pair<std::vector<int>, std::vector<int>> res_pair = utils::grid_operations::rayCasting(robot_grid, beam_grid);

        cells_x = res_pair.first;
        cells_y = res_pair.second;

        // Visiting the cells
        for (int j = 0; j < cells_x.size(); ++j)
        {
            // Inidividual cell limits
            kt_double limit_x = cells_x[j] * m_cell_resol;
            kt_double limit_y = cells_y[j] * m_cell_resol;

            std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections = utils::grid_operations::computeLineBoxIntersection(robot_pose.GetPosition(), laser_readings[i], robot_grid, beam_grid, limit_x, limit_y, m_cell_resol);

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
                calculateScanMassProbabilityBetween(distances[1], 5.0f),
                calculateScanMassProbabilityBetween(distances[0], distances[1]),
                calculateScanMassProbabilityBetween(0.0f, distances[0])
            };

            // Appending new measurement outcomes for the current cell
            appendCellProbabilities(probabilities, {cells_x[j], cells_y[j]});

            // Get all the measurement outcomes for the current cell
            std::vector<std::vector<kt_double>> cell_prob = retreiveCellProbabilities({cells_x[j], cells_y[j]});

            // Compute all the possible combinations for the current cell - algorithm 1
            std::unordered_map<map_tuple, kt_double, utils::tuple_hash::HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(cell_prob);
            
            kt_double cell_mutual_inf = 0.0f;
            for (auto& pair : meas_out_prob)
            {
                cell_mutual_inf +=  pair.second * measurementOutcomeEntropy(pair.first);
            }

            // Mutual information of cell x, y given a set of measurements                
            updateCellMutualInformation(0.5 - cell_mutual_inf, {cells_x[j], cells_y[j]});
        }
    }
    m_map_mutual_info = calculateMapMutualInformation(); 
    
    // Extract the mutual information provided by this laser scan
    updateLaserMutualInformation();
    return m_map_mutual_info;
}

void InformationEstimates::appendCellProbabilities(std::vector<kt_double>& measurements, std::vector<int> cell)
{
    /*
        To append a new measurement for a specific cell
    */
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>>::iterator it_cell;

    it_cell = m_cell_probabilities.find({cell[0], cell[1]});
    if (it_cell == m_cell_probabilities.end())
    {
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<std::vector<int>, std::vector<std::vector<kt_double>>>(
            {cell[0], cell[1]}, {{measurements[0], measurements[1], measurements[2]}}));
        m_visited_grid[cell[0]][cell[1]] = true;
    }
    else
    {
        if(m_visited_grid[cell[0]][cell[1]] == true)
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
            m_visited_grid[cell[0]][cell[1]] = true;
        }
    }
}

kt_double InformationEstimates::calculateInformationContent(kt_double prob)
{
    /*
        To calculate the information content or self-information
        based on the proability of being occupied
    */
    return - (prob * log2(prob)) -  ((1 - prob) * log2(1 - prob));
}

kt_double InformationEstimates::calculateMapMutualInformation()
{
    /*
        To calculate map mutual information, this is the summation
        of all cells mutual information
    */
    kt_double sum = 0.0f;
    for (int i = 0; i < m_num_cells; ++i)
    {
        for (int j = 0; j < m_num_cells; ++j)
        {
            sum += m_mutual_grid[i][j];
        }
    }
    return sum;
}

kt_double InformationEstimates::calculateProbabilityFromLogOdds(kt_double log)
{
    /*
        To transform the Log-odds into probability
        This should be a free function
    */
    return (exp(log) / (1 + exp(log)));
}

std::unordered_map<InformationEstimates::map_tuple, kt_double, utils::tuple_hash::HashTuple> InformationEstimates::computeMeasurementOutcomesHistogram(std::vector<std::vector<kt_double>>& meas_outcm)
{
    /*
        To compute all the possible combinations of a grid cell, given a set of measurement outcomes
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
    /*
        To calculate the measurement outcome entropy (Measurement outcome in the form <fr, oc, un>)
            - Calculate Log-Odds from initial probability guess
            - Calculate the probability from those logs
            - Calculate the entropy with the retrieved probability
    */
    int num_free, num_occ, num_unk; 
    std::tie(num_free, num_occ, num_unk) = meas_outcome;
    kt_double log_occ = (num_free * l_free) + (num_occ * l_occ) - ((num_free + num_occ - 1) * l_o);
    kt_double prob_occ = calculateProbabilityFromLogOdds(log_occ);
    return calculateInformationContent(prob_occ);
}

kt_double InformationEstimates::calculateScanMassProbabilityBetween(kt_double range_1, kt_double range_2)
{
    /*
        To calculate the mass probability of a cell being observed by a given measurement
    */

    range_1 = (range_1 > m_max_sensor_range) ? m_max_sensor_range : range_1;
    range_2 = (range_2 > m_max_sensor_range) ? m_max_sensor_range : range_2;

    return m_obs_nu * (exp(-m_obs_lambda*range_1) - exp(-m_obs_lambda*range_2));
}

std::vector<std::vector<kt_double>> InformationEstimates::retreiveCellProbabilities(std::vector<int> cell)
{
    /*
        To get all the cell probabilities
    */
    std::map<std::vector<int>, std::vector<std::vector<kt_double>>>::iterator it_cells;
    it_cells = m_cell_probabilities.find({cell[0], cell[1]});

    return it_cells->second; 
}

void InformationEstimates::updateLaserMutualInformation()
{
    m_laser_mutual_info = fabs(m_laser_mutual_info - m_map_mutual_info);
}

void InformationEstimates::updateCellMutualInformation(kt_double mut_inf, std::vector<int> cell)
{
    /*
        To update the mutual information for each individual cell
        This is the result of the summation of 3.12
    */
    m_mutual_grid[cell[0]][cell[1]] = mut_inf;
}

void InformationEstimates::setMaxSensorRange(kt_double const sensor_range)
{
    m_max_sensor_range = sensor_range;
}

void InformationEstimates::setObservationLambda(kt_double const lambda)
{
    m_obs_lambda = lambda;
}

void InformationEstimates::setObservationNu(kt_double const nu)
{
    m_obs_nu = nu;
}

void InformationEstimates::setCellResolution(kt_double const resolution)
{
    m_cell_resol = resolution;
}

void InformationEstimates::setMapDistance(kt_double const distance)
{
    m_map_dist = distance;
}


kt_double InformationEstimates::getMapMutualInformation()
{
    return m_map_mutual_info;
}

kt_double InformationEstimates::getLaserMutualInformation()
{
    return m_laser_mutual_info;
}


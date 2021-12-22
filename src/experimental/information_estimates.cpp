#include <iostream>
#include "slam_toolbox/experimental/information_estimates.hpp"

InformationEstimates::InformationEstimates()
{
    /*
        Need to add the new elements for the constructor
    */
    m_cell_resol = 0.5f; // Map resolution
    m_map_dist = 20.0f; // Total map distance
    m_num_cells = static_cast<int>(m_map_dist / m_cell_resol);

    // Grids initialization (Occupancy and Mutual information)
    m_grid.resize(m_num_cells);
    m_mutual_grid.resize(m_num_cells);
    visited_cells.resize(m_num_cells);
    initializeGrids();

    scannerTest();
}

void InformationEstimates::scannerTest()
{
    /*
        La lectura de menor ganancia

        Dado uno computamos el score - Esta informacion esta en el toolbox
        Buscar scans en el mismo rango de vision

        // Set -> information
        // Nodo -> Que tanta ganancia aporta (Ganancia de informacion)
        Batch elimination (Hiteresis) - Groups

        Lista de adyacencia
    */

    // Loop through the different robot poses
    for (int r = 0; r < robot_poses.size(); ++r)
    {
        std::cout << "-------------------- Robot pose --------------------: " << r << std::endl;

        // Angles {-50, -25, 0, 25, 50} in degrees
        std::vector<double> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

        // Initial point
        std::vector<int> robot_grid_pos = getGridPosition(robot_poses[r][0], robot_poses[r][1]);
        std::cout << "Robot position: " << robot_grid_pos[0] << ", " << robot_grid_pos[1] << std::endl;

        // Set as false the current boolean map
        clearVisitedCells();

        for (int i = 0; i < laser_ranges[r].size(); ++i)
        {
            std::cout << "................ New laser ................" << std::endl;
            std::cout << "Distance: " << laser_ranges[r][i] << ", Angle: " << angles[i] << std::endl;

            // Laser continuous distance
            std::vector<double> laser_end = laserHitDistance(robot_poses[r], laser_ranges[r][i], angles[i]);

            // Laser final cell
            std::vector<int> final_grid_pos = getGridPosition(laser_end[0], laser_end[1]);

            // Ray tracing for getting the visited cells
            std::vector<int> cells_x, cells_y;
            std::pair<std::vector<int>, std::vector<int>> res_pair = rayCasting(robot_grid_pos[0], robot_grid_pos[1], final_grid_pos[0], final_grid_pos[1]);
            cells_x = res_pair.first;
            cells_y = res_pair.second;

            // Visiting the cells
            for (int j = 0; j < cells_x.size(); ++j)
            {
                // Cells visualization
                std::cout << "Current cell: " << cells_x[j] << ", " << cells_y[j] << std::endl;

                // Inidividual cell limits
                double limit_x = cells_x[j] * m_cell_resol;
                double limit_y = cells_y[j] * m_cell_resol;

                std::pair<std::vector<double>, std::vector<double>> intersections = computeLineBoxIntersection(robot_poses[r], laser_end, robot_grid_pos, final_grid_pos, limit_x, limit_y);

                if (intersections.first.size() == 0)
                    continue;

                // Enter (d1) and Exit (d2) distances
                std::vector<double> distances;
                for (int k = 0; k < intersections.first.size(); ++k)
                {
                    // From robot position to intersection points
                    double dist_point = euclideanDistance(robot_poses[r][0], robot_poses[r][1], intersections.first[k], intersections.second[k]);
                    distances.push_back(dist_point);
                }

                // Integral 1: 0 to d1 which is distance from robot pose to first point where the cell is cut
                // Integral 2: d1 which is distance from robot pose to first point where the cell is cut to
                // d2 which is distance from robot pose to second point where the cell is cut
                // Integral 3: d2 which is distance from robot pose to second point where the cell is cut to z_max

                // Measurement outcomes vector {Pfree, Pocc, Pun}
                // Map of tuples
                std::vector<double> probabilities {
                    calculateScanMassProbabilityBetween(distances[1], 5.0f),
                    calculateScanMassProbabilityBetween(distances[0], distances[1]),
                    calculateScanMassProbabilityBetween(0.0f, distances[0])
                };

                // Appending new measurement outcomes for the current cell
                appendCellProbabilities(probabilities, {cells_x[j], cells_y[j]});

                // Get all the measurement outcomes for the current cell
                std::vector<std::vector<double>> cell_prob = retreiveCellProbabilities({cells_x[j], cells_y[j]});

                // Compute all the possible combinations for the current cell - algorithm 1
                std::unordered_map<map_tuple, double, HashTuple> meas_out_prob = computeMeasurementOutcomesHistogram(cell_prob);
                
                double cell_mutual_inf = 0.0f;
                for (auto& pair : meas_out_prob)
                {
                    cell_mutual_inf +=  pair.second * measurementOutcomeEntropy(pair.first);
                }

                // Mutual information of cell x, y given a set of measurements                
                updateCellMutualInformation(0.5 - cell_mutual_inf, {cells_x[j], cells_y[j]});
                std::cout << "++++++++++++++++++++++++" << std::endl;
            }
        }
        double mutual = calculateMapMutualInformation();
        std::cout << "Mutual information: " << mutual << std::endl;
    }
}

void InformationEstimates::clearVisitedCells()
{
    /*
        To clear the visited cell
    */

    std::cout << "Clearing cells " << std::endl;
    for (int i = 0; i < visited_cells.size(); ++i)
    {
        for (int j = 0; j < visited_cells[0].size(); ++j)
        {
            visited_cells[i][j] = false;
        }
    }
}

void InformationEstimates::appendCellProbabilities(std::vector<double>& measurements, std::vector<int> cell)
{
    /*
        To append a new measurement for a specific cell
    */
    std::map<std::vector<int>, std::vector<std::vector<double>>>::iterator it_cell;

    it_cell = m_cell_probabilities.find({cell[0], cell[1]});

    if (it_cell == m_cell_probabilities.end())
    {
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<std::vector<int>, std::vector<std::vector<double>>>(
        {cell[0], cell[1]},
        {{measurements[0], measurements[1], measurements[2]}}
        ));
        visited_cells[cell[0]][cell[1]] = true;
    }
    else
    {
        if(visited_cells[cell[0]][cell[1]] == true)
        {
            /*
                Compare the unknown probability, the smallest it is the most information we will have
                from the occupied or free state
            */
            int idx = it_cell->second.size() - 1;
            if(measurements[2] < it_cell->second[idx][2])
            {
                // std::cout << "Replacing:" << it_cell->second[idx][0] << ", " << it_cell->second[idx][1] << ", " << it_cell->second[idx][2] << std::endl;
                // std::cout << "With:" << measurements[0] << ", " << measurements[1] << ", " << measurements[2] << std::endl;

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
            visited_cells[cell[0]][cell[1]] = true;
        }
    }
}

std::pair<std::vector<int>, std::vector<int>> InformationEstimates::rayCasting(int x_1, int y_1, int x_2, int y_2)
{
    /*
        To find the set of cells hit by a laser beam
        This is based on Bresenham algorithm
    */
    std::vector<int> x_bres;
    std::vector<int> y_bres;

    int x = x_1;
    int y = y_1;

    int delta_x = abs(x_2 - x_1);
    int delta_y = abs(y_2 - y_1);

    int s_x = signum(x_2 - x_1);
    int s_y = signum(y_2 - y_1);
    bool interchange = false;

    if (delta_y > delta_x)
    {
        int temp = delta_x;
        delta_x = delta_y;
        delta_y = temp;
        interchange = true;
    }
    else
    {
        interchange = false;
    }

    int a_res = 2 * delta_y;
    int b_res = 2 * (delta_y - delta_x);
    int e_res = (2 * delta_y) - delta_x;

    x_bres.push_back(x);
    y_bres.push_back(y);

    for (int i = 1; i < delta_x; ++i)
    {
        if (e_res < 0)
        {
            if (interchange)
            {
                y += s_y;
            }
            else
            {
                x += s_x;
            }
            e_res += a_res;
        }
        else
        {
            y += s_y;
            x += s_x;
            e_res += b_res;
        }
        x_bres.push_back(x);
        y_bres.push_back(y);
    }
    // Delete the current robot cell
    x_bres.erase(x_bres.begin());
    y_bres.erase(y_bres.begin()); 
            
    // Adding last hit cell to the set
    x_bres.push_back(x_2);
    y_bres.push_back(y_2);
    return std::pair<std::vector<int>, std::vector<int>>{x_bres, y_bres};
}

double InformationEstimates::calculateInformationContent(double prob)
{
    /*
        To calculate the information content or self-information
        based on the proability of being occupied
    */
    return - (prob * log2(prob)) -  ((1 - prob) * log2(1 - prob));
}

std::vector<double> InformationEstimates::calculateCellIntersectionPoints(
    std::vector<double> const& laser_start, std::vector<double> const& laser_end,
    std::vector<double> cell_start, std::vector<double> cell_end)
{
    /*
        Initial point laser beam: laser_start
        Final point laser beam: laser_end
        Initial point cell: cell_start
        Final point cell: cell_end
    */
    double x1 = laser_start[0];
    double x2 = laser_end[0];
    double x3 = cell_start[0];
    double x4 = cell_end[0];

    double y1 = laser_start[1];
    double y2 = laser_end[1];
    double y3 = cell_start[1];
    double y4 = cell_end[1];

    double den = ((x2-x1)*(y4-y3) - (x4-x3)*(y2-y1));
    if (den == 0.0f)
    {
        // Parallel lines or not intersection at all
        return {};
    }
    else
    {
        double x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2-x1)) / den;
        double y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2-y1)) / den;
        return {x, y};
    }
}

std::pair<std::vector<double>, std::vector<double>> InformationEstimates::computeLineBoxIntersection(
    std::vector<double> const& laser_start, std::vector<double> const& laser_end, 
    std::vector<int> const& robot_grid_pos, std::vector<int> const& final_grid_pos,
    double limit_x, double limit_y)
{
    // Cell limits: min_x, max_x, min_y, max_y
    std::vector<double> cell_limits {limit_x, limit_x + m_cell_resol, limit_y, limit_y + m_cell_resol};

    // Initial points for each of the 4 corners
    std::vector<double> initial_x {limit_x, limit_x, limit_x + m_cell_resol, limit_x + m_cell_resol};
    std::vector<double> initial_y {limit_y, limit_y, limit_y + m_cell_resol, limit_y + m_cell_resol};
    
    // Final points for each of the 4 corners
    std::vector<double> final_x {limit_x + m_cell_resol, limit_x, limit_x + m_cell_resol, limit_x};
    std::vector<double> final_y {limit_y, limit_y + m_cell_resol, limit_y, limit_y + m_cell_resol};

    // Set the new cell limits
    updateCellLimits(initial_x, initial_y, final_x, final_y, limit_x, limit_y, cell_limits, robot_grid_pos, final_grid_pos);

    std::vector<double> inter_x, inter_y;

    for (int k = 0; k < 4; ++k)
    {
        std::vector<double> intersection = calculateCellIntersectionPoints(laser_start, laser_end, {initial_x[k], initial_y[k]}, {final_x[k], final_y[k]});
        if(intersection.size() != 0)
        {
            if ((fabs(intersection[0]) >= (fabs(cell_limits[0]) - 0.001)) &&
            (fabs(intersection[0]) <= (fabs(cell_limits[1]) + 0.001)) &&
            (fabs(intersection[1]) >= (fabs(cell_limits[2]) - 0.001)) &&
            (fabs(intersection[1]) <= (fabs(cell_limits[3]) + 0.001)))
            {
                /*
                    Two points where the beam cuts the cell
                    - A laser beam can cut the cell at least 1 time (Enter)
                    - A laser beam can cut the cell at most 2 times (Enter an exit)
                */
                inter_x.push_back(intersection[0]);
                inter_y.push_back(intersection[1]);
            }
        }
    }
    return std::pair<std::vector<double>, std::vector<double>>{inter_x, inter_y}; 
}

void InformationEstimates::updateCellLimits(
    std::vector<double>& initial_x, std::vector<double>& initial_y, std::vector<double>& final_x, std::vector<double>& final_y,
    double limit_x, double limit_y, std::vector<double>& cell_limits, std::vector<int> const& robot_grid_pos, std::vector<int> const& final_grid_pos)
{
    /*
        To calculate grid grid limits for intersection
    */
    if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] >= robot_grid_pos[1])
    {
        // X greater and Y greater. WRO final points
        final_x[0] = limit_x + m_cell_resol;
        final_x[2] = limit_x + m_cell_resol;

        cell_limits[2] = limit_y;
        cell_limits[3] = limit_y + m_cell_resol;
    }

    if (final_grid_pos[0] >= robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
    {
        // X greater and Y minor. WRO final points
        initial_y[2] = limit_y - m_cell_resol;
        initial_y[3] = limit_y - m_cell_resol;

        final_y[1] = limit_y - m_cell_resol;
        final_y[3] = limit_y - m_cell_resol;

        cell_limits[2] = limit_y - m_cell_resol;
        cell_limits[3] = limit_y;
    }

    if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
    {
        // X minor and Y minor. WRO final points
        initial_x[2] = limit_x - m_cell_resol;
        initial_x[3] = limit_x - m_cell_resol;
        initial_y[2] = limit_y - m_cell_resol;
        initial_y[3] = limit_y - m_cell_resol;

        final_x[0] = limit_x - m_cell_resol;
        final_x[2] = limit_x - m_cell_resol;
        final_y[1] = limit_y - m_cell_resol;
        final_y[3] = limit_y - m_cell_resol;

        cell_limits[0] = limit_x - m_cell_resol;
        cell_limits[1] = limit_x;
        cell_limits[2] = limit_y - m_cell_resol;
        cell_limits[3] = limit_y;
    }
}

double InformationEstimates::calculateLogOddsFromProbability(double probability)
{
    /*
        To calculate the log-odds
        This should be a free function
    */
    return log(probability / (1 - probability));
}

double InformationEstimates::calculateMapMutualInformation()
{
    /*
        To calculate map mutual information, this is the summation
        of all cells mutual information
    */
    double sum = 0.0f;
    for (int i = 0; i < m_num_cells; ++i)
    {
        for (int j = 0; j < m_num_cells; ++j)
        {
            sum += m_mutual_grid[i][j];
        }
    }
    return sum;
}

double InformationEstimates::calculateProbabilityFromLogOdds(double log)
{
    /*
        To transform the Log-odds into probability
        This should be a free function
    */
    return (exp(log) / (1 + exp(log)));
}

std::unordered_map<InformationEstimates::map_tuple, double, InformationEstimates::HashTuple> InformationEstimates::computeMeasurementOutcomesHistogram(std::vector<std::vector<double>>& meas_outcm)
{
    /*
        To compute all the possible combinations of a grid cell, given a set of measurement outcomes
    */
    std::unordered_map<map_tuple, double, HashTuple> temp_map;

    // The number of measurements
    int k = meas_outcm.size(); 
    int r = 1;

    double p_free = meas_outcm[0][0];
    double p_occ = meas_outcm[0][1];
    double p_un = meas_outcm[0][2];

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
        std::vector<double> acc_prob;

        // Measurement outcome probability
        double free_prop = meas_outcm[r][0];
        double occ_prop = meas_outcm[r][1];
        double un_prop = meas_outcm[r][2];

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
    std::unordered_map<map_tuple, double, HashTuple> out_map;
    for (auto& pair : temp_map)
    {
        if (std::get<0>(pair.first) + std::get<1>(pair.first) + std::get<2>(pair.first) == k)
        {
            out_map[pair.first] = pair.second;
        }
    }
    return out_map;
}

double InformationEstimates::euclideanDistance(double x_1, double y_1, double x_2, double y_2)
{
    /*
        To calculate the euclidean distance between two points
    */
    double diff_x = x_2 - x_1;
    double diff_y = y_2 - y_1;

    return sqrt(diff_x*diff_x + diff_y*diff_y);
}

int InformationEstimates::signum(int num)
{
    /*
        To get the sign of an operation, used by Bresenham algorithm
    */
    if (num < 0) return -1; 
    if (num >= 1) return 1;
    return 0;
}

std::vector<int> InformationEstimates::getGridPosition(double x, double y)
{
    /*
        To maps the current position into grid coordinates
    */
    int x_cell = floor((x / m_cell_resol));
    int y_cell = floor((y / m_cell_resol));

    return {x_cell, y_cell};
}

std::vector<double> InformationEstimates::laserHitDistance(std::vector<double> const& robot_pose, double range, double angle)
{
    /*
        To get the distance where the laser beam hits something
            - Applying RBT from the global to the sensor frame
    */
    double x_tf = (range * cos(robot_pose[2] + angle)) + robot_pose[0];
    double y_tf = (range * sin(robot_pose[2] + angle)) + robot_pose[1];

    return {x_tf, y_tf};
}

void InformationEstimates::initializeGrids()
{
    /*
        To create the grid
    */
    for (int i = 0; i < m_num_cells; ++i)
    {
        // Adding columns
        m_grid[i].resize(m_num_cells);
        m_mutual_grid[i].resize(m_num_cells);
        visited_cells[i].resize(m_num_cells);
        for (int j = 0; j < m_num_cells; ++j)
        {
            m_grid[i][j] = 0;
            m_mutual_grid[i][j] = 0.0f;
            visited_cells[i][j] = false;
        }
    }
}

double InformationEstimates::measurementOutcomeEntropy(map_tuple const& meas_outcome)
{
    /*
        To calculate the measurement outcome entropy (Measurement outcome in the form <fr, oc, un>)
            - Calculate Log-Odds from initial probability guess
            - Calculate the probability from those logs
            - Calculate the entropy with the retrieved probability
    */
    int num_free, num_occ, num_unk; 
    std::tie(num_free, num_occ, num_unk) = meas_outcome;
    
    double log_occ = (num_free * l_free) + (num_occ * l_occ) - ((num_free + num_occ - 1) * l_o);
    double prob_occ = calculateProbabilityFromLogOdds(log_occ);
    return calculateInformationContent(prob_occ);
}

double InformationEstimates::calculateScanMassProbabilityBetween(double range_1, double range_2)
{
    /*
        To calculate the mass probability of a cell being observed by a given measurement
    */

    range_1 = (range_1 > m_max_sensor_range) ? m_max_sensor_range : range_1;
    range_2 = (range_2 > m_max_sensor_range) ? m_max_sensor_range : range_2;

    return m_obs_nu * (exp(-m_obs_lambda*range_1) - exp(-m_obs_lambda*range_2));
}

std::vector<std::vector<double>> InformationEstimates::retreiveCellProbabilities(std::vector<int> cell)
{
    /*
        To get all the cell probabilities
    */
    std::map<std::vector<int>, std::vector<std::vector<double>>>::iterator it_cells;
    it_cells = m_cell_probabilities.find({cell[0], cell[1]});

    return it_cells->second; 
}

void InformationEstimates::updateCellMutualInformation(double mut_inf, std::vector<int> cell)
{
    /*
        To update the mutual information for each individual cell
        This is the result of the summation of 3.12
    */
    m_mutual_grid[cell[0]][cell[1]] = mut_inf;
}

void InformationEstimates::setMaxSensorRange(double const sensor_range)
{
    m_max_sensor_range = sensor_range;
}

void InformationEstimates::setObservationLambda(double const lambda)
{
    m_obs_lambda = lambda;
}

void InformationEstimates::setObservationNu(double const nu)
{
    m_obs_nu = nu;
}

void InformationEstimates::setCellResolution(double const resolution)
{
    m_cell_resol = resolution;
}

void InformationEstimates::setMapDistance(double const distance)
{
    m_map_dist = distance;
}

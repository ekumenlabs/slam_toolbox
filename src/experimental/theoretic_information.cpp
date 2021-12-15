#include <iostream>
#include "slam_toolbox/experimental/theoretic_information.hpp"

TeorethicInformation::TeorethicInformation()
{
    std::cout << "Constructor of Theoretic information" << std::endl;
    m_resolution = 0.5f; // Map resolution
    m_map_dist = 20.0f; // Total map distance
    m_num_cells = static_cast<int>(m_map_dist / m_resolution);

    // Grids initialization (Occupancy and Mutual information)
    m_grid.resize(m_num_cells);
    m_mutual_grid.resize(m_num_cells);
    initializeGrids();

    scannerTest();
}

void TeorethicInformation::scannerTest()
{
    // Loop through the different robot poses
    for (int r = 0; r < robot_poses.size(); ++r)
    {
        std::cout << "---------- Robot pose ----------: " << r << std::endl;

        // Angles {-50, -25, 0, 25, 50} in degrees
        std::vector<float> angles{-0.87266f, -0.43633f, 0.0f, 0.43633f, 0.87266f};

        // Initial point
        std::vector<int> robot_grid_pos = getGridPosition(robot_poses[r][0], robot_poses[r][1]);
        std::cout << "Robot position: " << robot_grid_pos[0] << ", " << robot_grid_pos[1] << std::endl;

        // Current yaw + beam angle: -PI/2 (-1.570795) -0.87266 = 2.44345 (-55 degrees)
        for (int i = 0; i < laser_ranges[r].size(); ++i)
        {
            std::cout << "........ New laser ........" << std::endl;
            std::cout << "Distance: " << laser_ranges[r][i] << ", Angle: " << angles[i] << std::endl;

            // Laser continuous distance
            std::vector<float> laser_grid = laserHitDistance(robot_poses[r], laser_ranges[r][i], angles[i]);

            // Laser final cell
            std::vector<int> final_grid_pos = getGridPosition(laser_grid[0], laser_grid[1]);

            // robot_grid_pos[0] // X1 - robot_grid_pos[1] // Y1
            // final_grid_pos[0] // X2 - final_grid_pos[1] // Y2

            // Ray tracing for getting the visited cells
            std::vector<int> cells_x, cells_y;
            std::pair<std::vector<int>, std::vector<int>> res_pair = Bresenham(robot_grid_pos[0], robot_grid_pos[1], final_grid_pos[0], final_grid_pos[1]);
            cells_x = res_pair.first;
            cells_y = res_pair.second;

            // Deleting the current robot cell
            cells_x.erase(cells_x.begin());
            cells_y.erase(cells_y.begin()); 
            
            // Adding last hit cell to the set
            cells_x.push_back(final_grid_pos[0]);
            cells_y.push_back(final_grid_pos[1]);
            

            // Visiting the cells
            for (int j = 0; j < cells_x.size(); ++j)
            {
                // Cells visualization
                std::cout << "Current cell: " << cells_x[j] << ", " << cells_y[j] << std::endl;

                // Inidividual cell limits
                float limit_x = cells_x[j] * m_resolution;
                float limit_y = cells_y[j] * m_resolution;

                // Cell limits: min_x, max_x, min_y, max_y
                std::vector<float> cell_limits {limit_x, limit_x + m_resolution, limit_y, limit_y + m_resolution};

                // Initial points for each of the 4 corners
                std::vector<float> initial_x {limit_x, limit_x, limit_x + m_resolution, limit_x + m_resolution};
                std::vector<float> initial_y {limit_y, limit_y, limit_y + m_resolution, limit_y + m_resolution};
                
                // Final points for each of the 4 corners
                std::vector<float> final_x {limit_x + m_resolution, limit_x, limit_x + m_resolution, limit_x};
                std::vector<float> final_y {limit_y, limit_y + m_resolution, limit_y, limit_y + m_resolution};

                // Set the new cell limits
                updateCellLimits(initial_x, initial_y, final_x, final_y, limit_x, limit_y, cell_limits, robot_grid_pos, final_grid_pos);

                std::vector<float> inter_x, inter_y;
                for (int k = 0; k < 4; ++k)
                {
                    std::vector<float> intersection = calculateCellIntersectionPoints(robot_poses[r], laser_grid, {initial_x[k], initial_y[k]}, {final_x[k], final_y[k]});
                    if(intersection.size() != 0)
                    {
                        // If the laser and a cell intersects, we need to make sure it happens in the right bounds
                        if ((abs(intersection[0]) >= abs(cell_limits[0] - 0.01f)) &&
                        (abs(intersection[0]) <= abs(cell_limits[1] + 0.01f)) &&
                        (abs(intersection[1]) >= abs(cell_limits[2] - 0.01f)) &&
                        (abs(intersection[1]) <= abs(cell_limits[3] + 0.01f)))
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

                // When a cell is marked by Bresenham but there is not intersection points
                if (inter_x.size() == 0)
                    continue;

                // Enter (d1) and Exit (d2) distances
                std::vector<float> distances;
                for (int k = 0; k < inter_x.size(); ++k)
                {
                    // From robot position to intersection points
                    float dist_point = euclideanDistance(robot_poses[r][0], robot_poses[r][1], inter_x[k], inter_y[k]);
                    distances.push_back(dist_point);
                }

                // Integral 1: 0 to d1 which is distance from robot pose to first point where the cell is cut
                // Integral 2: d1 which is distance from robot pose to first point where the cell is cut to
                // d2 which is distance from robot pose to second point where the cell is cut
                // Integral 3: d2 which is distance from robot pose to second point where the cell is cut to z_max

                // Measurement outcomes vector {Pfree, Pocc, Pun}
                std::vector<float> probabilities {
                    probabilityFromObservation(distances[1], 5.0f),
                    probabilityFromObservation(distances[0], distances[1]),
                    probabilityFromObservation(0.0f, distances[0])
                };

                // Assigning the cells
                m_cell_x = cells_x[j];
                m_cell_y = cells_y[j];

                // Appending new measurement outcomes for the current cell
                appendCellProbabilities(probabilities, {cells_x[j], cells_y[j]});

                // Get all the measurement outcomes for the current cell
                std::vector<std::vector<float>> meas_outcomes = retreiveMeasurementOutcomes({cells_x[j], cells_y[j]});

                // Compute all the possible combinations for the current cell - algorithm 1
                std::unordered_map<map_tuple, float, HashTuple> meas_out_prob = computeProbabilities(meas_outcomes);

                // Calculate 3.12
                std::unordered_map<map_tuple, float, HashTuple>::iterator it_mutual;
                std::cout << "Number of measurements: " << meas_outcomes.size() << std::endl;
                float cell_mutual_inf = 0.0f;
                for (it_mutual = meas_out_prob.begin(); it_mutual != meas_out_prob.end(); ++it_mutual)
                {
                    // Interested in the final measurement outcomes
                    if (std::get<0>(it_mutual->first) + std::get<1>(it_mutual->first) + std::get<2>(it_mutual->first) == meas_outcomes.size())
                    {
                        // measurementOutcomeEntropy is negative
                        cell_mutual_inf +=  it_mutual->second * measurementOutcomeEntropy(it_mutual->first);
                        // std::cout << "----------Entropy---------- " << measurementOutcomeEntropy(it_mutual->first) << std::endl;
                        // std::cout << "----------Entropy---------- " << 0.3*log2(0.3) << std::endl;
                        // std::cout << "----------Entropy---------- " << calculateEntropy(0.7) << std::endl;
                        // std::cout << "----------Entropy---------- " << calculateEntropy(probabilityFromLogs(calculateLogs(0.7f))) << std::endl;
                    }

                }

                // Here should be the H(C) = 0.5 : 0.5 - SUM(P*H)
                std::cout << "Cell mutual information: " << 0.5 - cell_mutual_inf << std::endl;
                // Mutual information of cell x, y given a set of measurements
                updateCellMutualInformation(0.5 - cell_mutual_inf, {cells_x[j], cells_y[j]});
                std::cout << "++++++++++++++++++++++++" << std::endl;
            }
        }
        float mutual = calculateMapMutualInformation();
        std::cout << "Mutual information: " << mutual << std::endl;
    }
}

void TeorethicInformation::appendCellProbabilities(std::vector<float>& measurements, std::vector<int> cell)
{
    /*
        To append a new measurement for a specific cell
    */

    // Iterator for getting the cell
    std::map<std::vector<int>, std::vector<std::vector<float>>>::iterator it_cell;

    it_cell = m_cell_probabilities.find({cell[0], cell[1]});
    if (it_cell == m_cell_probabilities.end())
    {
        // Cell is not present in the map, so append it
        m_cell_probabilities.insert(std::pair<std::vector<int>, std::vector<std::vector<float>>>(
        {cell[0], cell[1]},
        {{measurements[0], measurements[1], measurements[2]}}
        ));
    }
    else
    {
        // Cell is already in the map, only add the next measurement outcome
        it_cell->second.push_back({measurements[0], measurements[1], measurements[2]});
    }
}


std::pair<std::vector<int>, std::vector<int>> TeorethicInformation::Bresenham(int x_1, int y_1, int x_2, int y_2)
{
    /*
        To find the set of cells hit by a laser beam
    */
    std::vector<int> x_bres;
    std::vector<int> y_bres;

    int x = x_1;
    int y = y_1;

    int delta_x = abs(x_2 - x_1);
    int delta_y = abs(y_2 - y_1);

    int s_x = getSign(x_1, x_2);
    int s_y = getSign(y_1, y_2);
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
    return std::pair<std::vector<int>, std::vector<int>>{x_bres, y_bres};
}

float TeorethicInformation::calculateEntropy(float probability)
{
    /*
        To calculate the entropy
    */
    return probability * log2(probability);
}


std::vector<float> TeorethicInformation::calculateCellIntersectionPoints(
    std::vector<float> & laser_start, std::vector<float> & laser_end,
    std::vector<float> cell_start, std::vector<float> cell_end)
{
    /*
        Initial point laser beam: laser_start
        Final point laser beam: laser_end
        Initial point cell: cell_start
        Final point cell: cell_end
    */
    float x1 = laser_start[0];
    float x2 = laser_end[0];
    float x3 = cell_start[0];
    float x4 = cell_end[0];

    float y1 = laser_start[1];
    float y2 = laser_end[1];
    float y3 = cell_start[1];
    float y4 = cell_end[1];

    float den = ((x2-x1)*(y4-y3) - (x4-x3)*(y2-y1));
    if (den == 0.0f)
    {
        // Parallel lines or not intersection at all
        return {};
    }
    else
    {
        float x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2-x1)) / den;
        float y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2-y1)) / den;
        return {x, y};
    }
}

void TeorethicInformation::updateCellLimits(
    std::vector<float> & initial_x, std::vector<float> & initial_y, std::vector<float> & final_x, std::vector<float> & final_y,
    float & limit_x, float & limit_y, std::vector<float> & cell_limits, std::vector<int> & robot_grid_pos, std::vector<int> & final_grid_pos)
{
    /*
        To calculate grid grid limits for intersection
    */
    if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] >= robot_grid_pos[1])
    {
        // X greater and Y greater. WRO final points
        final_x[0] = limit_x + m_resolution;
        final_x[2] = limit_x + m_resolution;

        cell_limits[2] = limit_y;
        cell_limits[3] = limit_y + m_resolution;
    }

    if (final_grid_pos[0] >= robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
    {
        // X greater and Y minor. WRO final points
        initial_y[2] = limit_y - m_resolution;
        initial_y[3] = limit_y - m_resolution;

        final_y[1] = limit_y - m_resolution;
        final_y[3] = limit_y - m_resolution;

        cell_limits[2] = limit_y - m_resolution;
        cell_limits[3] = limit_y;
    }

    if (final_grid_pos[0] < robot_grid_pos[0] && final_grid_pos[1] < robot_grid_pos[1])
    {
        // X minor and Y minor. WRO final points
        initial_x[2] = limit_x - m_resolution;
        initial_x[3] = limit_x - m_resolution;
        initial_y[2] = limit_y - m_resolution;
        initial_y[3] = limit_y - m_resolution;

        final_x[0] = limit_x - m_resolution;
        final_x[2] = limit_x - m_resolution;
        final_y[1] = limit_y - m_resolution;
        final_y[3] = limit_y - m_resolution;

        cell_limits[0] = limit_x - m_resolution;
        cell_limits[1] = limit_x;
        cell_limits[2] = limit_y - m_resolution;
        cell_limits[3] = limit_y;
    }
}

float TeorethicInformation::calculateLogs(float probability)
{
    /*
        To calculate the log-odds
    */
    return log(probability / (1 - probability));
}

float TeorethicInformation::calculateMapMutualInformation()
{
    /*
        To calculate map mutual information, this is the summation
        of all cells mutual information
    */
    float sum = 0.0f;
    for (int i = 0; i < m_num_cells; ++i)
    {
        for (int j = 0; j < m_num_cells; ++j)
        {
            sum += m_mutual_grid[i][j];
        }
    }
    return sum;
}

std::unordered_map<TeorethicInformation::map_tuple, float, TeorethicInformation::HashTuple> TeorethicInformation::computeProbabilities(std::vector<std::vector<float>>& meas_outcm)
{
    /*
        To compute all the possible combinations of a grid cell, given a set of measurement outcomes
    */
    // Cleaning measurement outcomes map
    std::unordered_map<map_tuple, float, HashTuple> map_out;
    std::unordered_map<map_tuple, float, HashTuple>::iterator it_out;

    // The number of measurements
    int k = meas_outcm.size(); 
    int r = 1;

    float p_free = meas_outcm[0][2];
    float p_occ = meas_outcm[0][1];
    float p_un = meas_outcm[0][0];

    // Root
    map_out.insert(map_pair(std::make_tuple(0, 0, 0), 1.0f));

    // First measurement
    map_out.insert(map_pair(std::make_tuple(1, 0, 0), p_free));
    map_out.insert(map_pair(std::make_tuple(0, 1, 0), p_occ));
    map_out.insert(map_pair(std::make_tuple(0, 0, 1), p_un));

    for (int i = r; r < k; ++r)
    {
        std::vector<map_tuple> tup_vct;
        std::vector<float> acc_prob;

        for (it_out = map_out.begin(); it_out != map_out.end(); ++it_out)
        {
            // Index
            int fr_idx = std::get<0>(it_out->first);
            int oc_idx = std::get<1>(it_out->first);
            int un_idx = std::get<2>(it_out->first);

            // Measurement outcome probability
            float free_prop = meas_outcm[r][0];
            float occ_prop = meas_outcm[r][1];
            float un_prop = meas_outcm[r][2];

            if (fr_idx + oc_idx + un_idx == r)
            {
                // Searching for the current combination in this level
                std::vector<map_tuple>::iterator it_comb;
                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(fr_idx + 1, oc_idx, un_idx));

                // Free
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += it_out->second * free_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(fr_idx + 1, oc_idx, un_idx));
                    acc_prob.push_back(it_out->second * free_prop);
                }

                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(fr_idx, oc_idx + 1, un_idx));

                // Occupied
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += it_out->second * occ_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(fr_idx, oc_idx + 1, un_idx));
                    acc_prob.push_back(it_out->second * occ_prop);
                }

                it_comb = std::find(tup_vct.begin(), tup_vct.end(), std::make_tuple(fr_idx, oc_idx, un_idx + 1));

                // Unobserved
                if (it_comb != tup_vct.end())
                {
                    acc_prob[it_comb - tup_vct.begin()] += it_out->second * un_prop;
                }
                else
                {
                    tup_vct.push_back(std::make_tuple(fr_idx, oc_idx, un_idx + 1));
                    acc_prob.push_back(it_out->second * un_prop);
                }
            }
        }
        // Inserting the elements into the map
        for (int k = 0; k < tup_vct.size(); ++k)
        {
            map_out.insert(map_pair(tup_vct[k], acc_prob[k]));
        }
    }
    return map_out;
}

float TeorethicInformation::euclideanDistance(float x_1, float y_1, float x_2, float y_2)
{
    /*
        To calculate the euclidean distance between two points
    */
    float diff_x = x_2 - x_1;
    float diff_y = y_2 - y_1;

    return sqrt(diff_x*diff_x + diff_y*diff_y);
}

int TeorethicInformation::getSign(int n_1, int n_2)
{
    /*
        To get the sign of an operation, used for Bresenham algorithm
    */
    int difference = n_2 - n_1;

    if (difference == 0) { return 0; }
    else if (difference < 0) { return -1; }
    else { return 1; }
}

std::vector<int> TeorethicInformation::getGridPosition(float x, float y)
{
    /*
        To maps the current position into grid coordinates
    */
    int x_cell = floor((1 / m_resolution) * x);
    int y_cell = floor((1 / m_resolution) * y);

    return {x_cell, y_cell};
}

std::vector<float> TeorethicInformation::laserHitDistance(std::vector<float> const& robot_pose, float range, float angle)
{
    /*
        To get the distance where the laser beam hits something
            - Applying RBT from the global to the sensor frame
    */
    float x_tf = (range * cos(robot_pose[2] + angle)) + robot_pose[0];
    float y_tf = (range * sin(robot_pose[2] + angle)) + robot_pose[1];

    return {x_tf, y_tf};
}

void TeorethicInformation::initializeGrids()
{
    /*
        To create the grid
    */
    for (int i = 0; i < m_num_cells; ++i)
    {
        // Adding columns
        m_grid[i].resize(m_num_cells);
        m_mutual_grid[i].resize(m_num_cells);
        for (int j = 0; j < m_num_cells; ++j)
        {
            m_grid[i][j] = 0;
            m_mutual_grid[i][j] = 0.0f;
        }
    }
}

float TeorethicInformation::measurementOutcomeEntropy(map_tuple const& meas_outcome)
{
    /*
        To calculate the measurement outcome entropy (Measurement outcome in the form <fr, oc, un>)
            - Calculate Log-Odds from initial probability guess
            - Calculate the probability from those logs
            - Calculate the entropy with the retrieved probability
    */
    float entropy = std::get<0>(meas_outcome) * calculateEntropy(probabilityFromLogs(calculateLogs(0.3f))) + 
                    std::get<1>(meas_outcome) * calculateEntropy(probabilityFromLogs(calculateLogs(0.7f))) + 
                    std::get<2>(meas_outcome) * calculateEntropy(probabilityFromLogs(calculateLogs(0.5f)));
    return entropy;
}

float TeorethicInformation::probabilityFromLogs(float log)
{
    /*
        To transform the Log-odds into probability
    */
    return (exp(log) / (1 + exp(log)));
}

float TeorethicInformation::probabilityFromObservation(float range_1, float range_2)
{
    /*
        To calculate the probability of a cell being observed by a given measurement
    */

    float max_range = 5.0f;
    float lambda = 0.35f;
    float nu = 0.28f;

    range_1 = (range_1 > max_range) ? max_range : range_1;
    range_2 = (range_2 > max_range) ? max_range : range_2;

    return nu * (exp(-lambda*range_1) - exp(-lambda*range_2));
}

std::vector<std::vector<float>> TeorethicInformation::retreiveMeasurementOutcomes(std::vector<int> cell)
{
    /*
        To get all the measurement outcomes for the current cell
    */

    std::vector<std::vector<float>> meas_outcomes;
    std::map<std::vector<int>, std::vector<std::vector<float>>>::iterator it_cells;
    it_cells = m_cell_probabilities.find({cell[0], cell[1]});

    if (it_cells != m_cell_probabilities.end())
    {
        // Exploring the measurement outcomes for the specific cell
        for (int i = 0; i < it_cells->second.size(); ++i)
        {
            // Append the measurement outcomes for the current cell
            meas_outcomes.push_back({it_cells->second[i][0], it_cells->second[i][1], it_cells->second[i][2]});
        }
    }
    return meas_outcomes;
}


void TeorethicInformation::updateCellMutualInformation(float mut_inf, std::vector<int> cell)
{
    /*
        To update the mutual information for each individual cell
        This is the result of the summation of 3.12
    */
    m_mutual_grid[cell[0]][cell[1]] = mut_inf;
}

#include "slam_toolbox/experimental/utils.hpp"
#include <iostream>

namespace utils
{
    namespace grid_operations
    {
        void updateCellLimits(std::vector<kt_double>& initial_x, std::vector<kt_double>& initial_y, std::vector<kt_double>& final_x, 
            std::vector<kt_double>& final_y, kt_double limit_x, kt_double limit_y, std::vector<kt_double>& cell_limits, karto::Vector2<int> const& robot_grid_pos, 
            karto::Vector2<int> const& final_grid_pos, kt_double resolution)
        {
            /*
                To calculate grid grid limits for intersection
            */
            if (final_grid_pos.GetX() < robot_grid_pos.GetX() && final_grid_pos.GetY() >= robot_grid_pos.GetY())
            {
                // X greater and Y greater. WRO final points
                final_x[0] = limit_x + resolution;
                final_x[2] = limit_x + resolution;

                cell_limits[2] = limit_y;
                cell_limits[3] = limit_y + resolution;
            }

            if (final_grid_pos.GetX() >= robot_grid_pos.GetX() && final_grid_pos.GetY() < robot_grid_pos.GetY())
            {
                // X greater and Y minor. WRO final points
                initial_y[2] = limit_y - resolution;
                initial_y[3] = limit_y - resolution;

                final_y[1] = limit_y - resolution;
                final_y[3] = limit_y - resolution;

                cell_limits[2] = limit_y - resolution;
                cell_limits[3] = limit_y;
            }

            if (final_grid_pos.GetX() < robot_grid_pos.GetX() && final_grid_pos.GetY() < robot_grid_pos.GetY())
            {
                // X minor and Y minor. WRO final points
                initial_x[2] = limit_x - resolution;
                initial_x[3] = limit_x - resolution;
                initial_y[2] = limit_y - resolution;
                initial_y[3] = limit_y - resolution;

                final_x[0] = limit_x - resolution;
                final_x[2] = limit_x - resolution;
                final_y[1] = limit_y - resolution;
                final_y[3] = limit_y - resolution;

                cell_limits[0] = limit_x - resolution;
                cell_limits[1] = limit_x;
                cell_limits[2] = limit_y - resolution;
                cell_limits[3] = limit_y;
            }
        }

        int signum(int num)
        {
            /*
                To get the sign of an operation, used by Bresenham algorithm
            */
            if (num < 0) return -1; 
            if (num >= 1) return 1;
            return 0;
        }

        std::pair<std::vector<int>, std::vector<int>> rayCasting(
            karto::Vector2<int> const& initial_pt, karto::Vector2<int> const& final_pt)
        {
            /*
                To find the set of cells hit by a laser beam
                This is based on Bresenham algorithm
            */
            std::vector<int> x_bres;
            std::vector<int> y_bres;

            int x = initial_pt.GetX();
            int y = initial_pt.GetY();

            int delta_x = abs(final_pt.GetX() - initial_pt.GetX());
            int delta_y = abs(final_pt.GetY() - initial_pt.GetY());

            int s_x = signum(final_pt.GetX() - initial_pt.GetX());
            int s_y = signum(final_pt.GetY() - initial_pt.GetY());
            bool interchange = false;

            if (delta_y > delta_x)
            {
                int temp = delta_x;
                delta_x = delta_y;
                delta_y = temp;
                interchange = true;
            }
            else { interchange = false; }

            int a_res = 2 * delta_y;
            int b_res = 2 * (delta_y - delta_x);
            int e_res = (2 * delta_y) - delta_x;

            x_bres.push_back(x);
            y_bres.push_back(y);

            for (int i = 1; i < delta_x; ++i)
            {
                if (e_res < 0)
                {
                    if (interchange) { y += s_y; }
                    else { x += s_x; }
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
            x_bres.push_back(final_pt.GetX());
            y_bres.push_back(final_pt.GetY());

            return std::pair<std::vector<int>, std::vector<int>>{x_bres, y_bres};
        }


        karto::Vector2<int> getGridPosition(karto::Vector2<kt_double> const& pose, kt_double resolution)
        {
            int x_cell = floor((pose.GetX() / resolution));
            int y_cell = floor((pose.GetY() / resolution));

            return karto::Vector2<int>{x_cell, y_cell};
        }

        std::vector<kt_double> calculateCellIntersectionPoints(karto::Vector2<kt_double> const & laser_start, 
            karto::Vector2<kt_double> const & laser_end, std::vector<kt_double> cell_start, std::vector<kt_double> cell_end)
        {
            /*
                Initial point laser beam: laser_start
                Final point laser beam: laser_end
                Initial point cell: cell_start
                Final point cell: cell_end
            */
            kt_double x1 = laser_start.GetX();
            kt_double x2 = laser_end.GetX();
            kt_double x3 = cell_start[0];
            kt_double x4 = cell_end[0];

            kt_double y1 = laser_start.GetY();
            kt_double y2 = laser_end.GetY();
            kt_double y3 = cell_start[1];
            kt_double y4 = cell_end[1];

            kt_double den = ((x2 - x1)*(y4 - y3) - (x4 - x3)*(y2 - y1));
            if (den == 0.0f)
            {
                // Parallel lines or not intersection at all
                return {};
            }
            else
            {
                kt_double x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2 - x1)) / den;
                kt_double y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2 - y1)) / den;
                return {x, y};
            }
        }

        std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
            karto::Vector2<kt_double> const & laser_start, karto::Vector2<kt_double> const & laser_end, 
            karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos,
            kt_double limit_x, kt_double limit_y, kt_double resolution)
        {
            // Cell limits: min_x, max_x, min_y, max_y
            std::vector<kt_double> cell_limits {limit_x, limit_x + resolution, limit_y, limit_y + resolution};

            // Initial points for each of the 4 corners
            std::vector<kt_double> initial_x {limit_x, limit_x, limit_x + resolution, limit_x + resolution};
            std::vector<kt_double> initial_y {limit_y, limit_y, limit_y + resolution, limit_y + resolution};
            
            // Final points for each of the 4 corners
            std::vector<kt_double> final_x {limit_x + resolution, limit_x, limit_x + resolution, limit_x};
            std::vector<kt_double> final_y {limit_y, limit_y + resolution, limit_y, limit_y + resolution};

            // Set the new cell limits
            updateCellLimits(initial_x, initial_y, final_x, final_y, limit_x, limit_y, cell_limits, robot_grid_pos, final_grid_pos, resolution);

            std::vector<kt_double> inter_x, inter_y;

            for (int k = 0; k < 4; ++k)
            {
                std::vector<kt_double> intersection = calculateCellIntersectionPoints(laser_start, laser_end, {initial_x[k], initial_y[k]}, {final_x[k], final_y[k]});
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
            return std::pair<std::vector<kt_double>, std::vector<kt_double>>{inter_x, inter_y}; 
        }

        // template<typename T>
        // void initializeGrid(std::vector<std::vector<T>> & grid, int num_rows, int num_columns)
        // {
        //     /*
        //         To create the grid
        //     */
        //     std::cout << "Grid Initializaion -----------------" << std::endl;
        //     for (int i = 0; i < num_rows; ++i)
        //     {
        //         // Adding columns
        //         grid[i].resize(num_columns);
        //         for (int j = 0; j < num_columns; ++j)
        //         {
        //             grid[i][j] = static_cast<T>(0);
        //         }
        //     }
        // }

        void test()
        {
            std::cout << "Test" << std::endl;
        }
    }// namespace grid_operations
} // namespace utils
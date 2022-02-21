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
            /**
             * To calculate grid limits for intersection
             * Arguments:
                * initial_x [std::vector<kt_double>]: Cell initial limits in x (4 possible limits)
                * initial_y [std::vector<kt_double>]: Cell initial limits in y (4 possible limits)
                * final_x [std::vector<kt_double>]: Cell final limits in x (4 possible limits)
                * final_y [std::vector<kt_double>]: Cell final limits in y (4 possible limits)
                * limit_x [kt_double]: Current cell position in x (Transformed from int to kt_double)
                * limit_y [kt_double]: Current cell position in y (Transformed from int to kt_double)
                * cell_limits [std::vector<kt_double>]: Cell final points for assertion in x and y
                * robot_grip_position [std::vector<kt_double>]: Initial laser beam position
                * final_grid_position [std::vector<kt_double>]: Final laser beam position
                * resolution [kt_double]: Cell resolution
             * Return:
                * Void
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
            /**
             * Get the sign of an operation, used by Bresenham algorithm
             * Arguments:
                * num [int]: Number for perform the sign operation
             * Return:
                * int: Sign
             */
            if (num < 0) return -1;
            if (num >= 1) return 1;
            return 0;
        }

        std::vector<karto::Vector2<int>> rayCasting(
            karto::Vector2<int> const& initial_pt, karto::Vector2<int> const& final_pt)
        {
            /**
             * Find the set of cells hit by a laser beam (Bresenham algorithm)
             * Arguments:
                * initial_pt [karto::Vector2<int>]: Laser beam initial position
                * final_pt [karto::Vector2<int>]: Laser beam final position
             * Return:
                * std::vector<karto::Vector2<int>>: Vector of cells visited by the given laser beam
             */
            std::vector<karto::Vector2<int>> cells;

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

            cells.push_back(karto::Vector2<int>{x, y});

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
                cells.push_back(karto::Vector2<int>{x, y});
            }
            // Delete the current robot cell
            cells.erase(cells.begin());

            // Adding last hit cell to the set
            cells.push_back(karto::Vector2<int>{final_pt.GetX(), final_pt.GetY()});

            return cells;
        }


        karto::Vector2<int> getGridPosition(karto::Vector2<kt_double> const& pose, kt_double resolution)
        {
            /**
             * Mapping a continuous position into a grid position
             * Arguments:
                * pose [karto::Vector2<kt_double>]: Continuos pose
                * resolution [kt_double]: Cell resolution
             * Return:
                * karto::Vector2<int>: Grid position
             */
            int x_cell = floor((pose.GetX() / resolution));
            int y_cell = floor((pose.GetY() / resolution));

            return karto::Vector2<int>{x_cell, y_cell};
        }

        karto::Vector2<kt_double> calculateCellIntersectionPoints(karto::Vector2<kt_double> const & laser_start,
            karto::Vector2<kt_double> const & laser_end, karto::Vector2<kt_double> const & cell_start, karto::Vector2<kt_double> const & cell_end)
        {
            /**
             * Find the intersection point between a cell line and a laser beam
             * Arguments:
                * laser_start [karto::Vector2<kt_double>]: Laser initial point in x and y
                * laser_end [karto::Vector2<kt_double>]: Laser final point in x and y
                * cell_start [karto::Vector2<kt_double>]: Cell initial point in x and y
                * cell_end [karto::Vector2<kt_double>]: Cell final point in x and y
             * Return:
                * karto::Vector2<kt_double>: Intersection point
             */
            kt_double x1 = laser_start.GetX();
            kt_double x2 = laser_end.GetX();
            kt_double x3 = cell_start.GetX();
            kt_double x4 = cell_end.GetX();

            kt_double y1 = laser_start.GetY();
            kt_double y2 = laser_end.GetY();
            kt_double y3 = cell_start.GetY();
            kt_double y4 = cell_end.GetY();

            kt_double den = ((x2 - x1)*(y4 - y3) - (x4 - x3)*(y2 - y1));

            karto::Vector2<kt_double> intersection;
            if (den == 0.0f)
            {
                // Parallel lines
                return {};
            }
            else
            {
                kt_double x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2 - x1)) / den;
                kt_double y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2 - y1)) / den;
                intersection.SetX(x);
                intersection.SetY(y);
            }
            return intersection;
        }

        std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
            karto::Vector2<kt_double> const & laser_start, karto::Vector2<kt_double> const & laser_end,
            karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos,
            kt_double limit_x, kt_double limit_y, kt_double resolution)
        {
            /**
             * Compute intersection between a cell and a laser beam
             * Arguments:
                * laser_start [karto::Vector2<kt_double>]: Laser initial point in x and y
                * laser_end [karto::Vector2<kt_double>]: Laser final point in x and y
                * robot_grid_pos [karto::Vector2<int>]: Initial grid position in x and y
                * final_grid_pos [karto::Vector2<int>]: Final grid position in x and y
                * limit_x [kt_double]: Current cell position in x (Transformed from int to kt_double)
                * limit_y [kt_double]: Current cell position in y (Transformed from int to kt_double)
                * resolution [kt_double]: Cell resolution
             * Return:
                * std::vector<kt_double>: Intersection point
             */

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
                karto::Vector2<kt_double> start{initial_x[k], initial_y[k]};
                karto::Vector2<kt_double> end{final_x[k], final_y[k]};
                karto::Vector2<kt_double> intersection = calculateCellIntersectionPoints(laser_start, laser_end, start, end);
                if(intersection.Length() != 0)
                {
                    if ((fabs(intersection.GetX()) >= (fabs(cell_limits[0]) - 0.001)) &&
                    (fabs(intersection.GetX()) <= (fabs(cell_limits[1]) + 0.001)) &&
                    (fabs(intersection.GetY()) >= (fabs(cell_limits[2]) - 0.001)) &&
                    (fabs(intersection.GetY()) <= (fabs(cell_limits[3]) + 0.001)))
                    {
                        // Two points where the beam cuts the cell
                        //  - A laser beam can cut the cell at least 1 time (Enter)
                        //  - A laser beam can cut the cell at most 2 times (Enter an exit)
                        inter_x.push_back(intersection.GetX());
                        inter_y.push_back(intersection.GetY());
                    }
                }
            }
            return std::pair<std::vector<kt_double>, std::vector<kt_double>>{inter_x, inter_y};
        }

    } // namespace grid_operations
} // namespace utils

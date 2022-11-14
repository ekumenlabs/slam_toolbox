#include "slam_toolbox/experimental/utils.hpp"
#include <iostream>


namespace utils
{
    namespace grid_operations
    {
        void updateCellLimits(
            std::array<karto::Vector2<kt_double>, 4> & initial_points,
            std::array<karto::Vector2<kt_double>, 4> & final_points,
            karto::Vector2<kt_double> const & current_point,
            std::array<kt_double, 4> & cell_limits,
            utils::Segment2<int> const & discretized_segment,
            kt_double const & resolution
        )
        {
            /**
             * To calculate grid limits for intersection
             * Arguments:
                * initial_points [std::array<karto::Vector2<kt_double>, 4>]: Cell initial point in x and y (4 corners)
                * final_points [std::array<karto::Vector2<kt_double>, 4>]: Cell final point in x and y (4 corners)
                * current_point [karto::Vector2<kt_double>]: Current cell position
                * cell_limits [std::array<kt_double, 4>]: Cell final points for assertion in x and y
                * robot_grip_position [std::vector<kt_double>]: Initial laser beam position
                * final_grid_position [std::vector<kt_double>]: Final laser beam position
                * resolution [kt_double]: Cell resolution
             * Return:
                * Void
             */
            if (discretized_segment.end.GetX() < discretized_segment.start.GetX()
                && discretized_segment.end.GetY() >= discretized_segment.start.GetY()
            )
            {
                final_points[0].SetX(current_point.GetX() + resolution);
                final_points[2].SetX(current_point.GetX() + resolution);

                cell_limits[2] = current_point.GetY();
                cell_limits[3] = current_point.GetY() + resolution;
            }

            if (discretized_segment.end.GetX() >= discretized_segment.start.GetX()
                && discretized_segment.end.GetY() < discretized_segment.start.GetY()
            )
            {
                initial_points[2].SetY(current_point.GetY() - resolution);
                initial_points[3].SetY(current_point.GetY() - resolution);

                final_points[1].SetY(current_point.GetY() - resolution);
                final_points[3].SetY(current_point.GetY() - resolution);

                cell_limits[2] = current_point.GetY() - resolution;
                cell_limits[3] = current_point.GetY();
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


        karto::Vector2<int> discretize(karto::Vector2<kt_double> const & position, kt_double const & resolution)
        {
            /**
             * Mapping a continuous position into a grid position
             * Arguments:
                * pose [karto::Vector2<kt_double>]: Continuos pose
                * resolution [kt_double]: Cell resolution
             * Return:
                * karto::Vector2<int>: Grid position
             */
            const int x_cell = floor((position.GetX() / resolution));
            const int y_cell = floor((position.GetY() / resolution));

            return karto::Vector2<int>{x_cell, y_cell};
        }

        std::optional<int> returnint(bool b)
        {
            if (b)
            {
                return 2;
            }
            return {};
        }


        karto::Vector2<kt_double> calculateCellIntersectionPoints(
            utils::Segment2<kt_double> const & segment_1,
            utils::Segment2<kt_double> const & segment_2
        )
        {
            /**
             * Find the intersection point between two segments (cell and laser beam)
             * Arguments:
                * segment_1 [utils::Segment2<kt_double>]: First segment points
                * segment_2 [utils::Segment2<kt_double>]: Second segment points
             * Return:
                * karto::Vector2<kt_double>: Intersection point
             */

            const kt_double x1 = segment_1.start.GetX();
            const kt_double x2 = segment_1.end.GetX();
            const kt_double x3 = segment_2.start.GetX();
            const kt_double x4 = segment_2.end.GetX();

            const kt_double y1 = segment_1.start.GetY();
            const kt_double y2 = segment_1.end.GetY();
            const kt_double y3 = segment_2.start.GetY();
            const kt_double y4 = segment_2.end.GetY();

            const kt_double den = ((x2 - x1)*(y4 - y3) - (x4 - x3)*(y2 - y1));

            if (den != 0.0f)
            {
                const kt_double x = ((x2*y1 - x1*y2)*(x4 - x3) - (x4*y3 - x3*y4)*(x2 - x1)) / den;
                const kt_double y = ((x2*y1 - x1*y2)*(y4 - y3) - (x4*y3 - x3*y4)*(y2 - y1)) / den;
                const karto::Vector2<kt_double> intersection { x, y };
                return intersection;
            }
            return {};
        }


        std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
            utils::Segment2<kt_double> const & segment,
            karto::Vector2<kt_double> const & current_point,
            kt_double const & resolution
        )
        {
            /**
             * Compute intersection between a cell and a laser beam
             * Arguments:
                * laser_start [karto::Vector2<kt_double>]: Laser initial point in x and y
                * laser_end [karto::Vector2<kt_double>]: Laser final point in x and y
                * robot_grid_pos [karto::Vector2<int>]: Initial grid position in x and y
                * final_grid_pos [karto::Vector2<int>]: Final grid position in x and y
                * current_point [kt_double]: Current cell position in x and y
                * resolution [kt_double]: Cell resolution
             * Return:
                * std::vector<kt_double>: Intersection point
             */
            utils::Segment2<int> discretized_segment {
                discretize(segment.start, resolution),
                discretize(segment.end, resolution)
            };

            std::array<kt_double, 4> cell_limits {
                current_point.GetX(),
                current_point.GetX() + resolution,
                current_point.GetY(),
                current_point.GetY() + resolution
            };

            std::array<karto::Vector2<kt_double>, 4> initial_points{{
                {current_point.GetX(), current_point.GetY()},
                {current_point.GetX(), current_point.GetY()},
                {current_point.GetX() + resolution, current_point.GetY() + resolution},
                {current_point.GetX() + resolution, current_point.GetY() + resolution}
            }};

            std::array<karto::Vector2<kt_double>, 4> final_points{{
                {current_point.GetX() + resolution, current_point.GetY()},
                {current_point.GetX(), current_point.GetY() + resolution},
                {current_point.GetX() + resolution, current_point.GetY()},
                {current_point.GetX(), current_point.GetY() + resolution}
            }};

            // Set the new cell limits
            updateCellLimits(initial_points, final_points, current_point, cell_limits, discretized_segment, resolution);

            std::vector<kt_double> inter_x, inter_y;

            /*
                initial {current_point.GetX(), current_point.GetY()}
                final {current_point.GetX() + res, current_point.GetY()}

                initial {current_point.GetX(), current_point.GetY()}
                final {current_point.GetX(), current_point.GetY() + res}

                initial {current_point.GetX() + res, current_point.GetY() + res}
                final {current_point.GetX() + res, current_point.GetY()}

                initial {current_point.GetX() + res, current_point.GetY() + res}
                final {current_point.GetX(), current_point.GetY() + res}
            */

            for (int k = 0; k < 4; ++k)
            {
                karto::Vector2<kt_double> start { initial_points[k].GetX(), initial_points[k].GetY() };
                karto::Vector2<kt_double> end { final_points[k].GetX(), final_points[k].GetY() };

                utils::Segment2<kt_double> cell_segment { start, end };

                karto::Vector2<kt_double> intersection = calculateCellIntersectionPoints(segment, cell_segment);

                if(intersection.Length() != 0)
                {
                    if ((fabs(intersection.GetX()) >= (fabs(cell_limits[0]) - 0.001)) &&
                        (fabs(intersection.GetX()) <= (fabs(cell_limits[1]) + 0.001)) &&
                        (fabs(intersection.GetY()) >= (fabs(cell_limits[2]) - 0.001)) &&
                        (fabs(intersection.GetY()) <= (fabs(cell_limits[3]) + 0.001))
                    )
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

#include "slam_toolbox/experimental/utils.hpp"
#include <iostream>

#include <iomanip>

namespace utils
{
    namespace grid_operations
    {
        void updateCellLimits(
            std::array<karto::Vector2<kt_double>, 4> & initial_points,
            std::array<karto::Vector2<kt_double>, 4> & final_points,
            karto::Vector2<kt_double> const & cell_reference_pos,
            std::array<kt_double, 4> & cell_limits,
            utils::Segment2<int> const & beam_discrete_segment,
            kt_double const & resolution
        )
        {
            /**
             * To calculate grid limits for intersection
             * Arguments:
                * initial_points [std::array<karto::Vector2<kt_double>, 4>]: Cell initial points in x and y (4 corners)
                * final_points [std::array<karto::Vector2<kt_double>, 4>]: Cell final points in x and y (4 corners)
                * cell_reference_pos [karto::Vector2<kt_double>]: Current cell position
                * cell_limits [std::array<kt_double, 4>]: Cell limits in x and y
                * beam_discrete_segment [utils::Segment2<int>]: Beam segment to identify its direction
                * resolution [kt_double]: Cell resolution
             * Return:
                * Void
             */
            if (beam_discrete_segment.end.GetX() < beam_discrete_segment.start.GetX()
                && beam_discrete_segment.end.GetY() >= beam_discrete_segment.start.GetY()
            )
            {
                final_points[0].SetX(cell_reference_pos.GetX() + resolution);
                final_points[2].SetX(cell_reference_pos.GetX() + resolution);

                cell_limits[2] = cell_reference_pos.GetY();
                cell_limits[3] = cell_reference_pos.GetY() + resolution;
            }

            if (beam_discrete_segment.end.GetX() >= beam_discrete_segment.start.GetX()
                && beam_discrete_segment.end.GetY() < beam_discrete_segment.start.GetY()
            )
            {
                initial_points[2].SetY(cell_reference_pos.GetY() - resolution);
                initial_points[3].SetY(cell_reference_pos.GetY() - resolution);

                final_points[1].SetY(cell_reference_pos.GetY() - resolution);
                final_points[3].SetY(cell_reference_pos.GetY() - resolution);

                cell_limits[2] = cell_reference_pos.GetY() - resolution;
                cell_limits[3] = cell_reference_pos.GetY();
            }
        }


        int signum(int num)
        {
            /**
             * Get the sign of an operation, used by Bresenham algorithm
             * Arguments:
                * num [int]: Number to perform the sign operation
             * Return:
                * int: Integer representing the sign
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
                * position [karto::Vector2<kt_double>]: Position to be discretized
                * resolution [kt_double]: Cell resolution
             * Return:
                * karto::Vector2<int>: Grid position
             */
            const int x_cell = floor((position.GetX() / resolution));
            const int y_cell = floor((position.GetY() / resolution));

            return karto::Vector2<int>{ x_cell, y_cell };
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
            karto::Vector2<kt_double> const & cell_reference_pos,
            kt_double const & resolution
        )
        {
            /**
             * Compute intersection between a cell and a laser beam
             * Arguments:
                * segment [karto::Segment2<kt_double>]: Laser beam initial and final points
                * cell_reference_pos [karto::Vector2<kt_double>]: Current cell position
                * resolution [kt_double]: Cell resolution
             * Return:
                * std::pair<std::vector<kt_double>, std::vector<kt_double>>: Intersection points for current beam and cell
             */
            utils::Segment2<int> beam_discrete_segment {
                discretize(segment.start, resolution),
                discretize(segment.end, resolution)
            };

            std::array<kt_double, 4> cell_limits {
                cell_reference_pos.GetX(),
                cell_reference_pos.GetX() + resolution,
                cell_reference_pos.GetY(),
                cell_reference_pos.GetY() + resolution
            };

            std::array<karto::Vector2<kt_double>, 4> initial_points {{
                {cell_reference_pos.GetX(), cell_reference_pos.GetY()},
                {cell_reference_pos.GetX(), cell_reference_pos.GetY()},
                {cell_reference_pos.GetX() + resolution, cell_reference_pos.GetY() + resolution},
                {cell_reference_pos.GetX() + resolution, cell_reference_pos.GetY() + resolution}
            }};

            std::array<karto::Vector2<kt_double>, 4> final_points {{
                {cell_reference_pos.GetX() + resolution, cell_reference_pos.GetY()},
                {cell_reference_pos.GetX(), cell_reference_pos.GetY() + resolution},
                {cell_reference_pos.GetX() + resolution, cell_reference_pos.GetY()},
                {cell_reference_pos.GetX(), cell_reference_pos.GetY() + resolution}
            }};

            // Set the new cell limits
            updateCellLimits(initial_points, final_points, cell_reference_pos, cell_limits, beam_discrete_segment, resolution);

            /*
                initial {cell_reference_pos.GetX(), cell_reference_pos.GetY()}
                final {cell_reference_pos.GetX() + res, cell_reference_pos.GetY()}

                initial {cell_reference_pos.GetX(), cell_reference_pos.GetY()}
                final {cell_reference_pos.GetX(), cell_reference_pos.GetY() + res}

                initial {cell_reference_pos.GetX() + res, cell_reference_pos.GetY() + res}
                final {cell_reference_pos.GetX() + res, cell_reference_pos.GetY()}

                initial {cell_reference_pos.GetX() + res, cell_reference_pos.GetY() + res}
                final {cell_reference_pos.GetX(), cell_reference_pos.GetY() + res}
            */

            std::vector<kt_double> inter_x, inter_y;
            for (int k = 0; k < 4; ++k)
            {
                karto::Vector2<kt_double> start { initial_points[k].GetX(), initial_points[k].GetY() };
                karto::Vector2<kt_double> end { final_points[k].GetX(), final_points[k].GetY() };
                utils::Segment2<kt_double> cell_segment { start, end };

                karto::Vector2<kt_double> intersection = calculateCellIntersectionPoints(segment, cell_segment);
                if(intersection.Length() != 0)
                {
                    if ((fabs(intersection.GetX()) >= (fabs(cell_limits[0]) - cell_limit_eps)) &&
                        (fabs(intersection.GetX()) <= (fabs(cell_limits[1]) + cell_limit_eps)) &&
                        (fabs(intersection.GetY()) >= (fabs(cell_limits[2]) - cell_limit_eps)) &&
                        (fabs(intersection.GetY()) <= (fabs(cell_limits[3]) + cell_limit_eps))
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
            return std::make_pair(inter_x, inter_y);
        }
    } // namespace grid_operations
} // namespace utils

#ifndef SLAM_TOOLBOX__EXPERIMENTAL__UTILS_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__UTILS_HPP_

#include <algorithm>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <unordered_map>
#include "lib/karto_sdk/include/karto_sdk/Karto.h"
#include "Eigen/Core"

namespace utils
{
    namespace grid_operations
    {
        void updateCellLimits(std::vector<kt_double>& initial_x, std::vector<kt_double>& initial_y, std::vector<kt_double>& final_x, 
            std::vector<kt_double>& final_y, kt_double limit_x, kt_double limit_y, std::vector<kt_double>& cell_limits, 
            karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos, kt_double resolution);
        int signum(int num);
        std::vector<karto::Vector2<int>> rayCasting(karto::Vector2<int> const& initial_pt, karto::Vector2<int> const& final_pt);
        karto::Vector2<int> getGridPosition(karto::Vector2<kt_double> const& pose, kt_double resolution);
        std::vector<kt_double> calculateCellIntersectionPoints(karto::Vector2<kt_double> const & laser_start, 
            karto::Vector2<kt_double> const & laser_end, std::vector<kt_double> cell_start, std::vector<kt_double> cell_end);
        std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
            karto::Vector2<kt_double> const & laser_start, karto::Vector2<kt_double> const & laser_end, 
            karto::Vector2<int> const& robot_grid_pos, karto::Vector2<int> const& final_grid_pos,
            kt_double limit_x, kt_double limit_y, kt_double resolution);
        void clearVisitedCells(Eigen::MatrixXd & grid);
        void clearVisitedCells(Eigen::MatrixXi & grid);
    } // namespace grid_operations

    namespace tuple_hash
    {
        struct HashTuple
        {
            std::size_t operator() (std::tuple<int, int, int> const& key) const
            {
                /**
                 * Tuple Hashing
                */
                std::size_t hash = 5381u;
                hash = (hash << 5) + hash + std::get<0>(key);
                hash = (hash << 5) + hash + std::get<1>(key);
                hash = (hash << 5) + hash + std::get<2>(key);
                return hash;
            }
        };
    } // namespace tuple_hash

} // namespace utils

#endif

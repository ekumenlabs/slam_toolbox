#ifndef SLAM_TOOLBOX__EXPERIMENTAL__UTILS_HPP_
#define SLAM_TOOLBOX__EXPERIMENTAL__UTILS_HPP_

#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <unordered_map>
#include "lib/karto_sdk/include/karto_sdk/Karto.h"
#include "Eigen/Core"

#include <optional>

namespace utils
{
    template<typename T>
    struct Segment2 {
        karto::Vector2<T> start;
        karto::Vector2<T> end;
    };

    template<typename T>
    struct Box2 {
        karto::Vector2<T> bl_corner;
        karto::Vector2<T> br_corner;
        karto::Vector2<T> tl_corner;
        karto::Vector2<T> tr_corner;
    };

    namespace grid_operations
    {
        void updateCellLimits(
            std::array<karto::Vector2<kt_double>, 4> & initial_points,
            std::array<karto::Vector2<kt_double>, 4> & final_points,
            karto::Vector2<kt_double> const & current_point,
            std::array<kt_double, 4> & cell_limits,
            utils::Segment2<int> const & discretized_segment,
            kt_double const & resolution
        );

        int signum(int num);

        karto::Vector2<int> discretize(
            karto::Vector2<kt_double> const & position,
            kt_double const & resolution
        );

        karto::Vector2<kt_double> calculateCellIntersectionPoints(
            utils::Segment2<kt_double> const & segment_1,
            utils::Segment2<kt_double> const & segment_2
        );

        std::optional<int> returnint(bool b);

        std::pair<std::vector<kt_double>, std::vector<kt_double>> computeLineBoxIntersection(
            utils::Segment2<kt_double> const & segment,
            karto::Vector2<kt_double> const & current_point,
            kt_double const & resolution
        );
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

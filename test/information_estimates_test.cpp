
#ifndef SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_
#define SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

#include <gtest/gtest.h>
#include "slam_toolbox/experimental/information_estimates.hpp"

/*
    Test comments
    List of possible new tests
        1. resizeGridFromScans
        2. findClosestLaserIndexToCell
        3. adjustBeamReadingDistance

    Requirements for tests:
        1. Laser with complete information
        2. Laser with complete information and able to prodive a laser range finder
        3. Laser with complete information and able to prodive a laser range finder
*/


TEST(UtilsInformationEstimatesTests, SigNumTest)
{
    InformationEstimates information_estimates;

    ASSERT_EQ(utils::grid_operations::signum(10), 1) << "FAIL in positive";
    ASSERT_EQ(utils::grid_operations::signum(-10), -1) << "FAIL in negative";
    ASSERT_EQ(utils::grid_operations::signum(0), 0) << "FAIL in zero";
}

// TEST(UtilsInformationEstimatesTests, InformationContentTest)
// {
//     InformationEstimates inf_estimates;
//     kt_double result = inf_estimates.calculateInformationContent(0.75);
//     ASSERT_DOUBLE_EQ(result, 0.81125);
// }

TEST(UtilsInformationEstimatesTests, ScanMassProbabilityBetweenTest)
{
    kt_double range_1 = 4.25;
    kt_double range_2 = 6.12;

    InformationEstimates inf_estimates;
    kt_double probability = inf_estimates.calculateScanMassProbabilityBetween(range_1, range_2);

    ASSERT_DOUBLE_EQ(probability, 0.07522);
}

TEST(UtilsInformationEstimatesTests, GridPositionsTest)
{
    karto::Vector2<kt_double> position{10.3, 52.7};
    kt_double resolution = 0.5;
    karto::Vector2<int> grid_position = utils::grid_operations::getGridPosition(position, resolution);
    ASSERT_EQ(grid_position[0], 20);
    ASSERT_EQ(grid_position[1], 105);
}

TEST(UtilsInformationEstimatesTests, LineBoxIntersectionTest)
{
    kt_double resolution = 0.5;

    karto::Vector2<kt_double> scan_position{ 18.3, 16.1 };
    karto::Vector2<kt_double> beam_end_point{ 22.2, 12.7 };

    karto::Vector2<int> scan_position_cell = utils::grid_operations::getGridPosition(scan_position, resolution);
    karto::Vector2<int> beam_end_cell = utils::grid_operations::getGridPosition(beam_end_point, resolution);

    karto::Vector2<int> const & cell { 11, 6 };

    // Result could be not intersection or intersection
    std::pair<std::vector<kt_double>, std::vector<kt_double>> intersections =
        utils::grid_operations::computeLineBoxIntersection(
            scan_position,
            beam_end_point,
            scan_position_cell,
            beam_end_cell
            cell.GetX() * resolution,
            cell.GetY() * resolution,
            resolution
        );

    ASSERT_TRUE(intersections.first.size() != 0)

    // Enter (d1) and Exit (d2) distances
    std::vector<kt_double> distances;
    for (int k = 0; k < intersections.first.size(); ++k)
    {
        // From robot position to intersection points
        karto::Vector2<kt_double> intersection{
            intersections.first[k],
            intersections.second[k]
        };
        kt_double distance = scan_position.Distance(intersection);
        distances.push_back(distance);
        std::sort(distances.begin(), distances.end());
    }

    ASSERT_DOUBLE_EQ(distance[0], 21.28235);
    ASSERT_DOUBLE_EQ(distance[1], 13.5);
}

TEST(UtilsInformationEstimatesTests, IntersectionPointsTest)
{
    karto::Vector2<kt_double> int_points;

    // Intersection
    karto::Vector2<kt_double> laser_start_i{1.5, 1.6};
    karto::Vector2<kt_double> laser_end_i{6.1, 8.4};
    karto::Vector2<kt_double> cell_start_i{3.6, 8.2};
    karto::Vector2<kt_double> cell_end_i{4.7, 1.6};
    int_points = utils::grid_operations::calculateCellIntersectionPoints(laser_start_i, laser_end_i, cell_start_i, cell_end_i);

    ASSERT_DOUBLE_EQ(int_points.GetX(), 4.06744) << "FAIL in X coordinate intersection";
    ASSERT_DOUBLE_EQ(int_points.GetY(), 5.39535) << "FAIL in Y coordinate intersection";

    // Parallel lines
    karto::Vector2<kt_double> laser_start_p{1.5, 1.5};
    karto::Vector2<kt_double> laser_end_p{6.1, 6.1};
    karto::Vector2<kt_double> cell_start_p{3.6, 3.6};
    karto::Vector2<kt_double> cell_end_p{8.7, 8.7};
    int_points = utils::grid_operations::calculateCellIntersectionPoints(laser_start_p, laser_end_p, cell_start_p, cell_end_p);

    // I should change here the longitud instead
    // .size() or.Lenght()
    ASSERT_DOUBLE_EQ(int_points.GetX(), 0.0) << "FAIL in X coordinate parallel";
    ASSERT_DOUBLE_EQ(int_points.GetY(), 0.0) << "FAIL in Y coordinate parallel";
}

TEST(InformationEstimatesTests, MutualInformationTest)
{
    InformationEstimates inf_estimates;

    std::unique_ptr<karto::LocalizedRangeScan> s1 = std::make_unique<karto::LocalizedRangeScan>();
    std::unique_ptr<karto::LocalizedRangeScan> s2 = std::make_unique<karto::LocalizedRangeScan>();
    std::unique_ptr<karto::LocalizedRangeScan> s3 = std::make_unique<karto::LocalizedRangeScan>();

    karto::Pose2 p1 = karto::Pose2(3.5, 4.0, 0.0);
    karto::Pose2 p2 = karto::Pose2(3.5, 5.5, 0.0);
    karto::Pose2 p3 = karto::Pose2(5.2, 7.3, 0.0);

    s1->SetCorrectedPose(p1);
    s2->SetCorrectedPose(p2);
    s3->SetCorrectedPose(p3);
    karto::BoundingBox2 bb1, bb2, bb3;
    bb1.SetMinimum(karto::Vector2<kt_double>(2.0, 2.0));
    bb1.SetMaximum(karto::Vector2<kt_double>(5.0, 6.0));
    bb2.SetMinimum(karto::Vector2<kt_double>(2.0, 4.0));
    bb2.SetMaximum(karto::Vector2<kt_double>(5.0, 7.0));
    bb3.SetMinimum(karto::Vector2<kt_double>(2.0, 5.0));
    bb3.SetMaximum(karto::Vector2<kt_double>(5.0, 9.0));
    s1->SetBoundingBox(bb1);
    s2->SetBoundingBox(bb2);
    s3->SetBoundingBox(bb3);
    karto::PointVectorDouble pts1, pts2, pts3;
    pts1.push_back(karto::Vector2<double>(3.0, 5.0));
    pts1.push_back(karto::Vector2<double>(3.0, 3.0));
    pts1.push_back(karto::Vector2<double>(1.5, 7.0));
    pts2.push_back(karto::Vector2<double>(4.0, 5.0));
    pts2.push_back(karto::Vector2<double>(4.0, 2.9));
    pts3.push_back(karto::Vector2<double>(3.0, 5.0));
    pts3.push_back(karto::Vector2<double>(5.0, 6.0));
    pts3.push_back(karto::Vector2<double>(3.5, 7.2));
    s1->SetPointReadings(pts1, true);
    s2->SetPointReadings(pts2, true);
    s3->SetPointReadings(pts3, true);
    bool dirty = false;
    s1->SetIsDirty(dirty);
    s2->SetIsDirty(dirty);
    s3->SetIsDirty(dirty);

    std::vector<karto::LocalizedRangeScan *> range_scan_vct;
    range_scan_vct.push_back(s1.get());
    range_scan_vct.push_back(s2.get());
    range_scan_vct.push_back(s3.get());

    // std::vector<kt_double> mut_inf_vct = inf_estimates.findMutualInfo(range_scan_vct);

    // kt_double prev_mutual_info = 0.0;

    // int idx = 0;
    // int final_idx = 0;

    // for (auto & mutual_info : mut_inf_vct)
    // {
    //     if (mutual_info > prev_mutual_info)
    //     {
    //         prev_mutual_info = mutual_info;
    //         final_idx = idx;
    //     }
    //     ++idx;
    // }
    // // Max value will be the one having the largest mutual information value
    // // Since we are extracting one element at a time:
    //     // So the element with the less mutual information will have a group with the largest result
    //     // In this case index 1 will have the less mutual information
    // EXPECT_EQ(final_idx, 1) << "FAIL in laser index";
    // EXPECT_NE(prev_mutual_info, 0.0) << "FAIL in mutual information equality";
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif // SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

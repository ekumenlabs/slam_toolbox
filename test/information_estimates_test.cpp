
#ifndef SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_
#define SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

#include <gtest/gtest.h>
#include "slam_toolbox/experimental/information_estimates.hpp"


TEST(UtilsInformationEstimatesTests, SigNumTest)
{
    InformationEstimates information_estimates;

    ASSERT_EQ(utils::grid_operations::signum(10), 1) << "FAIL in positive";
    ASSERT_EQ(utils::grid_operations::signum(-10), -1) << "FAIL in negative";
    ASSERT_EQ(utils::grid_operations::signum(0), 0) << "FAIL in zero";
}

TEST(UtilsInformationEstimatesTests, RayCastTest)
{
    InformationEstimates information_estimates;

    // First quadrant
    karto::Vector2<int> initial_pt_1{22,43};
    karto::Vector2<int> final_pt_1{29,38};

    std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(initial_pt_1, final_pt_1);

    ASSERT_EQ(cells[2].GetX(), 25) << "FAIL in X position for cell 3";
    ASSERT_EQ(cells[2].GetY(), 41) << "FAIL in Y position for cell 3";

    ASSERT_EQ(cells[4].GetX(), 27) << "FAIL in X position for cell 5";
    ASSERT_EQ(cells[4].GetY(), 39) << "FAIL in Y position for cell 5";

    // Second quadrant
    karto::Vector2<int> initial_pt_2{22,43};
    karto::Vector2<int> final_pt_2{15,38};

    karto::Vector2<int> initial_pt_3{22,43};
    karto::Vector2<int> final_pt_3{29,48};

    karto::Vector2<int> initial_pt_4{22,43};
    karto::Vector2<int> final_pt_4{15,48};

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

    ASSERT_FLOAT_EQ(int_points.GetX(), 4.06744) << "FAIL in X coordinate intersection";
    ASSERT_FLOAT_EQ(int_points.GetY(), 5.39535) << "FAIL in Y coordinate intersection";

    // Parallel lines
    karto::Vector2<kt_double> laser_start_p{1.5, 1.5};
    karto::Vector2<kt_double> laser_end_p{6.1, 6.1};
    karto::Vector2<kt_double> cell_start_p{3.6, 3.6};
    karto::Vector2<kt_double> cell_end_p{8.7, 8.7};
    int_points = utils::grid_operations::calculateCellIntersectionPoints(laser_start_p, laser_end_p, cell_start_p, cell_end_p);

    ASSERT_FLOAT_EQ(int_points.GetX(), 0.0) << "FAIL in X coordinate parallel";
    ASSERT_FLOAT_EQ(int_points.GetY(), 0.0) << "FAIL in Y coordinate parallel";
}

TEST(InformationEstimatesTests, MutualInformationTest)
{
    InformationEstimates inf_estimates;

    std::unique_ptr<karto::LocalizedRangeScan> s1 = std::make_unique<karto::LocalizedRangeScan>();
    std::unique_ptr<karto::LocalizedRangeScan> s2 = std::make_unique<karto::LocalizedRangeScan>();

    karto::Pose2 p1 = karto::Pose2(3.5, 4.0, 0.0);
    karto::Pose2 p2 = karto::Pose2(3.5, 5.5, 0.0);

    s1->SetCorrectedPose(p1);
    s2->SetCorrectedPose(p2);
    karto::BoundingBox2 bb1, bb2;
    bb1.SetMinimum(karto::Vector2<kt_double>(2.0, 2.0));
    bb1.SetMaximum(karto::Vector2<kt_double>(5.0, 6.0));
    bb2.SetMinimum(karto::Vector2<kt_double>(2.0, 4.0));
    bb2.SetMaximum(karto::Vector2<kt_double>(5.0, 7.0));
    s1->SetBoundingBox(bb1);
    s2->SetBoundingBox(bb2);
    karto::PointVectorDouble pts1, pts2;
    pts1.push_back(karto::Vector2<double>(3.0, 5.0));
    pts1.push_back(karto::Vector2<double>(3.0, 3.0));
    pts1.push_back(karto::Vector2<double>(1.5, 7.0));
    pts2.push_back(karto::Vector2<double>(4.0, 5.0));
    pts2.push_back(karto::Vector2<double>(4.0, 2.9));
    s1->SetPointReadings(pts1, true);
    s2->SetPointReadings(pts2, true);
    bool dirty = false;
    s1->SetIsDirty(dirty);
    s2->SetIsDirty(dirty);

    std::vector<karto::LocalizedRangeScan *> range_scan_vct;
    range_scan_vct.push_back(s1.get());
    range_scan_vct.push_back(s2.get());

    std::tuple<int, kt_double> min_inf = inf_estimates.calculateMutualInformation(range_scan_vct);

    int idx;
    kt_double mut_inf;
    std::tie(idx, mut_inf) = min_inf;

    EXPECT_EQ(idx, 1) << "FAIL in laser index";
    EXPECT_NE(mut_inf, 0.0) << "FAIL in mutual information equality";
}

int main(int argc, char ** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

#endif // SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

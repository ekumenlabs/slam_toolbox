
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

    karto::Vector2<int> initial_pt{22,43};
    karto::Vector2<int> final_pt{29,38};

    std::vector<karto::Vector2<int>> cells = utils::grid_operations::rayCasting(initial_pt, final_pt);

    ASSERT_EQ(cells[2].GetX(), 25) << "FAIL in X position for cell 3";
    ASSERT_EQ(cells[2].GetY(), 41) << "FAIL in Y position for cell 3";

    ASSERT_EQ(cells[4].GetX(), 27) << "FAIL in X position for cell 5";
    ASSERT_EQ(cells[4].GetY(), 39) << "FAIL in Y position for cell 5";
}

TEST(UtilsInformationEstimatesTests, IntersectionPointsTest)
{
    karto::Vector2<kt_double> laser_start{1.5, 1.6};
    karto::Vector2<kt_double> laser_end{6.1, 8.4};
    std::vector<kt_double> cell_start{3.6, 8.2};
    std::vector<kt_double> cell_end{4.7, 1.6};

    std::vector<kt_double> int_points = utils::grid_operations::calculateCellIntersectionPoints(laser_start, laser_end, cell_start, cell_end);
    
    ASSERT_FLOAT_EQ(int_points[0], 4.06744) << "FAIL in X coordinate";
    ASSERT_FLOAT_EQ(int_points[1], 5.39535) << "FAIL in Y coordinate";
}

TEST(InformationEstimatesTests, MutualInformationTest)
{
    InformationEstimates inf_estimates;

    karto::LocalizedRangeScan * scan_1 = new karto::LocalizedRangeScan();
    karto::LocalizedRangeScan * scan_2 = new karto::LocalizedRangeScan();

    karto::Pose2 pose_1 = karto::Pose2(3.5, 4.0, 0.0);
    karto::Pose2 pose_2 = karto::Pose2(3.5, 5.5, 0.0);

    scan_1->SetBarycenterPose(pose_1);
    scan_2->SetBarycenterPose(pose_2);

    karto::BoundingBox2 bb1, bb2;
    bb1.SetMinimum(karto::Vector2<kt_double>(2.0, 2.0));
    bb1.SetMaximum(karto::Vector2<kt_double>(5.0, 6.0));
    bb2.SetMinimum(karto::Vector2<kt_double>(2.0, 4.0));
    bb2.SetMaximum(karto::Vector2<kt_double>(5.0, 7.0));
    
    scan_1->SetBoundingBox(bb1);
    scan_2->SetBoundingBox(bb2);
    karto::PointVectorDouble points;
    points.push_back(karto::Vector2<double>(3.0, 5.0));
    points.push_back(karto::Vector2<double>(3.0, 3.1));
    
    scan_1->SetPointReadings(points, true);

    karto::PointVectorDouble laser_readings = scan_1->GetPointReadings(true);

    EXPECT_FLOAT_EQ(laser_readings[0].GetX(), 3.0) << "EXPECT FAIL in X coordinate";
    EXPECT_FLOAT_EQ(laser_readings[0].GetY(), 5.1) << "EXPECT FAIL in Y coordinate";

    ASSERT_FLOAT_EQ(laser_readings[0].GetX(), 3.0) << "FAIL in X coordinate";
    ASSERT_FLOAT_EQ(laser_readings[0].GetY(), 5.0) << "FAIL in Y coordinate";

    // bool dirty = false;
    // s1->SetIsDirty(dirty);
    // s2->SetIsDirty(dirty);
}

int main(int argc, char ** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

// colcon build --packages-select slam_toolbox  --cmake-args -DBUILD_TESTING=ON
// colcon test --packages-select slam_toolbox --event-handlers console_direct+


#endif // SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

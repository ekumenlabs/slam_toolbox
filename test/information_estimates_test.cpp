
#ifndef SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_
#define SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

#include <gtest/gtest.h>
#include "slam_toolbox/experimental/information_estimates.hpp"


TEST(InformationEstimatesTests, BareTest)
{
    InformationEstimates information_estimates;

    ASSERT_EQ(utils::grid_operations::signum(10), 0) << "FAIL in positive";
    ASSERT_EQ(utils::grid_operations::signum(-10), -1) << "FAIL in negative";
    ASSERT_EQ(utils::grid_operations::signum(0), 0) << "FAIL in zero";
}

int main(int argc, char ** argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


#endif // SLAM_TOOLBOX_SLAM_TOOLBOX_INFORMATION_ESTIMATES_TEST_H_

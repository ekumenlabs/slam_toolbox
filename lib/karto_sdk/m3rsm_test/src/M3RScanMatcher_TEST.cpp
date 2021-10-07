#include <gtest/gtest.h>
#include <karto_sdk/M3RScanMatcher.h>
#include <initializer_list>

namespace test {
namespace {

typedef std::initializer_list<uint8_t> MatrixRow;

GTEST_TEST(MockLookupTableTest, MockLookupTableTest_IsMaxResolutionWhenCreatedWithBaseScans) {
  std::vector<std::vector<uint8_t>> data;
  data.reserve(16);
  for (int i = 0; i < 16; ++i) {
    data.push_back(std::vector<uint8_t>(16, 0));
  }
  karto::MockLookupTable lookup_table(data, 1, 2);
  
  EXPECT_TRUE(lookup_table.isMaxResolution());
}

GTEST_TEST(MockLookupTableTest, MockLookupTableTest_IsCreatedFromMockLookupTableOk) {
  std::vector<std::vector<uint8_t>> data;
  int base_table_size = 16;
  data.reserve(base_table_size);
  for (int i = 0; i < base_table_size; ++i) {
    data.push_back(std::vector<uint8_t>(base_table_size, 0));
  }
  karto::MockLookupTable base_lookup_table(data, 2, 3);
  karto::MockLookupTable lookup_table = karto::MockLookupTable::GetLowerResTable(base_lookup_table);
  
  EXPECT_EQ(lookup_table.data_.size(), base_table_size / 2);
  EXPECT_EQ(lookup_table.decimation_factor_, base_lookup_table.decimation_factor_);
  EXPECT_EQ(lookup_table.kernel_width_, base_lookup_table.kernel_width_);
}

GTEST_TEST(MockLookupTableTest, fourByFourTableReduction) {
  std::vector<std::vector<uint8_t>> data;
  int base_table_size = 4;
  data.emplace_back(MatrixRow{0, 1,  5, 4});
  data.emplace_back(MatrixRow{1, 2,  3, 4});
  data.emplace_back(MatrixRow{0,  0, 0, 0});
  data.emplace_back(MatrixRow{10, 0, 3, 4});

  std::vector<std::vector<uint8_t>> expected_mid_res_data;
  expected_mid_res_data.emplace_back(MatrixRow{5, 5});
  expected_mid_res_data.emplace_back(MatrixRow{10, 4});

  std::vector<std::vector<uint8_t>> expected_lowest_res_data;
  expected_lowest_res_data.emplace_back(MatrixRow{10});

  karto::MockLookupTable highest_res_table(data, 2, 3);
  karto::MockLookupTable mid_res_table = karto::MockLookupTable::GetLowerResTable(highest_res_table);
  karto::MockLookupTable lowest_res_table = karto::MockLookupTable::GetLowerResTable(mid_res_table);

  EXPECT_EQ(mid_res_table.data_, expected_mid_res_data);  
  EXPECT_EQ(lowest_res_table.data_, expected_lowest_res_data);
}

GTEST_TEST(MockLookupTableTest, ComputeCostTest) {
  std::vector<std::vector<uint8_t>> data;
  int base_table_size = 4;
  data.emplace_back(MatrixRow{0, 1,  5, 4});
  data.emplace_back(MatrixRow{1, 2,  3, 4});
  data.emplace_back(MatrixRow{0,  0, 0, 0});
  data.emplace_back(MatrixRow{10, 0, 3, 4});
  karto::MockLookupTable table(data, 2, 3);

  std::vector<std::pair<int, int>> occupied_cells;
  occupied_cells.emplace_back(std::make_pair(0, 0));
  occupied_cells.emplace_back(std::make_pair(1, 0));
  occupied_cells.emplace_back(std::make_pair(2, 0));
  occupied_cells.emplace_back(std::make_pair(3, 0));

  EXPECT_EQ(table.computeCost(occupied_cells), 11);
}

GTEST_TEST(LookupTableManagerTest, LookupTableManagerIsCreatedOk) {
  std::vector<std::vector<uint8_t>> data;
  int base_table_size = 4;
  data.emplace_back(MatrixRow{0, 1,  5, 4});
  data.emplace_back(MatrixRow{1, 2,  3, 4});
  data.emplace_back(MatrixRow{0,  0, 0, 0});
  data.emplace_back(MatrixRow{10, 0, 3, 4});

  std::vector<std::vector<uint8_t>> expected_mid_res_data;
  expected_mid_res_data.emplace_back(MatrixRow{5, 5});
  expected_mid_res_data.emplace_back(MatrixRow{10, 4});

  std::vector<std::vector<uint8_t>> expected_lowest_res_data;
  expected_lowest_res_data.emplace_back(MatrixRow{10});

  karto::LookupTableManager manager(data);

  EXPECT_EQ(manager.getAmountOfResolutions(), 3);  
  EXPECT_EQ(manager.getTable(0).data_, data);
  EXPECT_EQ(manager.getTable(1).data_, expected_mid_res_data);
  EXPECT_EQ(manager.getTable(2).data_, expected_lowest_res_data);
}

}  // namespace
}  // namespace test

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

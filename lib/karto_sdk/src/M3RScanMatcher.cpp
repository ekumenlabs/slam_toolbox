#include "karto_sdk/Types.h"
#include <math.h>
#include <assert.h>
#include <boost/serialization/vector.hpp>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <queue>
#include <set>
#include <list>
#include <iterator>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>

#include "karto_sdk/M3RScanMatcher.h"


namespace {
  constexpr int kInitialLookupTableSize{1024};
  typedef std::vector<std::vector<uint8_t>> DataMatrix;

  uint8_t getMax(int row_init, int row_end, int col_init, int col_end, const std::vector<std::vector<uint8_t>>& matrix) {
    int16_t max = -1;
    for (int i = row_init; i <= row_end; i++) {
      for (int j = col_init; j <= col_end; j++) {
        if (max < matrix[i][j]) {
          max = matrix[i][j];
        }
      }
    }
    return max;
  }
}

namespace karto
{

/********************************* MockLookupTable *********************************************/

MockLookupTable::MockLookupTable(const std::vector<std::vector<uint8_t>>& data, int decimation_factor, int kernel_width) {
  is_max_res_ = true;
  decimation_factor_ = decimation_factor;
  kernel_width_ = kernel_width;
  data_ = data;
}

MockLookupTable::MockLookupTable(const MockLookupTable& table) {
  is_max_res_ = table.is_max_res_;
  decimation_factor_ = table.decimation_factor_;
  kernel_width_ = table.kernel_width_;
  data_ = table.data_;
}

MockLookupTable MockLookupTable::GetLowerResTable(const MockLookupTable& higher_res_table) {
  int decimation_factor = higher_res_table.decimation_factor_;
  int kernel_width = higher_res_table.kernel_width_;
  DataMatrix data;
  int new_table_size = higher_res_table.data_.size() / decimation_factor;
  data.reserve(new_table_size);
  for (int i = 0; i < new_table_size; i++) {
    data.push_back(std::vector<uint8_t>(new_table_size, 0));
  }

  for (int i = 0; i < new_table_size; i++) {
    for (int j = 0; j < new_table_size; j++) {
      int row_end = i * decimation_factor + kernel_width - 1;
      int col_end = j * decimation_factor + kernel_width - 1;
      if (row_end > higher_res_table.data_.size() - 1) {
        row_end = higher_res_table.data_.size() - 1;
      }
      if (col_end > higher_res_table.data_.size() - 1) {
        col_end = higher_res_table.data_.size() - 1;
      }
      data[i][j] = getMax(i*decimation_factor, row_end, j*decimation_factor, col_end, higher_res_table.data_);
    }
  }
  MockLookupTable new_table(data, decimation_factor, kernel_width);
  new_table.is_max_res_ = false;
  return new_table;
}

double MockLookupTable::computeCost(const std::vector<std::pair<int, int>>& occupied_cells) {
  // Rasterize scan into table and compute cost
  double cost = 0;
  for (int i = 0; i < occupied_cells.size(); i++) {
    std::pair<int, int> pair = occupied_cells[i];
    cost += data_[pair.first][pair.second];
  }
  return cost;
}

bool MockLookupTable::isMaxResolution() {
  return is_max_res_;
}

int MockLookupTable::getTableSize() {
  return data_.size();
}

/********************************* LookupTableManager **************************************/
LookupTableManager::LookupTableManager(const DataMatrix& full_res_matrix) {
  int decimation_factor = 2;
  int kernel_width = 3;

  // Pre-reserve to avoid reallocating data.
  n_resolutions_ = std::ceil(log2(full_res_matrix.size())) + 1;
  tables_.reserve(n_resolutions_);

  MockLookupTable table(full_res_matrix, decimation_factor, kernel_width);
  tables_.emplace_back(table);
  for (int resolution = 0; resolution < n_resolutions_; resolution++) {
    table = MockLookupTable::GetLowerResTable(table);
    tables_.emplace_back(table);
  }
}

MockLookupTable LookupTableManager::getTable(int resolution) {
  if (resolution >= n_resolutions_ || resolution < 0) {
    throw std::runtime_error("Invalid resolution: " + std::to_string(resolution));
  }
  return tables_[resolution];
}

int LookupTableManager::getAmountOfResolutions() {
  return n_resolutions_;
}

} // namespace karto
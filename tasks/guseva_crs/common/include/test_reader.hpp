#pragma once

#include <cstddef>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include "common.hpp"

namespace guseva_crs {

inline CRS ReadCRSFromFile(std::ifstream &file) {
  std::size_t nz = 0;
  std::size_t nrows = 0;
  std::size_t ncols = 0;
  std::vector<double> values;
  std::vector<std::size_t> cols;
  std::vector<std::size_t> row_ptrs;
  file >> nz >> nrows >> ncols;
  double tmp = {};
  for (std::size_t i = 0; i < nz; i++) {
    file >> tmp;
    values.push_back(tmp);
  }
  std::size_t second_tmp = 0;
  for (std::size_t i = 0; i < nz; i++) {
    file >> second_tmp;
    cols.push_back(second_tmp);
  }

  for (std::size_t i = 0; i < nrows + 1; i++) {
    file >> second_tmp;
    row_ptrs.push_back(second_tmp);
  }
  return CRS(nz, nrows, ncols, values, cols, row_ptrs);
}

inline std::tuple<CRS, CRS, CRS> ReadTestFromFile(const std::string &filename) {
  std::ifstream file(filename);
  return {ReadCRSFromFile(file), ReadCRSFromFile(file), ReadCRSFromFile(file)};
}

}  // namespace guseva_crs

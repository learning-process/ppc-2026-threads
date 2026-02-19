// #pragma once

// #include <climits>
// #include <cmath>
// #include <iostream>
// #include <utility>
// #include <vector>

// namespace guseva_crs {

// constexpr double kZERO = 10e-9;

// template <typename T>
// class CRS {
//  public:
//   size_t nz{};
//   size_t nrows{};
//   size_t ncols{};
//   std::vector<T> values;
//   std::vector<size_t> cols;
//   std::vector<size_t> row_ptrs;

//   // CRS() = default;

//   // CRS(size_t nz_p, size_t nrows_p, size_t ncols_p, std::vector<T> values_p, std::vector<size_t> cols_p,
//   //     std::vector<size_t> row_ptrs_p)
//   //     : nz(nz_p),
//   //       nrows(nrows_p),
//   //       ncols(ncols_p),
//   //       values(std::move(values_p)),
//   //       cols(std::move(cols_p)),
//   //       row_ptrs(std::move(row_ptrs_p)) {}

//   // CRS(const CRS<T> &other)
//   //     : nz(other.nz),
//   //       nrows(other.nrows),
//   //       ncols(other.ncols),
//   //       values(other.values),
//   //       cols(other.cols),
//   //       row_ptrs(other.row_ptrs) {}

//   // CRS(CRS &&other) noexcept
//   //     : nz(std::exchange(other.nz, 0)),
//   //       nrows(std::exchange(other.nrows, 0)),
//   //       ncols(std::exchange(other.ncols, 0)),
//   //       values(std::move(other.values)),
//   //       cols(std::move(other.cols)),
//   //       row_ptrs(std::move(other.row_ptrs)) {
//   //   other.nz = 0;
//   //   other.nrows = 0;
//   //   other.ncols = 0;
//   // }

//   // CRS<T> &operator=(const CRS &other) {
//   //   if (this != &other) {
//   //     nz = other.nz;
//   //     nrows = other.nrows;
//   //     ncols = other.ncols;
//   //     values = other.values;
//   //     cols = other.cols;
//   //     row_ptrs = other.row_ptrs;
//   //   }
//   //   return *this;
//   // }

//   // CRS<T> &operator=(CRS &&other) noexcept {
//   //   if (this != &other) {
//   //     nz = std::exchange(other.nz, 0);
//   //     nrows = std::exchange(other.nrows, 0);
//   //     ncols = std::exchange(other.ncols, 0);
//   //     values = std::move(other.values);
//   //     cols = std::move(other.cols);
//   //     row_ptrs = std::move(other.row_ptrs);
//   //   }
//   //   return *this;
//   // }

//   bool operator==(const CRS<T> &other) const {
//     return nz == other.nz && ncols == other.ncols && nrows == other.nrows &&
//            std::ranges::equal(values, other.values, [](T a, T b) { return std::fabs(a - b) > kZERO; });
//   }

//   void Print() const {
//     std::cout << "Nonzero: " << nz << " Dim: " << nrows << "x" << ncols << "\nValues: [ ";
//     for (const auto &i : values) {
//       std::cout << i << " ";
//     }
//     std::cout << "]\nCols: [";
//     for (const auto &i : cols) {
//       std::cout << i << " ";
//     }
//     std::cout << "]\nRows Pointers: [";
//     for (const auto &i : row_ptrs) {
//       std::cout << i << " ";
//     }
//     std::cout << "]\n";
//   }

//   void PrettyPrint() const {
//     std::vector<std::vector<T>> dense(nrows, std::vector<T>(ncols, T{0}));
//     for (size_t i = 0; i < nrows; i++) {
//       auto start = row_ptrs[i];
//       auto end = row_ptrs[i + 1];
//       for (auto j = start; j < end; j++) {
//         dense[i][cols[j]] = values[j];
//       }
//     }
//     for (const auto &vec : dense) {
//       for (const auto &elem : vec) {
//         std::cout << elem << " ";
//       }
//       std::cout << "\n";
//     }
//   }
// };

// }  // namespace guseva_crs
